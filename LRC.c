/*
# RegeneratingCodes/LRC.c

Copyright (c) 2014, AT&T Intellectual Property.  All other rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software must display the following acknowledgement:  This product includes software developed by the AT&T.
4. Neither the name of AT&T nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AT&T INTELLECTUAL PROPERTY ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL AT&T INTELLECTUAL PROPERTY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Chao Tian
# AT&T Labs-Research
# Bedminster, NJ 07943
# tian@research.att.com

# $Revision: 0.1 $
# $Date: 2014/02/25 $
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "jerasure.h"
#include "reed_sol.h"

#include "regenerating_codes.h"
#include "cauchy.h"


/*
The code is based on the paper "Locally repairable codes"
D. Papailiopoulos and A. Dimakis, ISIT 2012

The coding has two steps: the first step is a conventional Reed-Solomon code, and the second step involves
a coding on some partition of the coded symbols, after which they are carefully placed into each device.

Restriction: (f+1) has to be a multiple of n. The alphabet size 2^w>= n
*/


int encode_LRC(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info)
{
	int i,j,c1;	
	int n = info->req.n;
	int k = info->req.k;
	int f = info->req.f;	
	int w = info->req.w;	
	int subpacket_size = input_size/(k*f);
	int num_of_groups = n/(f+1);
	int base;
	char **data_ptrs = malloc(sizeof(char*)*f);   
	if(data_ptrs==NULL){
		printf("Out of memory.\n");
		return(-1);
	}	     

	for(i=0;i<k;i++)
		memcpy(output[i],input+i*f*subpacket_size, f*subpacket_size);
	jerasure_schedule_encode(k, n-k, w, info->schedule, output, (output+k), subpacket_size*f, ALIGNMENT);
	for(i=0;i<num_of_groups;i++){
		base = i*(f+1);
		for(j=0;j<f+1;j++){
			for(c1=0;c1<f;c1++)
				data_ptrs[c1] = output[base+(j+c1)%(f+1)]+c1*subpacket_size;
			jerasure_do_parity(f, data_ptrs, (output[base+(j+f)%(f+1)]+f*subpacket_size), subpacket_size);
		}
	}
	
	// clean-up
	free(data_ptrs);       
	return(1);
}

int decode_LRC(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info)
{
	int i, j, c1;
	int n = info->req.n;
	int k = info->req.k;
	int w = info->req.w;
	int f = info->req.f;
	int base;
	int num_of_groups = n/(f+1);
	int target_device;
	int subpacket_size = input_size/(f+1);	
	int* erased = jerasure_erasures_to_erased(k,n-k,erasures);
       
	jerasure_schedule_decode_lazy(k,n-k,w,info->bitmatrix,erasures,input,(input+k),subpacket_size*f,ALIGNMENT,0);
	
	for(i=0;i<k;i++)
		memcpy(output+i*f*subpacket_size,input[i], f*subpacket_size);

	// allocate memory for data arrangement
	char **data_ptrs = malloc(sizeof(char*)*f); 	
	if(data_ptrs==NULL){
		printf("Out of memory.\n");
		return(-1);
	}

	for(i=0;i<num_of_groups;i++){
		base = i*(f+1);
		for(j=0;j<f+1;j++){
			target_device = base+(j+f)%(f+1);
			if(erased[target_device]==1){
				for(c1=0;c1<f;c1++)
					data_ptrs[c1] = input[base+(j+c1)%(f+1)]+c1*subpacket_size;
				jerasure_do_parity(f, data_ptrs, (input[base+(j+f)%(f+1)]+f*subpacket_size), subpacket_size);
			}
		}
	}
	// clean-up
	free(data_ptrs);
	free(erased);
	return(1);

}

int repair_encode_LRC(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info)
{
	memcpy(output,input,input_size);		
	return(1);
}

int repair_decode_LRC(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info)
{
	int i,j,counter;	
	int f = info->req.f;	
	int subpacket_size = input_size/(f+1);
	char **data_ptrs = malloc(sizeof(void*)*f);   
	int *helpers_inv_ID = malloc(sizeof(int)*f);
	int base = to_device_ID/(f+1)*(f+1);
	if(data_ptrs==NULL||helpers_inv_ID==NULL){
		printf("Out of memory.\n");
		return(-1);
	}	     

	for(i=0,counter;i<f;i++){
		if(helpers[i]>=base&&helpers[i]<=base+f){
			helpers_inv_ID[(helpers[i]-to_device_ID+f)%(f+1)] = i;		
			counter ++;
		}		
	}
	if(counter<f){
		printf("Insufficient number of helpers\n");
		return(-1);
	}
	for(i=0;i<f+1;i++){
		for(j=1;j<f+1;j++)
			data_ptrs[j-1] = input[helpers_inv_ID[j-1]]+((i+j)%(f+1))*subpacket_size;		
		jerasure_do_parity(f, data_ptrs, (output+i*subpacket_size), subpacket_size);
	}
	
	// clean-up
	free(data_ptrs); 
	free(helpers_inv_ID);      
	return(1);
}
