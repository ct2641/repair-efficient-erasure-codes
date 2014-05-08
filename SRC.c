/*
# RegeneratingCodes/SRC.c

Copyright (c) 2014, AT&T Intellectual Property.  All other rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software must display the following acknowledgement:  This product includes software developed by the AT&T.
4. Neither the name of AT&T nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AT&T INTELLECTUAL PROPERTY ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL AT&T INTELLECTUAL PROPERTY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# 
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
The code is based on the paper "Simple Regenerating Codes: Network Coding for Cloud Storage"
by Dimitris S. Papailiopoulos, Jianqiang Luo, Alexandros G. Dimakis, Cheng Huang, Jin Li 

disk indices:  1      2      3    4  ....  n
symbols:     | X1     X2     X3   X4 ....  Xn    |
             | Yn     Y1     Y2   Y3 ....  Y(n-1)|
             | Z(n-1) Zn     Z1      ....        |     
             | S(n-2) S(n-1) Sn   S1 ....        |

Here f=3 and S is the summation of a parity group (Xi,Yi,Zi). The value of f is a parameter of choice, but f<=n.

However, we will group the information symbols in a different way: along (X1,Yn,Z(n-1))->(X2,Y1,Zn)->... is the file stream.
This arrangement implies that on each device, the systematic file content is in a single contiguous block. Another benefit is that
the decoding matrix/bitmatrix is now fixed, and we can treat data block such as (X1,Yn,Z(n-1)) as a single piece of data, and need 
to call the decoding routine just once.

Restriction: alphabet size 2^w>=n.
*/

int encode_SRC(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info)
{
	int i,j;	
	int n = info->req.n;
	int k = info->req.k;
	int f = info->req.f;
	int subpacket_size = input_size/(k*f);	
  
	// rearrange the memory pointers in preparation for encoding
	char **data_plus_coding_ptrs = malloc(sizeof(char*)*n); 
	if(data_plus_coding_ptrs==NULL){
		printf("Out of memory.\n");
		return(-1);
	}	
	for(i=0;i<k;i++) // i-th device
		memcpy(output[i],input+i*f*subpacket_size,subpacket_size*f);			

	// call jerasure routine for encoding;
	jerasure_schedule_encode(k, n-k, info->req.w, 
			info->schedule, output, output+k, 
			subpacket_size*f, ALIGNMENT);

	for(i=0;i<n;i++){
		for(j=0; j<f;j++)// reuse these pointers for the correct data blocks before XORs vertically/diagonally			
			data_plus_coding_ptrs[j] = output[(i+j)%n]+j*subpacket_size;
		jerasure_do_parity(f, data_plus_coding_ptrs, output[(i+f)%n] + f*subpacket_size, subpacket_size);
	}
	
	// clean up
	free(data_plus_coding_ptrs);       
	return(1);
}

int decode_SRC(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info)
{
	int i, j;
	int n = info->req.n;
	int k = info->req.k;
	int f = info->req.f;
	int subpacket_size = output_size/k;	       
       
	// allocate memory for data arrangement
	char** data_plus_coding_ptrs = malloc(sizeof(char*)*n);
	if(data_plus_coding_ptrs==NULL){
		printf("Out of memory.\n");
		return(-1);
	}
	
	jerasure_schedule_decode_lazy(k, n-k, info->req.w, 
			info->bitmatrix, erasures, 
			input, input+k, 	
		        subpacket_size, ALIGNMENT, 0);

	for(i=0;i<k;i++) // i-th device
		memcpy(output+i*subpacket_size,input[i],subpacket_size);	

	
	int* erased = jerasure_erasures_to_erased(n, n-k, erasures);
	subpacket_size = subpacket_size/f;
	for(i=0;i<n;i++){
		if(erased[(i+f)%n]!=1)
			continue;
		for(j=0; j<f;j++)// reuse these pointers for the correct data blocks before XORs vertically/diagonally			
			data_plus_coding_ptrs[j] = input[(i+j)%n]+j*subpacket_size;
		jerasure_do_parity(f, data_plus_coding_ptrs, input[(i+f)%n] + f*subpacket_size, subpacket_size);
	}
	free(data_plus_coding_ptrs);
	free(erased);
	return(1);
}

int repair_encode_SRC(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info)
{
	int f = info->req.f;
	int n = info->req.n;
	int subpacket_size = input_size/(f+1);		

	// repair encoding is a simple copy operation
	int distance = (from_device_ID-to_device_ID+n)%n;
	int shift = 0;
	if(distance<=f){
		shift = subpacket_size*(f+1-distance);
		memcpy(output,input+subpacket_size*(distance),shift);
	}
	distance = (to_device_ID-from_device_ID+n)%n;
	if(distance<=f)
		memcpy(output+shift,input,subpacket_size*(f+1-distance));	
	
	return(1);
}

int repair_decode_SRC(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info)
{
	int i,j, counter;
	int f = info->req.f;
	int n = info->req.n;
	int d = info->req.d;
	int subpacket_size = input_size/(f+1);
	char **data_ptrs = malloc(sizeof(char*)*f);	
	int distance,shift;
	
	if(data_ptrs==NULL){
		printf("Can not allocate memory.\n");
		return(-1);
	}
	for(i=0,counter=0;i<n&&helpers[i]>=0;i++){
		if((helpers[i]-to_device_ID+n)%n>f&&(to_device_ID-helpers[i]+n)%n>f)
			continue;
		counter++;
	}
	if(counter!=d){
		printf("Insufficient number of helpers.\n");
		return(-1);
	}
	memset(output,0,output_size);
	for(j=0;j<f+1;j++){			
		for(i=0,counter=0;i<d;i++){
			distance = (helpers[i]-to_device_ID+n)%n;
			shift = 0;
			if(distance<=f)
				shift = subpacket_size*(f+1-distance);			
			if(j+distance<=f){				
				data_ptrs[counter] = input[i]+subpacket_size*j;
				counter++;
			}
			distance = (to_device_ID-helpers[i]+n)%n;
			if(j-distance>=0){
				data_ptrs[counter] = input[i]+subpacket_size*(j-distance)+shift;
				counter++;
			}
			
		}		
		jerasure_do_parity(f, data_ptrs, (output+j*subpacket_size), subpacket_size);
	}
	
	free(data_ptrs);
	
	return(1);
}
