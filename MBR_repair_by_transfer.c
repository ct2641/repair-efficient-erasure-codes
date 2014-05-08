/*
# RegeneratingCodes/MBR_Repair_by_transfer.c

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
#include "jerasure_add.h"
#include "reed_sol.h"

#include "regenerating_codes.h"
#include "cauchy.h"


/*
The code is based on the paper "Distributed Storage Codes With Repair-by-Transfer and Nonachievability of Interior Points..."
by Shah, Rashmi, Kumar and Ramchandran, Trans. IT, Mar. 2012.

The coding has two steps: the first step is a conventional Reed-Solomon code, and the second step involves
a careful symbol placement arrangement. The placement is not unique, and here we use the pattern
disk indices:  1    2    3    4  ....  n
symbols:     | c1   c1   c2   c3       c(n-1)|
             | c2   cn   cn   c(n+1) ....
             | c3   c(n+1)...
             | c4   ...

             ...
     
             | c(n-1) ....

In other words, the matrix is determined by the upper trianglar matrix excluding the diagonal elements,
the the first row is tranposed and placed in the first column, the second row is place into the second column...
This pattern has the advantage of being easily mapped back and forth between symbols.

Moreover, to simplify disk I/O, we group the symobls with the same subsript together in a packet-group: each packet is of size
packetsize/(k(n-1)-k(k-1)/2)*(n-1), and each subpacket is of size=packetsize/(k(n-1)-k(k-1)/2). 

We first initialize these subpackets, call the conventional RS encoder routine, and then replicate within the packet group. 

Restriction: this construction works only for d=n-1. The alphabet size 2^w >= n*(n-1)/2 
*/


int encode_MBR_repair_by_transfer(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info)
{
	int i,j,counter;	
	int n = info->req.n;
	int subpacket_size = input_size/info->req.inner_k;
	char **data_plus_coding_ptrs = malloc(sizeof(char*)*info->req.inner_n);        
  
	// rearrange the memory pointers in preparation for encoding
	if(data_plus_coding_ptrs==NULL){
		printf("Out of memory.\n");
		return(-1);
	}	
	for(counter=0,i=0; i<n;i++){ // i-th device, or i-th column	
		for(j=i;j<n-1;j++){ // j-th row, or j-th subpacket on the i-th device	 		
			data_plus_coding_ptrs[counter] = output[i]+j*subpacket_size;	
			counter++;				
		}
	}	
        
	// now copy the data content into the output buffer
	for(counter=0;counter<info->req.inner_k;counter++)
		memcpy(data_plus_coding_ptrs[counter],input+subpacket_size*counter,subpacket_size);	
        
	// call jerasure routine for encoding;
	jerasure_schedule_encode(info->req.inner_k, info->req.inner_n-info->req.inner_k, info->req.w, 
				info->schedule, data_plus_coding_ptrs, 
				(data_plus_coding_ptrs+info->req.inner_k), 
				subpacket_size, ALIGNMENT);
	
	// now replicate data using the symbol placement pattern specified above	
	for(counter=0,j=0; j<n-1;j++){ // j-th subpacket, or j-th row
		for(i=j+1;i<n;i++){ // i-th column, or i-th device				
			memcpy(output[i]+j*subpacket_size,data_plus_coding_ptrs[counter],subpacket_size);	
			counter++;				
		}
	}
	// clean up
	free(data_plus_coding_ptrs);       
	return(1);
}

int decode_MBR_repair_by_transfer(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info)
{
	int i, j, counter = 0, num_erasures=0;
	int n = info->req.n;
	int subpacket_size = input_size/(n-1);	
	int* erased;	
	int* pseudo_erasures=calloc(info->req.inner_n,sizeof(int));
	if(pseudo_erasures==NULL){		
		printf("Out of memory.\n");
		return(-1);
	} 
       
	// allocate memory for data arrangement
	char** data_plus_coding_ptrs = malloc(sizeof(void*)*info->req.inner_n);
	if(data_plus_coding_ptrs==NULL){
		printf("Out of memory.\n");
		return(-1);
	}
	// convert erasures into erased format for easy processing 
	erased=jerasure_erasures_to_erased(n, n-info->req.k, erasures);
	if(erased==NULL){
		printf("Too many erasures, can not recover.\n");
		goto complete;
	}
	// assign buffer pointers in preparation for decoding	
	for(counter=0,i=0;i<n;i++){  //i-th device	 		
		if(erased[i]==1){
			for (j=i;j<n-1;j++){	// j-th row			
				data_plus_coding_ptrs[counter] = input[i]+j*subpacket_size;
				if(erased[j+1]==1) // yes, this subpacket is indeed lost
				{
					pseudo_erasures[num_erasures] = counter;
					num_erasures++;					
				}
				else
					memcpy(data_plus_coding_ptrs[counter],input[j+1]+i*subpacket_size,subpacket_size);
				counter++;
			}
		}
		else{
			for (j=i;j<n-1;j++){	// j-th row			
				data_plus_coding_ptrs[counter] = input[i]+j*subpacket_size;			
				counter ++;
			}
		}
	}
	pseudo_erasures[num_erasures] = -1;

	// call jerasure routine for decoding
	jerasure_schedule_decode_lazy(info->req.inner_k, info->req.inner_n-info->req.inner_k, info->req.w, 
				info->bitmatrix, pseudo_erasures, 
				data_plus_coding_ptrs, 
				(data_plus_coding_ptrs+info->req.inner_k), 	
			        subpacket_size, ALIGNMENT, 0);

	// replicate data using the symbol placement pattern, which repairs all the lost devices also	
	for(j=0, counter = 0; j<n-1;j++){ // j-th subpacket, or j-th row
		for(i=j+1;i<n;i++){ // i-th column, or i-th device				
			memcpy(input[i]+j*subpacket_size,data_plus_coding_ptrs[counter],subpacket_size);			
			counter++;				
		}
	}	
	for(counter=0;counter<info->req.inner_k;counter++)
		memcpy(output+subpacket_size*counter,data_plus_coding_ptrs[counter],subpacket_size);	
	
	// clean up
	free(erased);
complete:
	free(data_plus_coding_ptrs);
	free(pseudo_erasures);
	
	return(1);
}

int repair_encode_MBR_repair_by_transfer(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info)
{
	int n = info->req.n;
	int subpacket_size = input_size/(n-1);	
	if(subpacket_size!=output_size){
		printf("Incorrect buffer size.\n");
		return(-1);
	}
	// repair encoding is a simple copy operation
	if(from_device_ID<to_device_ID)
		memcpy(output,input+subpacket_size*(to_device_ID-1),subpacket_size);	
	else
		memcpy(output,input+subpacket_size*to_device_ID,subpacket_size);		
	
	return(1);
}

int repair_decode_MBR_repair_by_transfer(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info)
{
	int i,counter;
	int n = info->req.n;
	int subpacket_size = output_size/(n-1);
	int* helpers_inv_ID = malloc(sizeof(int)*n);
	for(i=0,counter=0;i<n-1&&helpers[i]>=0;i++){		
		counter++;
		if(helpers[i]<to_device_ID)
			helpers_inv_ID[helpers[i]] = i;
		else
			helpers_inv_ID[helpers[i]-1] = i;
	}
	if(counter<info->req.d){
		printf("Insufficient number of helpers.\n");
		return(-1);
	}

	// repair decoding is also a simple copy operation following the right order
	for(i=0;i<n-1;i++)        
		memcpy(output+subpacket_size*i,input[helpers_inv_ID[i]],subpacket_size);	
	free(helpers_inv_ID);
	return(1);
}
