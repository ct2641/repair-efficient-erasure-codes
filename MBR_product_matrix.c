/*
# RegeneratingCodes/MBR_product_matrix.c

Copyright (c)2014, AT&T Intellectual Property.  All other rights reserved.

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
The code is based on the paper "Optimal Exact-Regenerating Codes for Distributed Storage at the MSR and MBR Points via a Product-Matrix Construction"
by K. V. Rashmi, Nihar B. Shah, P. Vijay Kumar, Trans. IT, Aug. 2011.

The message matrix is arranged as M =[S  T
                                      T' 0],
where S is a symmetric matrix of size k-by-k, and T is a matrix of size k-by-(d-k). It has a total of k(k+1)/2+k(d-k) info symbols. 
This message matrix is not explicitly formed, but only used through arrangement of pointers. We fill the top-triangle of S and T matrix together row-wise,

The encoding matrix is G=[I 0 
                          A B], 
where I is a k-by-k matrix and [A B] together is a (n-k)-by-d Cauchy matrix.

The coded symbol is then G*M, i.e., an n-by-d matrix, each row corresponds to a device.

We group the subpackets together in order to simplify the disk io where each device has d subpackets.  

Restriction: alphabet size 2^w>=(n-k+d)

*/


int encode_MBR_product_matrix(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info)
{
	int i,j, counter;	
	int n = info->req.n;
	int d = info->req.d;
	int k = info->req.k;
	int subpacket_size = input_size/(k*(k+1)/2+k*(d-k));
	char **data_ptrs = malloc(sizeof(char*)*d*d); // this is the pointer matrix for the message matrix M
	char **coding_ptrs = malloc(sizeof(char*)*(n-k)); // this is the pointers for holding encoded pieces
	if(data_ptrs==NULL||coding_ptrs==NULL){
		printf("Cannot allocate memory!\n");
		return(-1);
	}

	for(counter=0,i=0;i<d;i++){ // this is the d column in the matrix M,
		for(j=0;j<i;j++)
			data_ptrs[j+i*d] = data_ptrs[i+j*d];
		for(j=i;j<d;j++)
			data_ptrs[j+i*d] = input+subpacket_size*(counter+j-i);	
		counter += d-i;		
	}	

	for(i=0;i<k;i++){
		for(j=0;j<n-k;j++)
			coding_ptrs[j] = output[j+k]+subpacket_size*i;
		jerasure_schedule_encode(d, n-k, info->req.w, 
				info->schedule,  data_ptrs+d*i, 
				coding_ptrs, 
				subpacket_size, ALIGNMENT);		
		for(j=0;j<k;j++)
			memcpy(output[j]+subpacket_size*i,data_ptrs[d*i+j],subpacket_size);
	}
	for(i=k;i<d;i++){
		for(j=0;j<n-k;j++)
			coding_ptrs[j] = output[j+k]+subpacket_size*i;		
		jerasure_schedule_encode(k, n-k, info->req.w, 
				info->subschedule_array[0],  data_ptrs+d*i, coding_ptrs, 
				subpacket_size, ALIGNMENT);
		for(j=0;j<k;j++)
			memcpy(output[j]+subpacket_size*i,data_ptrs[d*i+j],subpacket_size);
	}	
	// clean up
	free(data_ptrs);       
	free(coding_ptrs);       
	return(1);
}

int decode_MBR_product_matrix(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info)
{
	int i, j,counter = 0;
	int n = info->req.n;
	int d = info->req.d;
	int k = info->req.k;
	int subpacket_size = input_size/d;
	int num_of_long = subpacket_size/sizeof(long);	
	long *src_pos,*des_pos;
       
	// allocate memory for data arrangement
	char** data_plus_coding_ptrs = malloc(sizeof(char*)*n);
	char** to_be_XORed = malloc(sizeof(char*)*(n-k));
	for(i=0;i<n-k;i++)
		to_be_XORed[i] = malloc(subpacket_size);
	if(data_plus_coding_ptrs==NULL||to_be_XORed==NULL||to_be_XORed[0]==NULL){
		printf("Out of memory.\n");
		return(-1);
	}
	// convert erasures into erased format for easy processing 
	int* erased=jerasure_erasures_to_erased(n, n-info->req.k, erasures);
	if(erased==NULL){
		printf("Too many erasures, can not recover.\n");
		goto complete;
	}
	// first decode the T portion of the matrix M, this also repairs the T portion of the coded info
	for(i=0; i<d-k ; i++){
		for(j=0;j<n;j++)
			data_plus_coding_ptrs[j] = input[j]+(i+k)*subpacket_size;
		jerasure_schedule_decode_lazy(k, n-k, info->req.w, 
				info->subbitmatrix_array[0], erasures, 
				data_plus_coding_ptrs, 
				data_plus_coding_ptrs+k, 	
			        subpacket_size, ALIGNMENT, 0);
	}
	
	// next decode the S portion of the matrix M
	for(i=0;i<k;i++){ // i-th column of S		
		// need to complete an XOR operation B*T' to cancel out the decoded part
		// view the matrix B*T' as encoding, so first set up the T information arrays
		for(j=0;j<d-k;j++)
			data_plus_coding_ptrs[j] = input[i]+subpacket_size*(k+j);
		
		jerasure_schedule_encode(d-k, n-k, info->req.w, 
			info->subschedule_array[1], data_plus_coding_ptrs, to_be_XORed, 
			subpacket_size, ALIGNMENT);
		
		for(j=0;j<n-k;j++){
			des_pos = (long*)(input[k+j]+i*subpacket_size);
			src_pos = (long*)(to_be_XORed[j]);
			for(counter=0;counter<num_of_long;counter++,des_pos++,src_pos++)
				*des_pos ^= *src_pos;							
		}

		for(j=0;j<n;j++)
			data_plus_coding_ptrs[j] = input[j]+i*subpacket_size;

		// then we can decode 
		jerasure_schedule_decode_lazy(k, n-k, info->req.w, 
				info->subbitmatrix_array[0], erasures, 
				data_plus_coding_ptrs, 
				data_plus_coding_ptrs+k, 	
			        subpacket_size, ALIGNMENT, 0);
		//write the XOR part back to the correct coded info
		for(j=0;j<n-k;j++){
			des_pos = (long*)(input[k+j]+i*subpacket_size);
			src_pos = (long*)(to_be_XORed[j]);
			for(counter=0;counter<num_of_long;counter++,des_pos++,src_pos++)
				*des_pos ^= *src_pos;	
		}
	}

	// extract data from M matrix and copy to output buffer
	for(counter=0,i=0;i<k;i++){ 
		memcpy(output+counter*subpacket_size,input[i]+i*subpacket_size,subpacket_size*(d-i));
		counter += d-i;		
	}

	// clean up
complete:
	free(data_plus_coding_ptrs);
	free(erased);
	for(i=0;i<n-k;i++)
		free(to_be_XORed[i]);
	free(to_be_XORed);	
	return(1);
}

int repair_encode_MBR_product_matrix(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info)
{
	int d = info->req.d;
	int k = info->req.k;
	int w = info->req.w;
	int subpacket_size = input_size/d;	
	int i;
	if(subpacket_size!=output_size){
		printf("Incorrect buffer size.\n");
		return(-1);
	}
	char** data_ptrs = malloc(sizeof(char*)*d);
	if(data_ptrs==NULL)
		return(-1);
	// repair here needs an encoding step
	if(to_device_ID<k)//just copy that single position
		memcpy(output,input+subpacket_size*to_device_ID,subpacket_size);
	else{ // otherwise need do real computation, but it can be thought as an encoding step 
		for(i=0;i<d;i++)
			data_ptrs[i] = input+i*subpacket_size;
		jerasure_bitmatrix_encode(d,1,w,info->bitmatrix+d*w*(to_device_ID-k)*w,data_ptrs,&output,subpacket_size,ALIGNMENT);
	}	
	free(data_ptrs);	
	return(1);
}
int repair_decode_MBR_product_matrix(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info)
{

	int i,counter;
	int d = info->req.d;	
	int k = info->req.k;
	int w = info->req.w;
	int subpacket_size = output_size/d;
	char** coding_ptrs = malloc(sizeof(void*)*d);
	int *repair_matrix = malloc(sizeof(int)*d*d);
	int *repair_matrix_inv = malloc(sizeof(int)*d*d);

	memset(repair_matrix,0,sizeof(int)*d*d);
	for(i=0,counter=0;i<d&&helpers[i]>=0;i++){
		counter ++;
		if(helpers[i]<k)
			repair_matrix[i*d+helpers[i]]=1;
		else			
			memcpy(repair_matrix+d*i,info->matrix+(helpers[i]-k)*d,sizeof(int)*d);
	}	
	if(counter<d){
		printf("Insufficient number of helpers.\n");
		return(-1);
	}
	jerasure_invert_matrix(repair_matrix,repair_matrix_inv,d,w);

	for(i=0; i<d ;i++)
		coding_ptrs[i] = output + i*subpacket_size;
	int *coding_bitmatrix = jerasure_matrix_to_bitmatrix(d,d,w,repair_matrix_inv);
	jerasure_bitmatrix_encode(d,d,w,coding_bitmatrix,(char**)input,coding_ptrs,subpacket_size,ALIGNMENT);

	free(coding_ptrs);
	free(repair_matrix);
	free(repair_matrix_inv);
	free(coding_bitmatrix);
	return(1);
}
