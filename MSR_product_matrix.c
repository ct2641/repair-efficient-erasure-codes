/*
# RegeneratingCodes/MBR_ProductMatrix.c

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
#include "include/jerasure.h"
#include "jerasure_add.h"
#include "include/reed_sol.h"
#include <time.h>

#include "regenerating_codes.h"
#include "include/cauchy.h"
/*
The code is based on the paper "A Unified Form of Exact-MSR Codes via Product-Matrix Framework" 
by S.J. Lin and W.H. Chung, in Proc. 2013 IEEE 24th Annual International Symposium on Personal, 
Indoor, and Mobile Radio Communications, which provides a more implementation-friendly version  
of the construction in "Optimal Exact-Regenerating Codes for Distributed Storage at the MSR 
and MBR Points via a Product-Matrix Construction" by K. V. Rashmi, Nihar B. Shah, P. Vijay Kumar, 
Trans. IT, Aug. 2011.

The message matrix is arranged as M =[S1  0
                                      S2  T
                                      T'  Z],
where S1 and S2 are a symmetric matrix of size (k-1)-by-(k-1), and T is a matrix of size (k-1)-by-(d-2k+2), 
Z is a symmetric matrix with zeros in the lower-right (d-2k+1)-by-(d-2k+1) part of the matrix. 
The matrix M has a total of k(d-k+1) info symbols. 

The encoding matrix is G=[\Lambda \Phi, \Phi, \Delta], where \Phi is an n-by-(k-1) matrix, \Lambda is a diagonal
matrix, and \Delta is a n-by-(d-2k+2) matrix.

The coded symbol is then G*M, i.e., an n-by-(d-k+1) matrix, each row corresponds to a device. 
We group the subpackets together in order to simplify the disk io where each device has d subpackets.  

Restriction: alphabet size 2^w>=n.
*/

int decode_MSR_product_matrix_no_output(char **input, size_t input_size, int* erasures, struct coding_info *info);

int encode_MSR_product_matrix(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info)
{
	int i;	
	int n = info->req.n;	
	int k = info->req.k;	
	int *erasures = malloc(sizeof(int)*n);
	for(i=0;i<n-k;i++)	
		erasures[i] = i+k;
	erasures[n-k] = -1;

	for(i=0;i<k;i++)
		memcpy(output[i],input+output_size*i,output_size);		
	if(decode_MSR_product_matrix_no_output(output,output_size,erasures,info)<0){
		free(erasures);	return(-1);
	}
	free(erasures);
	return(1);
}

int decode_MSR_product_matrix_no_output(char **input, size_t input_size, int* erasures, struct coding_info *info)
{
	clock_t clk, tclk;
	int i, j,c1,c2,tdone,inv;
	int n = info->req.n;
	int d = info->req.d;
	int k = info->req.k;
	int w = info->req.w;
	int subpacket_size = input_size/(d-k+1);
	int num_of_long = subpacket_size/sizeof(long);	
	long *src_pos,*des_pos;
	int *vector_A=NULL;
 	char **ptrs;
	int **decode_schedule=NULL;
	int *erased = NULL;	
	int *bitmatrix_temp = malloc(sizeof(int)*(k-1)*(k-1)*w*w*4);
	char *data_transformed = malloc(input_size*k);    // this is the buffer for tranformed data, i.e., matrix M. The output is the systematic part of 
							  // codingmatrix*M, and 
					                  // we will regenerate the erased data from M using the encoding matrix, which will be written to *output.
	int *pseudo_erasures = malloc(sizeof(int)*(n*2));  // in order to recompute M, we view it as a systematic MDS code, 
							  // sometimes (n+k,k), sometimes (n+d-k+1,d-k+1), thus we over allocate
	char **M_ptrs = malloc(sizeof(void*)*d*(d-k+1));  // this is the pointer matrix to elements in M.		
	char **data_ptrs = malloc(sizeof(void*)*n);
	char **coding_ptrs = malloc(sizeof(void*)*n);
	int* remaining = malloc(sizeof(int)*k);           // not erased devices
	char *buffer1 = malloc(subpacket_size*k*(k-1)*2);  // need subpacketsize*k>=(max(4,k-1)+k-1)*sizeof(int), thus put factor of 2 to guarentee it
	int *buffer1_int = (int*)buffer1;                  // alternative pointer for buffer1
	char *buffer2 = malloc(subpacket_size*k*(k-1));

	if(data_transformed==NULL||pseudo_erasures==NULL||data_ptrs==NULL
		||coding_ptrs==NULL||buffer1==NULL||buffer2==NULL||M_ptrs==NULL||remaining==NULL){
		printf("Can not allocate memory\n");
		jerasure_free_schedule(decode_schedule);
		if(data_ptrs!=NULL)free(data_ptrs);
		if(coding_ptrs!=NULL)free(coding_ptrs);
		if(pseudo_erasures!=NULL)free(pseudo_erasures);
		if(data_transformed!=NULL)free(data_transformed);
		if(M_ptrs!=NULL)free(M_ptrs);	
		if(buffer1!=NULL)free(buffer1);
		if(buffer2!=NULL)free(buffer2);
		if(remaining!=NULL)free(remaining);
		return(-1);
	}
	//set up pointers for matrix M
	// first k-1 rows, only have S1
	for(i=0,c2=0;i<k-1;i++){
		c1 = i*(d-k+1);
		for(j=i;j<k-1;j++,c2++)
			M_ptrs[c1+j] = data_transformed + subpacket_size*c2;
		for(j=k-1;j<d-k+1;j++)
			M_ptrs[c1+j] = NULL;			
		for(j=0;j<i;j++)
			M_ptrs[c1+j] = M_ptrs[j*(d-k+1)+i]; // symmetric matrix, thus there is (i,j) and (j,i) point to the same pointer.
	}
	// next k-1 rows
	for(i=0;i<k-1;i++){
		c1 = (k-1+i)*(d-k+1);
		for(j=i;j<d-k+1;j++,c2++)
			M_ptrs[c1+j] = data_transformed + subpacket_size*c2;
		for(j=0;j<i;j++)
			M_ptrs[c1+j] = M_ptrs[(j+k-1)*(d-k+1)+i];		
	}
	// next 1 and (d-2k+1) rows, may not exist
	if(d>2*k-2)
	{
		c1 = (2*k-2)*(d-k+1);
		for(j=k-1;j<d-k+1;j++,c2++)
			M_ptrs[c1+j] = data_transformed + subpacket_size*c2;
		for(j=0;j<i;j++)
			M_ptrs[c1+j] = M_ptrs[(j+k-1)*(d-k+1)+k-1];			
		for(i=2*k-1;i<d;i++){
			c1 = i*(d-k+1);
			for(j=0;j<k;j++)
				M_ptrs[c1+j] = M_ptrs[(j+k-1)*(d-k+1)+i-k+1];
			for(j=k;j<d-k+1;j++)
				M_ptrs[c1+j] = NULL;
		}
	}
	
        // first decode the last d-2k+1 columns of T and Z: view it as an (n+k,k) MDS code. Note strictly speaking this might 
	// not be a real (n+k,k) MDS code, but it hardly matters.
	// before decoding operation, prepare for the pseudoerasure location array
	if(d>2*k-2){ // only when d>2k-2, T and Z exist
		for(i=0;i<k;i++)
			pseudo_erasures[i] = i;
		for(i=0;i<n;i++){
			if(erasures[i]==-1){
				pseudo_erasures[i+k] = -1;
				break;
			}
			pseudo_erasures[i+k] = erasures[i] + k;		
		}

		// make the decoding schedule: we can save the trouble of manually generate this schedule, at the expense of
		// more computation
		decode_schedule = jerasure_generate_decoding_schedule(k, n, w, info->subbitmatrix_array[0], pseudo_erasures, 1);

		for(i=0;i<d-2*k+1;i++){
			for(j=0;j<k;j++)
				data_ptrs[j] = M_ptrs[i+k+(d-k+1)*(k-1+j)];		
			for(j=0;j<n;j++)
				coding_ptrs[j] = input[j]+subpacket_size*(k+i);
			ptrs = set_up_ptrs_for_scheduled_decoding(k, n, pseudo_erasures, data_ptrs,coding_ptrs);
		  	if (ptrs == NULL){
				printf("Can not allocate memory\n");
				goto complete;
			}
			// assume packetsize = ALIGNMENT
			for (tdone = 0; tdone < subpacket_size; tdone += ALIGNMENT*w) {
				jerasure_do_scheduled_operations(ptrs, decode_schedule, ALIGNMENT);
				for (c1 = 0; c1 < k+n; c1++) ptrs[c1] += (ALIGNMENT*w);
			}
			free(ptrs);
		} 
		//next decode the first column of T and Z: we view this as an (n+d-k+1,d-k+1) erasure codes
		// first setup the pseudo erasure location vector	
		for(i=0;i<k;i++)
			pseudo_erasures[i] = i; // in the first columne of the Z matrix, the last d-2k+1 elements are known, but the first k elements are not
		for(i=0;i<n;i++)
		{
			if(erasures[i]==-1){
				pseudo_erasures[i+k] = -1;
				break;
			}
			pseudo_erasures[i+k] = erasures[i]+d-k+1;		
		}
		for(j=0;j<d-k+1;j++)
			data_ptrs[j] = M_ptrs[(d-k+1)*(k-1+j)+k-1];
		for(j=0;j<n;j++)
			coding_ptrs[j] = input[j]+subpacket_size*(k-1);
		jerasure_schedule_decode_lazy(d-k+1,n,w,
					info->subbitmatrix_array[1],pseudo_erasures,data_ptrs,
					coding_ptrs,subpacket_size,ALIGNMENT,0);
	}
	//clk = clock();

	// now this is the hard part: to decode S1 and S2. The algorithm used here is slightly different from that in the paper:
	// instead of right multiply \Phi_{DC}', we only right multiply the sub-matrix of \Phi_{DC}' without of the last row

	// setting up remaining device array to facilitate decoding
	erased = jerasure_erasures_to_erased(k, n-k, erasures);
	for(i=0,c1=0;i<n&&c1<k;i++){
		if(erased[i]==0){
			remaining[c1] = i;
			c1++;
		}	
	}	
	//compute C_{DC}-\Delta_{DC}*T'
	for(i=0;i<k-1;i++){ // has k-1 columns
		for(j=0;j<d-2*k+2;j++)
			data_ptrs[j] = M_ptrs[(2*k-2+j)*(d-k+1)+i];
		for(j=0;j<k;j++)
			coding_ptrs[j] = buffer1+(j*(k-1)+i)*subpacket_size;
		for(j=0;j<k;j++){
			jerasure_bitmatrix_dotprod(d-2*k+2, w, info->subbitmatrix_array[2]+remaining[j]*(d-2*k+2)*w*w, NULL, j+d-2*k+2,
        	                data_ptrs, coding_ptrs, subpacket_size, ALIGNMENT);
			src_pos = (long*)(input[remaining[j]]+i*subpacket_size);
			des_pos	= (long*)(coding_ptrs[j]);
			for(c2=0;c2<num_of_long;c2++)	
				des_pos[c2] ^= src_pos[c2];
		}	
	}
	
	// right multiply \Phi_{DC}': this will be P. 
	// result is in buffer2
	for(j=0;j<k;j++){ //j-th row
		for(i=0;i<k-1;i++)
			data_ptrs[i] = buffer1+(j*(k-1)+i)*subpacket_size;
		for(i=0;i<k-1;i++)
			coding_ptrs[i] = buffer2 +(j*(k-1)+i)*subpacket_size;
		for(i=0;i<k-1;i++){
			jerasure_bitmatrix_dotprod(k-1, w, info->subbitmatrix_array[3]+remaining[i]*(k-1)*w*w, NULL, i+k-1,
        	                data_ptrs, coding_ptrs, subpacket_size, ALIGNMENT);
		}
	}
	// now solve for the off-diagonal terms
	for(i=0;i<k-1;i++){
		for(j=i+1;j<k-1;j++){
			// solve for S1 tilde off-diagonal 
			// here we directly use the fact that Lambda = [0 1 2 3 ....];	
			int** temp_schedule;			
			inv = galois_single_divide(1,remaining[i]^remaining[j],w);					
			buffer1_int[0] = buffer1_int[1] = inv;			
			data_ptrs[0] = buffer2+(i*(k-1)+j)*subpacket_size;
			data_ptrs[1] = buffer2+(j*(k-1)+i)*subpacket_size;
			coding_ptrs[0] = M_ptrs[i*(d-k+1)+j];
			jerasure_matrix_to_bitmatrix_noallocate(2,1,w,buffer1_int,bitmatrix_temp);
			temp_schedule = jerasure_smart_bitmatrix_to_schedule(2, 1, w, bitmatrix_temp);
			jerasure_schedule_encode(2, 1, w, temp_schedule, data_ptrs, coding_ptrs,subpacket_size, ALIGNMENT);	
			if(temp_schedule!=NULL){
				jerasure_free_schedule(temp_schedule);
				temp_schedule = NULL;
			}			
			// solve for S2 tilde off-diagonal
			buffer1_int[0] = galois_single_multiply(remaining[j],inv,w);
			buffer1_int[1] = galois_single_multiply(remaining[i],inv,w);
			coding_ptrs[0] = M_ptrs[(i+k-1)*(d-k+1)+j];
			jerasure_matrix_to_bitmatrix_noallocate(2,1,w,buffer1_int,bitmatrix_temp);
			temp_schedule = jerasure_smart_bitmatrix_to_schedule(2, 1, w, bitmatrix_temp);
			jerasure_schedule_encode(2, 1, w, temp_schedule, data_ptrs, coding_ptrs,subpacket_size, ALIGNMENT);	
			if(temp_schedule!=NULL){
				jerasure_free_schedule(temp_schedule);
				temp_schedule = NULL;
			}			
		}	
	}
	//tclk = clock()-clk;
	//printf("~S1 and ~S2 off-diagonal decoded %.3e clocks \n", (double)tclk);	
	//clk = clock();
	// compute the A vector: A*\Phi_{DC1} = \Phi_{DC2}, note this is always possible because \Phi_{DC1} is alway full rank by construction
	// we first reuse buffer1 here to form \Phi_{DC1} matrix and then the compute its inverse
	for(i=0;i<k-1;i++){
		memcpy(buffer1_int+(k-1)*(k-1+i),(void*)(info->matrix+remaining[i]*d+k-1),(k-1)*sizeof(int));
	}

	jerasure_invert_matrix(buffer1_int+(k-1)*(k-1),buffer1_int,k-1,w);	
	vector_A = jerasure_matrix_multiply(info->matrix+remaining[k-1]*d+k-1,buffer1_int,1,k-1,k-1,k-1,w);	
	if(vector_A==NULL)
		goto complete;	
	
	for(i=0;i<k-1;i++){
		buffer1_int[i+(k-1)*(k-1)] = galois_single_multiply(vector_A[i],remaining[k-1],w);		
		buffer1_int[i+k*(k-1)] = vector_A[i];
	}		

	for(i=0;i<k-1;i++){
		memset(buffer1_int+(k-1)*(k+1),0,sizeof(int)*2*(k-1));
		buffer1_int[(k-1)*(k+1)+i] = remaining[i];
		buffer1_int[(k-1)*(k+2)+i] = 1;

		pseudo_erasures[0] = i;
		pseudo_erasures[1] = k-1+i;
		pseudo_erasures[2] = -1;
		for(j=0;j<2*k-2;j++)
			data_ptrs[j] = M_ptrs[i+j*(d-k+1)];
		coding_ptrs[0] = buffer2+((k-1)*(k-1)+i)*subpacket_size;
		coding_ptrs[1] = buffer2+(i*(k-1)+i)*subpacket_size;

		jerasure_matrix_to_bitmatrix_noallocate(2*k-2,2,w,buffer1_int+(k-1)*(k-1),bitmatrix_temp);				
		jerasure_schedule_decode_lazy(2*k-2,2,w,bitmatrix_temp,pseudo_erasures,data_ptrs,coding_ptrs,subpacket_size,ALIGNMENT,0);
      	}
	//tclk = clock()-clk;
	//printf("~S1 and ~S2 decoded %.3e clocks \n", (double)tclk);	
	// now we have both \tilde{S_1} and \tilde{S_2} in M_ptrs, need to recover S1 and S2 from them
	// this is done by multiply \tilde{S_1} left and right by inv(\Phi_{DC1}).
	// right-multiply for S1 
	int* bitmatrix_inv = jerasure_matrix_to_bitmatrix(k-1,k-1,w,buffer1_int);
	int** inv_schedule = jerasure_smart_bitmatrix_to_schedule(k-1, k-1, w, bitmatrix_inv);

	// right-multiply for S1
	for(i=0;i<k-1;i++){
		for(j=0;j<k-1;j++)
			data_ptrs[j] = M_ptrs[i*(d-k+1)+j];		
		for(j=0;j<k-1;j++)
			coding_ptrs[j] = buffer2+(i*(k-1)+j)*subpacket_size;	
		jerasure_schedule_encode(k-1, k-1, w, inv_schedule, data_ptrs, coding_ptrs, subpacket_size, ALIGNMENT);	
	}
	// left-multiply for S1 
	for(j=0;j<k-1;j++){
		for(i=0;i<k-1;i++)
			data_ptrs[i] = buffer2+(i*(k-1)+j)*subpacket_size;
		for(i=0;i<k-1;i++)
			coding_ptrs[i] = M_ptrs[i*(d-k+1)+j];

		jerasure_schedule_encode(k-1, k-1, w, inv_schedule, data_ptrs, coding_ptrs, subpacket_size, ALIGNMENT);

	}
	// right-multiply for S2
	for(i=0;i<k-1;i++){
		for(j=0;j<k-1;j++)
			data_ptrs[j] = M_ptrs[(i+k-1)*(d-k+1)+j];
		for(j=0;j<k-1;j++)
			coding_ptrs[j] = buffer2+(i*(k-1)+j)*subpacket_size;

		jerasure_schedule_encode(k-1, k-1, w, inv_schedule, data_ptrs, coding_ptrs, subpacket_size, ALIGNMENT);
	}
	// left-multiply for S2 
	for(j=0;j<k-1;j++){
		for(i=0;i<k-1;i++)
			data_ptrs[i] = buffer2+(i*(k-1)+j)*subpacket_size;
		//for(i=j;i<k-1;i++)
		for(i=0;i<k-1;i++)
			coding_ptrs[i] = M_ptrs[(i+k-1)*(d-k+1)+j];

		jerasure_schedule_encode(k-1, k-1, w, inv_schedule, data_ptrs, coding_ptrs, subpacket_size, ALIGNMENT);
	}
	// having S1,S2,T, now can also fill the first k-1 column of the output		
	for(i=0;i<k-1;i++){
		for(j=0;j<d;j++)
			data_ptrs[j] = M_ptrs[j*(d-k+1)+i];
		for(j=0;j<n;j++)
			coding_ptrs[j] = input[j]+i*subpacket_size;
		for(j=0;j<n;j++){
			if(erased[j]==1)
				jerasure_bitmatrix_encode(d,1,w,info->bitmatrix+(j*d*w*w),data_ptrs,coding_ptrs+j,subpacket_size,ALIGNMENT);
		}
	}

	// clean up
complete:
	if(decode_schedule)
		jerasure_free_schedule(decode_schedule);
	if(inv_schedule)
		jerasure_free_schedule(inv_schedule);
	free(data_ptrs);
	free(coding_ptrs);
	if(erased)free(erased);
	free(pseudo_erasures);
	free(data_transformed);
	free(M_ptrs);	
	free(buffer1);
	free(buffer2);
	free(remaining);
	free(bitmatrix_temp);
	if(vector_A!=NULL)free(vector_A);
	return(1);
}

int decode_MSR_product_matrix(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info)
{
	int i;
	int k = info->req.k;
	if(decode_MSR_product_matrix_no_output(input, input_size, erasures, info)<0)
		return(-1);
	for(i=0;i<k;i++,output+=input_size)
		memcpy(output,input[i],input_size);
	return(1);
}

int repair_encode_MSR_product_matrix(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info)
{
	int d = info->req.d;
	int k = info->req.k;
	int w = info->req.w;
	int subpacket_size = input_size/(d-k+1);	
	int i;
	if(subpacket_size!=output_size){
		printf("Incorrect buffer size.\n");
		return(-1);
	}
	char** data_ptrs = malloc(sizeof(void*)*d);
	if(data_ptrs==NULL)
		return(-1);
	
	for(i=0;i<d-k+1;i++)
		data_ptrs[i] = input+i*subpacket_size;
	jerasure_bitmatrix_encode(d-k+1,1,w,info->subbitmatrix_array[1]+(d-k+1)*w*to_device_ID*w,data_ptrs,(char**)(&output),subpacket_size,ALIGNMENT);
	
	free(data_ptrs);	
	return(1);
}
int repair_decode_MSR_product_matrix(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info)
{
	int i,counter;
	int d = info->req.d;	
	int k = info->req.k;
	int w = info->req.w;
	int subpacket_size = input_size;
	char** coding_ptrs = malloc(sizeof(void*)*(d-k+1));
	int *repair_matrix = malloc(sizeof(int)*d*d);
	int *repair_matrix_inv = malloc(sizeof(int)*d*d);
	int *combination_matrix = malloc(sizeof(int)*(d-k+1)*d);

	for(i=0,counter=0;i<d&&helpers[i]>=0;i++,counter++)
		memcpy(repair_matrix+d*i,info->matrix+d*helpers[i],sizeof(int)*d);
	if(counter<d){
		printf("Insufficient number of helpers.\n");
		return(-1);
	}
	jerasure_invert_matrix(repair_matrix,repair_matrix_inv,d,w);
	memset(combination_matrix,0,sizeof(int)*(d-k+1)*d);
	for(i=0;i<k-1;i++){
		*(combination_matrix+(d*i)+i) = to_device_ID;
		*(combination_matrix+(d*i)+i+k-1) = 1;
	}
	for(i=k-1;i<d-k+1;i++)
		*(combination_matrix+(d*i)+i+k-1) = 1;

	int *coding_matrix = jerasure_matrix_multiply(combination_matrix,repair_matrix_inv,d-k+1,d,d,d,w);
		
	for(i=0; i<d-k+1 ;i++)
		coding_ptrs[i] = output + i*subpacket_size;
	int *coding_bitmatrix = jerasure_matrix_to_bitmatrix(d,d-k+1,w,coding_matrix);
	jerasure_bitmatrix_encode(d,d-k+1,w,coding_bitmatrix,(char**)input,coding_ptrs,subpacket_size,ALIGNMENT);

	free(coding_ptrs);
	free(coding_matrix);
	free(repair_matrix);
	free(repair_matrix_inv);
	free(combination_matrix);
	free(coding_bitmatrix);
	return(1);
}
