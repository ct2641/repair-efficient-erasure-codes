/*
# RegeneratingCodes/RegeneratingCodes.c

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
#include <math.h>
#include "jerasure.h"
#include "reed_sol.h"
#include "cauchy.h"

#include "regenerating_codes.h"


int make_coding_matrics(struct coding_info *info)
{
	int n,k,d,w;
	int i,j;
	int* pointer;
	int beta;

	switch (info->req.type)
	{
		case MBR_PRODUCTMATRIX:
			n = info->req.n;
			d = info->req.d;
			k = info->req.k;	
			w = info->req.w;	
			info->matrix = cauchy_good_general_coding_matrix(d, n-k, w);
  			if (info->matrix == NULL) {
				printf("couldn't make coding matrix.\n");
				return(-1);
			}
			info->bitmatrix = jerasure_matrix_to_bitmatrix(d, n-k, w, info->matrix);
			info->schedule = jerasure_smart_bitmatrix_to_schedule(d, n-k, w,  info->bitmatrix);

			// special for the product matrix based scheme, we need to generate two sub-coding matrices
			info->num_of_submatrices = 2; 
			info->submatrix_array = malloc(sizeof(int*)*info->num_of_submatrices);
			info->subbitmatrix_array = malloc(sizeof(int*)*info->num_of_submatrices);
			info->subschedule_array = malloc(sizeof(int**)*info->num_of_submatrices);
			info->submatrix_array[0] = malloc(sizeof(int)*(n-k)*k);
			info->submatrix_array[1] = malloc(sizeof(int)*(n-k)*(d-k));			
			// now fill the these new matrices			
			for(i=0;i<n-k;i++){
				memcpy(info->submatrix_array[0]+i*k,info->matrix+d*i,sizeof(int)*k);
				memcpy(info->submatrix_array[1]+i*(d-k),info->matrix+d*i+k,sizeof(int)*(d-k));
			}
			// generate bitmatrices and schedules
			info->subbitmatrix_array[0] = jerasure_matrix_to_bitmatrix(k, n-k, w, info->submatrix_array[0]); 
			info->subschedule_array[0] = jerasure_smart_bitmatrix_to_schedule(k, n-k, w,  info->subbitmatrix_array[0]);
			info->subbitmatrix_array[1] = jerasure_matrix_to_bitmatrix(d-k, n-k, w, info->submatrix_array[1]);
			info->subschedule_array[1] = jerasure_smart_bitmatrix_to_schedule(d-k, n-k, w,  info->subbitmatrix_array[1]);			
			break;
		case MBR_REPAIRBYTRANSFER:
			k = info->req.inner_k;
			n = info->req.inner_n;
			w = info->req.w;
			info->matrix = cauchy_good_general_coding_matrix(k, n-k, w);
			if (info->matrix == NULL) {
				printf("couldn't make coding matrix.\n");
			}
			info->bitmatrix = jerasure_matrix_to_bitmatrix(k, n-k, w, info->matrix);
			info->schedule = jerasure_smart_bitmatrix_to_schedule(k, n-k, w, info->bitmatrix);
			info->num_of_submatrices = 0;
			break;
		case MSR_PRODUCTMATRIX:
			n = info->req.n;
			d = info->req.d;
			k = info->req.k;
			w = info->req.w;
			// this is the whole encoding matrix
			pointer = malloc(sizeof(int)*n*d); 
			info->matrix = pointer;
			if(info->matrix==NULL)
				return (-1);
			// the first row is [ 0 0 ... 0 1 0 0 ..], i.e., the row of an extended Vandermonde matrix
			pointer[k-1] = 1;
			for(j=0 ; j<k-2; j++)
				pointer[j] = 0;
			for(j=k ; j<d ; j++)
				pointer[j] = 0;
			// generate the other rows of the encoding matrix 
			for(i=1 ; i < n; i++){
				pointer += d;
				beta = galois_single_multiply(i,i,w);
				// THE \Delta\Phi part
				pointer[0] = i;
				for(j=1;j<k-1;j++)
					pointer[j] = galois_single_multiply(pointer[j-1],beta,w);								
				// the \Phi part with the first column of \Delta
				pointer[k-1] = 1;
				for(j=k;j<2*k-1;j++)
					pointer[j] = galois_single_multiply(pointer[j-1],beta,w);
				// the \Delta part without the first column
				for(j=2*k-1;j<d;j++)
					pointer[j] = galois_single_multiply(pointer[j-1],i,w);
			}			
			info->bitmatrix = jerasure_matrix_to_bitmatrix(d, n, w, info->matrix);
			info->schedule = jerasure_smart_bitmatrix_to_schedule(d, n, w, info->bitmatrix);

			info->num_of_submatrices = 4;
			info->submatrix_array = malloc(sizeof(int*)*info->num_of_submatrices);
			info->subbitmatrix_array = malloc(sizeof(int*)*info->num_of_submatrices);
			info->subschedule_array = malloc(sizeof(int**)*info->num_of_submatrices);
			// submatrix0 is the submatrix with the left k columns of the matrix [\Phi \Delta].
			info->submatrix_array[0] = malloc(sizeof(int)*k*n);
			for(i=0; i<n;i++)
				memcpy(info->submatrix_array[0]+i*k,info->matrix+i*d+k-1,sizeof(int)*k);			
			info->subbitmatrix_array[0] = jerasure_matrix_to_bitmatrix(k, n, w, info->submatrix_array[0]);
			info->subschedule_array[0] = jerasure_smart_bitmatrix_to_schedule(k, n, w, info->subbitmatrix_array[0]);			
			// submatrix1 is the submatrix with the columns of the matrix [\Phi \Delta].
			info->submatrix_array[1] = malloc(sizeof(int)*n*(d-k+1));
			for(i=0; i<n;i++)
				memcpy(info->submatrix_array[1]+i*(d-k+1),info->matrix+i*d+k-1,sizeof(int)*(d-k+1));			
			info->subbitmatrix_array[1] = jerasure_matrix_to_bitmatrix(d-k+1, n, w, info->submatrix_array[1]);
			info->subschedule_array[1] = jerasure_smart_bitmatrix_to_schedule(d-k+1, n, w, info->subbitmatrix_array[1]);
			// submatrix2 is the submatrix with the columns of the matrix [\Delta].
			info->submatrix_array[2] = malloc(sizeof(int)*n*(d-2*k+2));
			for(i=0; i<n;i++)
				memcpy(info->submatrix_array[2]+i*(d-2*k+2),info->matrix+i*d+2*k-2,sizeof(int)*(d-2*k+2));			
			info->subbitmatrix_array[2] = jerasure_matrix_to_bitmatrix(d-2*k+2, n, w, info->submatrix_array[2]);
			info->subschedule_array[2] = jerasure_smart_bitmatrix_to_schedule(d-2*k+2, n, w, info->subbitmatrix_array[2]);
			// submatrix3 is the submatrix with the columns of the matrix [\Phi].
			info->submatrix_array[3] = malloc(sizeof(int)*n*(k-1));
			for(i=0; i<n;i++)
				memcpy(info->submatrix_array[3]+i*(k-1),info->matrix+i*d+k-1,sizeof(int)*(k-1));			
			info->subbitmatrix_array[3] = jerasure_matrix_to_bitmatrix(k-1, n, w, info->submatrix_array[3]);
			info->subschedule_array[3] = jerasure_smart_bitmatrix_to_schedule(k-1, n, w, info->subbitmatrix_array[3]);
			break;
		case SRC:
			n = info->req.n;
			k = info->req.k;	
			w = info->req.w;	
			info->matrix = cauchy_good_general_coding_matrix(k, n-k, w);
			if (info->matrix == NULL) {
				printf("couldn't make coding matrix.\n");
			}
  			info->bitmatrix = jerasure_matrix_to_bitmatrix(k, n-k, w, info->matrix);
			info->schedule = jerasure_smart_bitmatrix_to_schedule(k, n-k, w, info->bitmatrix);
			info->num_of_submatrices = 0;
			break;
		case LRC:
			n = info->req.n;
			k = info->req.k;	
			w = info->req.w;	
			info->matrix = cauchy_good_general_coding_matrix(k, n-k, w);
			if (info->matrix == NULL) {
				printf("couldn't make coding matrix.\n");
			}
			info->bitmatrix = jerasure_matrix_to_bitmatrix(k, n-k, w, info->matrix);
			info->schedule = jerasure_smart_bitmatrix_to_schedule(k, n-k, w, info->bitmatrix);
			info->num_of_submatrices = 0;
			break;
		case STEINERCODE:
			break;
		default: 
			printf("This type of regenerating code is not supported. \n");
			return(-1);	
	}
	return(1);
}

void cleanup_matrics(struct coding_info *info)
{
	int i;

	if(info->matrix!=NULL)
		free(info->matrix);
	if(info->bitmatrix!=NULL)
		free(info->bitmatrix);
	if(info->schedule!=NULL)
		jerasure_free_schedule(info->schedule);

	if(info->num_of_submatrices>0){
		for(i=0;i<info->num_of_submatrices;i++){
			free(info->submatrix_array[i]);
			free(info->subbitmatrix_array[i]);
			jerasure_free_schedule(info->subschedule_array[i]);
		}
		free(info->submatrix_array);
		free(info->subbitmatrix_array);
		free(info->subschedule_array);	
	}

}

int get_requirement(enum codetype type, struct requirement *req, int n, int k, int d, int w)
{
	if(n<3||k>=n||k<3)
	{
		printf("This type of regenerating code is not supported. \n");
		return(-1);
	}
	switch (type)
	{
		case MBR_REPAIRBYTRANSFER: 
			if(d>=n||k>=d){
				printf("This type of regenerating code is not supported. \n");
				return(-1);
			}
			req->inner_n = n*(n-1)/2;
			if((1<<w)<req->inner_n)
			{
				printf("invalid w values.\n");
				return(-1);
			}
			req->type = type;
			req->inner_k = k*n-k*(k+1)/2;
			req->multiple_of = req->inner_k*ALIGNMENT*w;	
			req->min_size = req->multiple_of;
			req->max_size = MIN(req->multiple_of*1024*1024*8,MAXPACKETSIZE);		
			req->n = n;
			req->k = k;
			req->d = n-1;
			req->w = w;
			break;						
		case MSR_PRODUCTMATRIX:
			if(n<=d||d<2*k-2||d<k||(1<<w)<n||k<=2)
			{
				printf("invalid n=%d,k=%d,d=%d,w=%d values.\n",n,k,d,w);	
				return(-1);
			}
			req->n = n;
			req->k = k;
			req->d = d;				
			req->type = type;			
			req->multiple_of = k*(d-k+1)*ALIGNMENT*w;	
			req->min_size = req->multiple_of;
			req->max_size = MIN(req->multiple_of*1024*1024*8,MAXPACKETSIZE);		
			req->w = w;
			break;
		case MBR_PRODUCTMATRIX:
			req->n = n;
			req->k = k;
			req->d = d;	
			if((1<<w)<req->n)
			{
				printf("invalid w values.\n");
				return(-1);
			}
			req->type = type;			
			req->multiple_of = ((k+1)*k/2+k*(d-k))*ALIGNMENT*w;	
			req->min_size = req->multiple_of;
			req->max_size = MIN(req->multiple_of*1024*1024*8,MAXPACKETSIZE);		
			req->w = w;
			break;
		case SRC:
			req->n = n;
			req->k = k;
			req->f = d;			
			if((1<<w)<n||req->f+1>n)
			{
				printf("invalid w values.\n");
				return(-1);
			}
			req->type = type;			
			req->multiple_of = k*req->f*ALIGNMENT*w;	
			req->min_size = req->multiple_of;
			req->max_size = MIN(req->multiple_of*1024*1024*8,MAXPACKETSIZE);		
			req->d = MIN(2*req->f,n-1);
			req->w = w;
			break;	
		case LRC:	
			req->n = n;
			req->k = k;
			req->f = d;			
			if((1<<w)<n||req->f+1>n)
			{
				printf("invalid w values.\n");
				return(-1);
			}
			req->type = type;			
			req->multiple_of = k*(req->f)*ALIGNMENT*w;	
			req->min_size = req->multiple_of;
			req->max_size = MIN(req->multiple_of*1024*1024*8,MAXPACKETSIZE);		
			req->d = d;
			req->w = w;
			break;		
		case STEINERCODE:
		default: 
			printf("This type of regenerating code is not supported. \n");
			return(-1);			
	}
		
	return(req->multiple_of);
}

int compute_coded_packet_size(struct requirement *req, int data_size)
{
	int packet_size=-1;		
	switch (req->type)
	{
		case MBR_REPAIRBYTRANSFER: 
			packet_size = data_size/req->inner_k*(req->n-1);	         
			break;					
		case MSR_PRODUCTMATRIX:
			packet_size = data_size/(req->k);
			break;
		case MBR_PRODUCTMATRIX:
			packet_size = data_size/((req->k+1)*req->k/2+req->k*(req->d-req->k))*req->d;			
			break;
		case SRC:
			packet_size = data_size/(req->k*req->f)*(req->f+1);
			break;
		case LRC:
			packet_size = data_size/(req->k*req->f)*(req->f+1);
			break;			
		case STEINERCODE:
		default: 
			printf("This type of regenerating code is not supported. \n");
			return(-1);			
	}	
	return(packet_size);
}

int compute_repair_packet_size(struct requirement *req, int data_size)
{
	int packet_size=-1;		
	switch (req->type)
	{
		case MBR_REPAIRBYTRANSFER: 
			packet_size = data_size/req->inner_k;	         
			break;					
		case MSR_PRODUCTMATRIX:
			packet_size = data_size/(req->k*(req->d-req->k+1));
			break;
		case MBR_PRODUCTMATRIX:
			packet_size = data_size/((req->k+1)*req->k/2+req->k*(req->d-req->k));
			break;
		case SRC:
			packet_size = data_size/(req->k*(req->f))*(req->f+1);
			break;
		case LRC:
			packet_size = data_size/(req->k*(req->f))*(req->f+1);
			break;
		case STEINERCODE:
		default: 
			printf("This type of regenerating code is not supported. \n");
			return(-1);			
	}	
	return(packet_size);
}
/*
int checkValid(size_t inputsize, size_t outputsize, int packetsize, struct Requirement *requirement, int w)
{

	size_t wordsize = w/8;		
	if(inputsize/(wordsize*req->multipleof)!=0||inputsize>req->maxSize*wordsize||inputsize<req->maxSize*wordsize)
	{
		printf("In correct input packet size. \n");
		return(-1);
	}
	return(1);
	
}*/

int encode_rc(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info)
{
	switch (info->req.type)
	{
		case MBR_REPAIRBYTRANSFER: 			
			encode_MBR_repair_by_transfer(input, input_size, output, output_size, info);
			break;						
		case MSR_PRODUCTMATRIX:
			encode_MSR_product_matrix(input, input_size, output, output_size, info);
			break;						
		case MBR_PRODUCTMATRIX:
			encode_MBR_product_matrix(input, input_size, output, output_size, info);
			break;								
		case SRC:
			encode_SRC(input, input_size, output, output_size, info);
			break;
		case LRC:
			encode_LRC(input, input_size, output, output_size, info);
			break;
		case STEINERCODE:
		default: 
			printf("This type of regenerating code is not supported. \n");
			return(-1);			
	}
	return(1);	
}
int decode_rc(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info)
{
	switch (info->req.type)
	{
		case MBR_REPAIRBYTRANSFER: 			
			decode_MBR_repair_by_transfer(input, input_size, output, output_size, erasures, info);
			break;						
		case MSR_PRODUCTMATRIX:
			decode_MSR_product_matrix(input, input_size, output, output_size, erasures, info);
			break;						
		case MBR_PRODUCTMATRIX:
			decode_MBR_product_matrix(input, input_size, output, output_size, erasures, info);
			break;								
		case SRC:
			decode_SRC(input, input_size, output, output_size, erasures, info);
			break;
		case LRC:
			decode_LRC(input, input_size, output, output_size, erasures, info);
			break;
		case STEINERCODE:
		default: 
			printf("This type of regenerating code is not supported. \n");
			return(-1);			
	}
	return(1);	
}
int repair_encode_rc(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info)
{
	switch (info->req.type)
	{
		case MBR_REPAIRBYTRANSFER: 			
			repair_encode_MBR_repair_by_transfer(input, input_size, output, output_size, from_device_ID, to_device_ID, info);
			break;						
		case MSR_PRODUCTMATRIX:
			repair_encode_MSR_product_matrix(input, input_size, output, output_size, from_device_ID, to_device_ID, info);
			break;						
		case MBR_PRODUCTMATRIX:
			repair_encode_MBR_product_matrix(input, input_size, output, output_size, from_device_ID, to_device_ID, info);
			break;								
		case SRC:
			repair_encode_SRC(input, input_size, output, output_size, from_device_ID, to_device_ID, info);
			break;
		case LRC:
			repair_encode_LRC(input, input_size, output, output_size, from_device_ID, to_device_ID, info);
			break;
		case STEINERCODE:
		default: 
			printf("This type of regenerating code is not supported. \n");
			return(-1);			
	}
	return(1);	
}
int repair_decode_rc(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info)
{
	switch (info->req.type)
	{
		case MBR_REPAIRBYTRANSFER: 			
			repair_decode_MBR_repair_by_transfer(input, input_size, output, output_size, to_device_ID, helpers, info);
			break;						
		case MSR_PRODUCTMATRIX:
			repair_decode_MSR_product_matrix(input, input_size, output, output_size, to_device_ID, helpers, info);
			break;						
		case MBR_PRODUCTMATRIX:
			repair_decode_MBR_product_matrix(input, input_size, output, output_size, to_device_ID, helpers, info);
			break;								
		case SRC:
			repair_decode_SRC(input, input_size, output, output_size, to_device_ID, helpers, info);
			break;
		case LRC:
			repair_decode_LRC(input, input_size, output, output_size, to_device_ID, helpers, info);
			break;
		case STEINERCODE:
		default: 
			printf("This type of regenerating code is not supported. \n");
			return(-1);			
	}
	return(1);	
}


