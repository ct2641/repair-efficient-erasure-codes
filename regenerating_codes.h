/*
# RegeneratingCodes/RegeneratingCodes.h
# header for regenerating code functions

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

#ifndef CODING_REGENERATING
#define CODING_REGENERATING

#include "galois.h"
#define MAXPACKETSIZE (67108864)
#define ALIGNMENT 512

enum codetype{
	MBR_REPAIRBYTRANSFER,
	MSR_PRODUCTMATRIX,
	MBR_PRODUCTMATRIX,
	SRC,
	LRC,
	STEINERCODE
};

struct requirement
{
	int min_size;
	int max_size;
	int multiple_of;
	int n,k,d,w;

	//extended fields
	int inner_n, inner_k;
	int f; //used by SRC and LRC

	enum codetype type;
};

struct coding_info
{
	struct requirement req;
	int* matrix;
	int* bitmatrix;
	int** schedule;
	// extended fields of submatrix for more sophisticated coding algorithms
	int num_of_submatrices;
	int** submatrix_array; 
	int** subbitmatrix_array;
	int*** subschedule_array;
};

	

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

int get_requirement(enum codetype type, struct requirement *req, int n, int k, int d, int w);
int compute_coded_packet_size(struct requirement *req, int data_size);
int compute_repair_packet_size(struct requirement *req, int data_size);
int make_coding_matrics(struct coding_info *info);
void cleanup_matrics(struct coding_info *info);

//int encode(void *input, size_t inputsize, void **output, size_t outputsize, struct CodingInfo *info);
//int CheckValid(size_t inputsize, size_t outputsize, int packetsize, struct PacketSizeRequirement *requirement, int w);

//specific coding functions, if you don't want to use the wrapper of Encode and Decode.

// MBR repair by transfer code
int encode_MBR_repair_by_transfer(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info);
int decode_MBR_repair_by_transfer(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info);
int repair_encode_MBR_repair_by_transfer(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info);
int repair_decode_MBR_repair_by_transfer(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info);

// Simple regenerating code
int encode_SRC(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info);
int decode_SRC(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info);
int repair_encode_SRC(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info);
int repair_decode_SRC(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info);

// Local regenerating code
int encode_LRC(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info);
int decode_LRC(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info);
int repair_encode_LRC(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info);
int repair_decode_LRC(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info);

// MBR code based on product matrix
int encode_MBR_product_matrix(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info);
int decode_MBR_product_matrix(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info);
int repair_encode_MBR_product_matrix(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info);
int repair_decode_MBR_product_matrix(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info);

// MSR code based on product matrix
int encode_MSR_product_matrix(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info);
int decode_MSR_product_matrix(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info);
int repair_encode_MSR_product_matrix(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info);
int repair_decode_MSR_product_matrix(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info);

int encode_rc(char *input, size_t input_size, char **output, size_t output_size, struct coding_info *info);
int decode_rc(char **input, size_t input_size, char *output, size_t output_size, int* erasures, struct coding_info *info);
int repair_encode_rc(char *input, size_t input_size, char *output, size_t output_size, int from_device_ID, int to_device_ID, struct coding_info *info);
int repair_decode_rc(char **input, size_t input_size, char *output, size_t output_size, int to_device_ID, int* helpers, struct coding_info *info);


#endif //CODING_REGENERATING
