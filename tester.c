/*
# RegeneratingCodes/tester.c

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
#include <time.h>
#include "jerasure.h"
#include "cauchy.h"
#include "regenerating_codes.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))
#define NUM_REPEAT 50

//#define DEBUG
#ifdef DEBUG
static void print_data(int start,int length, char* data)
{
	int i;
	for(i=start;i<start+length;i++)
		printf("%hhx ",*(data+i));
	printf("\n");

}

static void print_data_and_coding(int n, int size, void **data_plus_coding) 
{
	int i, j;
	for(i=0;i<n;i++)
	{
		printf("device %d\n",i);
		for(j=0;j<size;j++)		
			printf("%hhx ",*((char*)(data_plus_coding[i]+j)));		
		printf("\n");
	}
}
#endif /*BEBUG*/

// global pointers to various memory blocks for simplicity
char *data=NULL,*decoded_data=NULL, *repaired=NULL;  
char **coded=NULL, **repair_data=NULL;  
int *erasures=NULL, *erased=NULL;

int alloc_all_buffers(int size_of_data, int coded_packet_size, int repair_packet_size, struct coding_info *info)
{
	int n = info->req.n;
	int k = info->req.k;
	int d = info->req.d;
	int i;   
	long l;
	data = malloc(size_of_data); 
	repaired = malloc(coded_packet_size);
	decoded_data = malloc(size_of_data);  	
	coded = malloc(sizeof(char*)*n);
	repair_data = malloc(sizeof(char*)*(n-1));  
	if(coded==NULL||repair_data==NULL){
		printf("Out of memory!\n"); return(-1);
	}
	for (i = 0; i < n; i++) {
		coded[i] = malloc(coded_packet_size);
	}
	for (i = 0; i < d; i++)
		repair_data[i] = malloc(repair_packet_size);
	erasures = talloc(int, (n+1));
	erased = talloc(int, (n));
	
	// initialization
	for(i=0 ; i<size_of_data/sizeof(long); i++){
		l = lrand48(); l = (l << 32) | lrand48();
		memcpy((long*)(data+i*sizeof(long)), &l, sizeof(long));
	}
	for (i = 0; i < n; i++) 
		erased[i] = 0;
	for (i = 0; i < n-k; ){
		erasures[i] = lrand48()%(n);
		if (erased[erasures[i]] == 0) {
			erased[erasures[i]] = 1;
			memset(coded[erasures[i]],0,coded_packet_size);
			i++;
		}
	}
	erasures[i] = -1;
	return(1);
}

void free_all_buffers(struct coding_info *info)
{
	int i;
	int n = info->req.n;
	int d = info->req.d;
	if(data!=NULL)
		free(data);  
	if(decoded_data!=NULL)
		free(decoded_data);
	if(repaired!=NULL)
		free(repaired);
	if(erased!=NULL)
		free(erased);
	if(erasures!=NULL)
		free(erasures);

	if(coded!=NULL){  
		for(i=0;i<n;i++)
			free(coded[i]);
		free(coded);
	}
	if(repair_data!=NULL){
		for(i=0;i<d;i++)
			free(repair_data[i]);
		free(repair_data);
	}
}

void usage(char *s)
{
	printf("Usage: tester type_number size_of_data n k w v\n");
	printf("  type_number: \n");
	printf("	0: local regenerating code\n");
	printf("	1: simple regenerating codes\n");
	printf("	2: MSR codes using product matrix\n");
	printf("	3: MBR codes using product matrix\n");
	printf("	4: MBR codes with repair by transfer\n");
	printf("  parameters: \n");
	printf("	n: number of devices.\n");
	printf("	k: number of information devices\n");	
	printf("	w: number of bits per word (alphabet size) is the range of 1 and 32\n");
	printf("	v: type=0 or 1-->value of f; type=2 or 3-->value of d; type=4 or 5-->null.\n");
	
	if (s != NULL) fprintf(stderr, "\n Error: %s\n", s);
	exit(1);
}


int main(int argc, char **argv)
{
	int k, w, i, j, n, v, base, counter;  
	struct coding_info info;
	enum codetype type;
	int* helpers;
	int size_of_data;
	clock_t enc_clk=0, dec_clk=0, rep_enc_clk=0, rep_dec_clk=0, clk;

	srand48((unsigned)time(NULL) ); 
	if(argc<6)
		usage(NULL);
	
	
	switch(atoi(argv[1])){
		case 0: 	
			type = LRC;
			if(argc<6)
				usage(NULL);
			break;
		case 1: 	
			type = SRC;
			if(argc<6)
				usage(NULL);
			break;
		case 2: 	
			type = MSR_PRODUCTMATRIX;
			if(argc<6)
				usage(NULL);
			break;
		case 3: 	
			type = MBR_PRODUCTMATRIX;
			if(argc<6)
				usage(NULL);
			break;
		case 4: 	
			type = MBR_REPAIRBYTRANSFER;
			break;
		default: usage("unrecognized code type.");		
	}
	if (sscanf(argv[2], "%d", &size_of_data) == 0 || size_of_data <= 0)
		usage(NULL);	

	if (sscanf(argv[3], "%d", &n) == 0 || n <= 0)
		usage(NULL);

	if (sscanf(argv[4], "%d", &k) == 0 || k <= 0)
		usage(NULL);
	
	if (sscanf(argv[5], "%d", &w) == 0 || w <= 0 || w > 32)
		usage(NULL);
	v = n-1;	
	
	if(type == LRC||type == SRC){
		if (sscanf(argv[6], "%d", &v) == 0 || v <= 0|| (type == LRC&&n%(v+1)>0))
			usage(NULL);
	}	
	else if(type == MSR_PRODUCTMATRIX||type == MBR_PRODUCTMATRIX){	
		if (sscanf(argv[6], "%d", &v) == 0 || v <= 0)
			usage(NULL);
	}

	
	if (w < 30 && (n) > (1 << w)){
		usage("w is too small\n");
	}
	// setup parameters
	if(get_requirement(type, &(info.req), n, k, v, w)<0){
		printf("can not get coding requirements. Check parameters.\n");
		exit(1);
	}
	int coded_packet_size, repair_packet_size;	
	make_coding_matrics(&info);
	size_of_data = (int)(size_of_data/info.req.multiple_of)*info.req.multiple_of;
	coded_packet_size = compute_coded_packet_size(&(info.req),size_of_data);
	repair_packet_size = compute_repair_packet_size(&info.req,size_of_data);

	int erased_ID;
	int repeat_count;

	for(repeat_count=0; repeat_count< NUM_REPEAT ; repeat_count++)
	{
		erased_ID = lrand48()%n;	
	
		if(coded_packet_size<0||repair_packet_size<0){
			printf("Can not compute coded or repair packet size\n"); exit(-1);
		}  
		//setup data buffer and generate random data
		if(alloc_all_buffers(size_of_data,coded_packet_size,repair_packet_size,&info)<0)
			goto complete;

		clk = clock();
		//start testing: encode+decode
		if(encode_rc(data, size_of_data, coded, coded_packet_size, &info)<0){
		 	printf("Failed to encode"); exit(1);
		} 
		enc_clk += clock()-clk;
		//print_data_and_coding(n, codedPacketSize, coded);
		clk = clock();
		if(decode_rc(coded, coded_packet_size, decoded_data, size_of_data, erasures, &info)<0){
			printf("Failed to encode"); goto complete;
		}  	
		dec_clk += clock()-clk;
		if(memcmp(data,decoded_data, size_of_data))
			printf("Incorrected decoded.\n");
		//else
		//	printf("Complete testing encoding/decoding with no error.\n");
		
		clk = clock();
		switch(type){
			case LRC: 	
				helpers = erasures; // reuse erasures buffer for simplicity
				base = (erased_ID/(v+1))*(v+1);
				for(j=base, counter=0; j<base+v+1; j++){
					if(j!=erased_ID){
						helpers[counter] = j; 			  
						if(repair_encode_rc(coded[j], coded_packet_size, repair_data[counter], repair_packet_size, j, erased_ID, &info)<0){
							printf("Can not generate repair data"); goto complete;
						}
						counter++;
					}
				}
				helpers[info.req.f] = -1;
				break;
			case SRC: 		
				helpers = erasures; // reuse erasures buffer for simplicity
				for(i=0,counter=0; i<n;i++){
					if(i==erased_ID)
						continue;
					if(i<erased_ID&&(erased_ID-i>info.req.f&&i+n-erased_ID>info.req.f))
						continue;
					if(i>erased_ID&&(i-erased_ID>info.req.f&&erased_ID+n-i>info.req.f))
						continue;
					helpers[counter] = i;
					if(repair_encode_rc(coded[i], coded_packet_size, repair_data[counter], repair_packet_size, i, erased_ID, &info)<0){
						printf("Can not generate repair data"); goto complete;
					}
					counter++;
				}	
				helpers[info.req.d] = -1;			
				break;
			case MSR_PRODUCTMATRIX: 	
				//for (i = 0; i < n; i++) 
				//	erased[i] = 0;
				memset(erased,0,sizeof(int)*n);
				erased[erased_ID] = 1;
				for (i = 0; i < n-info.req.d-1; ){
					erasures[i] = lrand48()%(n);
					if (erased[erasures[i]] == 0) {
						erased[erasures[i]] = 1; i++;
						}
				}
				helpers = erasures;  // reuse the erasures buffer, but rename it for better clarity
				for(j=0,i=0;i<n&&j<info.req.d;i++){
					if(erased[i]==0){
						helpers[j] = i; j++;
					}		
				}
				helpers[j] = -1;		
		
				for(i=0;i<info.req.d;i++){        
					if(repair_encode_rc(coded[helpers[i]], coded_packet_size, repair_data[i], repair_packet_size, helpers[i], erased_ID, &info)<0){
						printf("Can not generate repair data"); goto complete;
					}
				}  
				break;
			case MBR_PRODUCTMATRIX: 	
				/*for (i = 0; i < n; i++) 
					erased[i] = 0;*/
				memset(erased,0,sizeof(int)*n);
				erased[erased_ID] = 1;
				for (i = 0; i < n-info.req.d-1; ){
					erasures[i] = lrand48()%(n);
					if (erased[erasures[i]] == 0) {
						erased[erasures[i]] = 1; i++;
					}
				}
				helpers = erasures;  // reuse the erasures buffer, but rename it for better clarity
				for(j=0,i=0;i<n&&j<info.req.d;i++){
					if(erased[i]==0){
						helpers[j] = i; j++;
					}		
				}		
				helpers[j] = -1;
	
				for(i=0;i<info.req.d;i++){        
					if(repair_encode_rc(coded[helpers[i]], coded_packet_size, repair_data[i], repair_packet_size, helpers[i], erased_ID, &info)<0){
						printf("Can not generate repair data"); goto complete;
					}
				}  
				break;
			case MBR_REPAIRBYTRANSFER: 	
				helpers = erasures; // reuse erasure buffer for simplicity
				for(j=0, i=(erased_ID+1)%n ;i!=erased_ID;i=(i+1)%n,j++){
					helpers[j] = i;
					if(repair_encode_rc(coded[i], coded_packet_size, repair_data[j], repair_packet_size, i, erased_ID, &info)<0){
						printf("Can not generate repair data"); goto complete;
					}
				}
				break;
			default: 
				usage("unrecognized code type.");		
				exit(1);
		}	
		rep_enc_clk += clock()-clk;
		clk = clock();
		if(repair_decode_rc(repair_data, repair_packet_size, repaired, coded_packet_size, erased_ID, helpers, &info)<0){
			printf("Can not generate repair data"); goto complete;
		}  
		rep_dec_clk += clock()-clk;
		if(memcmp(coded[erased_ID], repaired, coded_packet_size))
			printf("Incorrected repaired.\n");
		//else
		//	printf("Complete testing repair with no error.\n");
	
	complete:		
		free_all_buffers(&info);		
	}
	cleanup_matrics(&info);
/*
	printf("Running codec=%d, n=%d, k=%d, w=%d, v=%d\n", atoi(argv[1]),n,k,w,v);
	printf("encode average:      "
            " %.3e clocks %.3e sec %.3e bytes/sec\n",
            (double)enc_clk, (double)enc_clk/CLOCKS_PER_SEC,
            size_of_data * NUM_REPEAT* (double)CLOCKS_PER_SEC /
                ((double)enc_clk > 0 ? (double)enc_clk : 1e-10));
	printf("decode average:      "
            " %.3e clocks %.3e sec %.3e bytes/sec\n",
            (double)dec_clk, (double)dec_clk/CLOCKS_PER_SEC,
            size_of_data* NUM_REPEAT* (double)CLOCKS_PER_SEC /
                ((double)dec_clk > 0 ? (double)dec_clk : 1e-10));
	printf("repair encode average:      "
            " %.3e clocks %.3e sec %.3e bytes/sec\n",
            (double)rep_enc_clk, (double)rep_enc_clk/CLOCKS_PER_SEC,
            coded_packet_size * NUM_REPEAT* (double)CLOCKS_PER_SEC /
                ((double)rep_enc_clk > 0 ? (double)rep_enc_clk : 1e-10));
	printf("repair decode average:      "
            " %.3e clocks %.3e sec %.3e bytes/sec\n",
            (double)rep_dec_clk, (double)rep_dec_clk/CLOCKS_PER_SEC,
            coded_packet_size* NUM_REPEAT* (double)CLOCKS_PER_SEC /
                ((double)rep_dec_clk > 0 ? (double)rep_dec_clk : 1e-10));*/
	printf("%.3e, %.3e, %.3e, %.3e;\n", 
		size_of_data * NUM_REPEAT* (double)CLOCKS_PER_SEC /((double)enc_clk > 0 ? (double)enc_clk : 1e-10),  
		size_of_data * NUM_REPEAT* (double)CLOCKS_PER_SEC /((double)dec_clk > 0 ? (double)dec_clk : 1e-10),
		coded_packet_size * NUM_REPEAT* (double)CLOCKS_PER_SEC /((double)rep_enc_clk > 0 ? (double)rep_enc_clk : 1e-10),
		coded_packet_size* NUM_REPEAT* (double)CLOCKS_PER_SEC /((double)rep_dec_clk > 0 ? (double)rep_dec_clk : 1e-10));
	return 0;
	
}	
