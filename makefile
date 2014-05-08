# RegeneratingCodes/makefile
#Copyright (c) 2014, AT&T Intellectual Property.  All other rights reserved.

#Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#3. All advertising materials mentioning features or use of this software must display the following acknowledgement:  This product includes software developed by the AT&T.
#4. Neither the name of AT&T nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

#THIS SOFTWARE IS PROVIDED BY AT&T INTELLECTUAL PROPERTY ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL AT&T INTELLECTUAL PROPERTY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# Chao Tian
# AT&T Labs-Research
# Bedminster, NJ 07943
# tian@research.att.com

# $Revision: 0.1 $
# $Date: 2014/02/25 $

CC = gcc  
#CFLAGS = -g -Wall -I$(HOME)/include
CFLAGS = -O3 -mmmx -msse -DINTEL_SSE -msse2 -DINTEL_SSE2 -msse3 -DINTEL_SSE3 -mssse3 -msse4.1 -DINTEL_SSE4 -msse4.2 -DINTEL_SSE4 -fPIC -I$(HOME)/include -I./ -g -O2
INCLUDE = ./include
LIBDIR = ~/usr/local/lib
ALL =	tester

all: $(ALL)

clean:
	rm -f core *.o $(ALL) a.out 

install: $(ALL)
	rm *.o

.SUFFIXES: .c .o
.c.o:
	$(CC) $(CFLAGS) -c -I$(INCLUDE) $*.c

MBR_repair_by_transfer.o: regenerating_codes.h
SRC.o: regenerating_codes.h
LRC.o: regenerating_codes.h
MBR_product_matrix.o: regenerating_codes.h
MSR_product_matrix.o: regenerating_codes.h jerasure_add.h
regenerating_codes.o: regenerating_codes.h MSR_product_matrix.c MBR_product_matrix.c LRC.c SRC.c MBR_repair_by_transfer.c -lJerasure -lgf_complete
jerasure_add.o: jerasure_add.h

tester.o: regenerating_codes.h jerasure_add.h
tester: tester.o LRC.o SRC.o MBR_repair_by_transfer.o MBR_product_matrix.o MSR_product_matrix.o regenerating_codes.o jerasure_add.o
	$(CC) $(CFLAGS) -L$LIBDIR -o tester tester.o LRC.o regenerating_codes.o SRC.o MBR_repair_by_transfer.o MBR_product_matrix.o MSR_product_matrix.o jerasure_add.o -lJerasure -lgf_complete


