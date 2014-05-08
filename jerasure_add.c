/*
# RegeneratingCodes/jerasure_add.c

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
#include "jerasure.h"
#include "jerasure_add.h"

//added this function for cases where bitmatrix does not needs to be allocated each time

void jerasure_matrix_to_bitmatrix_noallocate(int k, int m, int w, int *matrix, int *bitmatrix) 
{

  int rowelts, rowindex, colindex, elt, i, j, l, x;
  
  if (matrix == NULL) { return; }

  rowelts = k * w;
  rowindex = 0;

  for (i = 0; i < m; i++) {
    colindex = rowindex;
    for (j = 0; j < k; j++) {
      elt = matrix[i*k+j];
      for (x = 0; x < w; x++) {
        for (l = 0; l < w; l++) {
          bitmatrix[colindex+x+l*rowelts] = ((elt & (1 << l)) ? 1 : 0);
        }
        elt = galois_single_multiply(elt, 2, w);
      }
      colindex += w;
    }
    rowindex += rowelts * w;
  }
}

