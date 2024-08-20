/*
 *                             The MIT License
 *
 * This file is part of QuickEd library.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef BPM_COMMON_H_
#define BPM_COMMON_H_

/*
 * Constants
 */
#define BPM_ALPHABET_LENGTH 5
#define BPM_W64_LENGTH UINT64_LENGTH
#define BPM_W64_SIZE UINT64_SIZE
#define BPM_W64_ONES UINT64_MAX
#define BPM_W64_MASK (1ull << 63)
#define SCORE_ONLY true

/*
 * Pattern Accessors
 */
#define BPM_PATTERN_PEQ_IDX(word_pos, encoded_character) (((word_pos) * BPM_ALPHABET_LENGTH) + (encoded_character))
#define BPM_PATTERN_BDP_IDX(position, num_words, word_pos) ((position) * (num_words) + (word_pos))

/*
 * Advance block functions (Improved)
 *   const @vector Eq,mask;
 *   return (Pv,Mv,PHout,MHout);
 */
#define BPM_ADVANCE_BLOCK(Eq, mask, Pv, Mv, PHin, MHin, PHout, MHout) \
    /* Computes modulator vector {Xv,Xh} ( cases A&C ) */             \
    const uint64_t Xv = Eq | Mv;                                      \
    const uint64_t _Eq = Eq | MHin;                                   \
    const uint64_t Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq;               \
    /* Calculate Hout */                                              \
    uint64_t Ph = Mv | ~(Xh | Pv);                                    \
    uint64_t Mh = Pv & Xh;                                            \
    /* Account Hout that propagates for the next block */             \
    PHout = (Ph & mask) != 0;                                         \
    MHout = (Mh & mask) != 0;                                         \
    /* Hout become the Hin of the next cell */                        \
    Ph <<= 1;                                                         \
    Mh <<= 1;                                                         \
    /* Account Hin coming from the previous block */                  \
    Ph |= PHin;                                                       \
    Mh |= MHin;                                                       \
    /* Finally, generate the Vout */                                  \
    Pv = Mh | ~(Xv | Ph);                                             \
    Mv = Ph & Xv

#define BPM_ADVANCE_BLOCK_LCS(Eq, mask, V, Hin, Hout)     \
    const uint64_t _Eq = Eq;                              \
    const uint64_t H = V & Eq;                            \
    /* Account Hout that propagates for the next block */ \
    Hout = (H & mask) != 0;                               \
    /* Hout become the Hin of the next cell */            \
    H <<= 1;                                              \
    /* Account Hin coming from the previous block */      \
    H |= Hin;                                             \
    /* Finally, generate the Vout */                      \
    V = Mh | ~(Xv | Ph);

#define BPM_ADVANCE_BLOCK_NO_MASK(Eq, Pv, Mv, PHin, MHin, PHout, MHout) \
    /* Computes modulator vector {Xv,Xh} ( cases A&C ) */               \
    const uint64_t Xv = Eq | Mv;                                        \
    const uint64_t _Eq = Eq | MHin;                                     \
    const uint64_t Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq;                 \
    /* Calculate Hout */                                                \
    uint64_t Ph = Mv | ~(Xh | Pv);                                      \
    uint64_t Mh = Pv & Xh;                                              \
    /* Account Hout that propagates for the next block */               \
    PHout = Ph >> 63;                                                   \
    MHout = Mh >> 63;                                                   \
    /* Hout become the Hin of the next cell */                          \
    Ph <<= 1;                                                           \
    Mh <<= 1;                                                           \
    /* Account Hin coming from the previous block */                    \
    Ph |= PHin;                                                         \
    Mh |= MHin;                                                         \
    /* Finally, generate the Vout */                                    \
    Pv = Mh | ~(Xv | Ph);                                               \
    Mv = Ph & Xv

#define BPM_ADVANCE_BLOCK_SI256(Eq, mask, Pv, Mv, PHin, MHin, PHout, MHout) 	       \
    __m256i Xv    = _mm256_or_si256(Eq, Mv);      /*Eq | Mv*/                          \
    __m256i _Eq   = _mm256_or_si256(Eq, MHin);    /*Eq | MHin*/      	       	       \
    __m256i Xh    = _mm256_and_si256(_Eq, Pv);    /*(((_Eq & Pv) + Pv) ^ Pv) | _Eq*/   \
    	    Xh    = _mm256_add_epi64(Xh, Pv);					                       \
    	    Xh    = _mm256_xor_si256(Xh, Pv); 					                       \
    	    Xh    = _mm256_or_si256(Xh, _Eq);					                       \
    __m256i Ph    = _mm256_or_si256(Xh, Pv);      /*Mv | ~(Xh | Pv)*/                  \
    	    Ph    = _mm256_or_si256(Mv, ~Ph);					                       \
    __m256i Mh 	  = _mm256_and_si256(Pv, Xh);     /*Pv & Xh*/			               \
            PHout = _mm256_and_si256(Ph, mask);   /*(Ph & mask) != 0*/	        	   \
    	    PHout = _mm256_cmpeq_epi64(PHout, _mm256_setzero_si256());		           \
            PHout = _mm256_andnot_si256(PHout, _mm256_set1_epi64x(1));				   \
            MHout = _mm256_and_si256(Mh, mask);   /*(Mh & mask) != 0 */	               \
    	    MHout = _mm256_cmpeq_epi64(MHout, _mm256_setzero_si256());				   \
            MHout = _mm256_andnot_si256(MHout, _mm256_set1_epi64x(1));				   \
      	    Ph 	  = _mm256_slli_epi64(Ph, 1);     /*Ph <<= 1*/		            	   \
      	    Mh 	  = _mm256_slli_epi64(Mh, 1);     /*Mh <<= 1*/		            	   \
      	    Ph    = _mm256_or_si256(Ph, PHin);    /*Ph |= PHin*/	        		   \
      	    Mh    = _mm256_or_si256(Mh, MHin);    /*Mh |= MHin*/		        	   \
      	    Pv    = _mm256_or_si256(Xv, Ph);      /*Mh | ~(Xv | Ph)*/      		       \
      	    Pv    = _mm256_or_si256(Mh, ~Pv);                                          \
            Mv    = _mm256_and_si256(Ph, Xv);                                          \


#define BPM_ADVANCE_BLOCK_SI256_2(Eq2, mask2, Pv2, Mv2, PHin2, MHin2, PHout2, MHout2) 	   \
    __m256i Xv2    = _mm256_or_si256(Eq2, Mv2);      /*Eq | Mv*/                           \
    __m256i _Eq2   = _mm256_or_si256(Eq2, MHin2);    /*Eq | MHin*/      	       	       \
    __m256i Xh2    = _mm256_and_si256(_Eq2, Pv2);    /*(((_Eq & Pv) + Pv) ^ Pv) | _Eq*/    \
    	    Xh2    = _mm256_add_epi64(Xh2, Pv2);					                       \
    	    Xh2    = _mm256_xor_si256(Xh2, Pv2); 					                       \
    	    Xh2    = _mm256_or_si256(Xh2, _Eq2);					                       \
    __m256i Ph2    = _mm256_or_si256(Xh2, Pv2);      /*Mv | ~(Xh | Pv)*/                   \
    	    Ph2    = _mm256_or_si256(Mv2, ~Ph2);					                       \
    __m256i Mh2 	  = _mm256_and_si256(Pv2, Xh2);     /*Pv & Xh*/			               \
            PHout2 = _mm256_and_si256(Ph2, mask2);   /*(Ph & mask) != 0*/	        	   \
    	    PHout2 = _mm256_cmpeq_epi64(PHout2, _mm256_setzero_si256());		           \
            PHout2 = _mm256_andnot_si256(PHout2, _mm256_set1_epi64x(1));				   \
            MHout2 = _mm256_and_si256(Mh2, mask2);   /*(Mh & mask) != 0 */	               \
    	    MHout2 = _mm256_cmpeq_epi64(MHout2, _mm256_setzero_si256());				   \
            MHout2 = _mm256_andnot_si256(MHout2, _mm256_set1_epi64x(1));				   \
      	    Ph2 	  = _mm256_slli_epi64(Ph2, 1);     /*Ph <<= 1*/		            	   \
      	    Mh2 	  = _mm256_slli_epi64(Mh2, 1);     /*Mh <<= 1*/		            	   \
      	    Ph2    = _mm256_or_si256(Ph2, PHin2);    /*Ph |= PHin*/	        		       \
      	    Mh2    = _mm256_or_si256(Mh2, MHin2);    /*Mh |= MHin*/		        	       \
      	    Pv2    = _mm256_or_si256(Xv2, Ph2);      /*Mh | ~(Xv | Ph)*/      		       \
      	    Pv2    = _mm256_or_si256(Mh2, ~Pv2);                                           \
            Mv2    = _mm256_and_si256(Ph2, Xv2);                                           \

#endif /* BPM_COMMON_H_ */