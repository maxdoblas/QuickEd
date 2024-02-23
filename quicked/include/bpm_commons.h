/*
 *                             The MIT License
 *
 * This file is part of QuickEdit library.
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

#endif /* BPM_COMMON_H_ */