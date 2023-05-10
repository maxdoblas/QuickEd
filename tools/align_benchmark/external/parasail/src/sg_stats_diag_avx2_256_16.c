/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>

#include <immintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_avx.h"

#define SG_STATS
#define SG_SUFFIX _diag_avx2_256_16
#include "sg_helper.h"

#define NEG_INF (INT16_MIN/(int16_t)(2))

#define _mm256_cmplt_epi16_rpl(a,b) _mm256_cmpgt_epi16(b,a)

#if HAVE_AVX2_MM256_INSERT_EPI16
#define _mm256_insert_epi16_rpl _mm256_insert_epi16
#else
static inline __m256i _mm256_insert_epi16_rpl(__m256i a, int16_t i, int imm) {
    __m256i_16_t A;
    A.m = a;
    A.v[imm] = i;
    return A.m;
}
#endif

#if HAVE_AVX2_MM256_EXTRACT_EPI16
#define _mm256_extract_epi16_rpl _mm256_extract_epi16
#else
static inline int16_t _mm256_extract_epi16_rpl(__m256i a, int imm) {
    __m256i_16_t A;
    A.m = a;
    return A.v[imm];
}
#endif

#define _mm256_srli_si256_rpl(a,imm) _mm256_or_si256(_mm256_slli_si256(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(3,0,0,1)), 16-imm), _mm256_srli_si256(a, imm))


#ifdef PARASAIL_TABLE
static inline void arr_store_si256(
        int *array,
        __m256i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (0 <= i+0 && i+0 < s1Len && 0 <= j-0 && j-0 < s2Len) {
        array[1LL*(i+0)*s2Len + (j-0)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 15);
    }
    if (0 <= i+1 && i+1 < s1Len && 0 <= j-1 && j-1 < s2Len) {
        array[1LL*(i+1)*s2Len + (j-1)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 14);
    }
    if (0 <= i+2 && i+2 < s1Len && 0 <= j-2 && j-2 < s2Len) {
        array[1LL*(i+2)*s2Len + (j-2)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 13);
    }
    if (0 <= i+3 && i+3 < s1Len && 0 <= j-3 && j-3 < s2Len) {
        array[1LL*(i+3)*s2Len + (j-3)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 12);
    }
    if (0 <= i+4 && i+4 < s1Len && 0 <= j-4 && j-4 < s2Len) {
        array[1LL*(i+4)*s2Len + (j-4)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 11);
    }
    if (0 <= i+5 && i+5 < s1Len && 0 <= j-5 && j-5 < s2Len) {
        array[1LL*(i+5)*s2Len + (j-5)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 10);
    }
    if (0 <= i+6 && i+6 < s1Len && 0 <= j-6 && j-6 < s2Len) {
        array[1LL*(i+6)*s2Len + (j-6)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 9);
    }
    if (0 <= i+7 && i+7 < s1Len && 0 <= j-7 && j-7 < s2Len) {
        array[1LL*(i+7)*s2Len + (j-7)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 8);
    }
    if (0 <= i+8 && i+8 < s1Len && 0 <= j-8 && j-8 < s2Len) {
        array[1LL*(i+8)*s2Len + (j-8)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 7);
    }
    if (0 <= i+9 && i+9 < s1Len && 0 <= j-9 && j-9 < s2Len) {
        array[1LL*(i+9)*s2Len + (j-9)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 6);
    }
    if (0 <= i+10 && i+10 < s1Len && 0 <= j-10 && j-10 < s2Len) {
        array[1LL*(i+10)*s2Len + (j-10)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 5);
    }
    if (0 <= i+11 && i+11 < s1Len && 0 <= j-11 && j-11 < s2Len) {
        array[1LL*(i+11)*s2Len + (j-11)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 4);
    }
    if (0 <= i+12 && i+12 < s1Len && 0 <= j-12 && j-12 < s2Len) {
        array[1LL*(i+12)*s2Len + (j-12)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 3);
    }
    if (0 <= i+13 && i+13 < s1Len && 0 <= j-13 && j-13 < s2Len) {
        array[1LL*(i+13)*s2Len + (j-13)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 2);
    }
    if (0 <= i+14 && i+14 < s1Len && 0 <= j-14 && j-14 < s2Len) {
        array[1LL*(i+14)*s2Len + (j-14)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 1);
    }
    if (0 <= i+15 && i+15 < s1Len && 0 <= j-15 && j-15 < s2Len) {
        array[1LL*(i+15)*s2Len + (j-15)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_rowcol(
        int *row,
        int *col,
        __m256i vWH,
        int32_t i,
        int32_t s1Len,
        int32_t j,
        int32_t s2Len)
{
    if (i+0 == s1Len-1 && 0 <= j-0 && j-0 < s2Len) {
        row[j-0] = (int16_t)_mm256_extract_epi16_rpl(vWH, 15);
    }
    if (j-0 == s2Len-1 && 0 <= i+0 && i+0 < s1Len) {
        col[(i+0)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 15);
    }
    if (i+1 == s1Len-1 && 0 <= j-1 && j-1 < s2Len) {
        row[j-1] = (int16_t)_mm256_extract_epi16_rpl(vWH, 14);
    }
    if (j-1 == s2Len-1 && 0 <= i+1 && i+1 < s1Len) {
        col[(i+1)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 14);
    }
    if (i+2 == s1Len-1 && 0 <= j-2 && j-2 < s2Len) {
        row[j-2] = (int16_t)_mm256_extract_epi16_rpl(vWH, 13);
    }
    if (j-2 == s2Len-1 && 0 <= i+2 && i+2 < s1Len) {
        col[(i+2)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 13);
    }
    if (i+3 == s1Len-1 && 0 <= j-3 && j-3 < s2Len) {
        row[j-3] = (int16_t)_mm256_extract_epi16_rpl(vWH, 12);
    }
    if (j-3 == s2Len-1 && 0 <= i+3 && i+3 < s1Len) {
        col[(i+3)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 12);
    }
    if (i+4 == s1Len-1 && 0 <= j-4 && j-4 < s2Len) {
        row[j-4] = (int16_t)_mm256_extract_epi16_rpl(vWH, 11);
    }
    if (j-4 == s2Len-1 && 0 <= i+4 && i+4 < s1Len) {
        col[(i+4)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 11);
    }
    if (i+5 == s1Len-1 && 0 <= j-5 && j-5 < s2Len) {
        row[j-5] = (int16_t)_mm256_extract_epi16_rpl(vWH, 10);
    }
    if (j-5 == s2Len-1 && 0 <= i+5 && i+5 < s1Len) {
        col[(i+5)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 10);
    }
    if (i+6 == s1Len-1 && 0 <= j-6 && j-6 < s2Len) {
        row[j-6] = (int16_t)_mm256_extract_epi16_rpl(vWH, 9);
    }
    if (j-6 == s2Len-1 && 0 <= i+6 && i+6 < s1Len) {
        col[(i+6)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 9);
    }
    if (i+7 == s1Len-1 && 0 <= j-7 && j-7 < s2Len) {
        row[j-7] = (int16_t)_mm256_extract_epi16_rpl(vWH, 8);
    }
    if (j-7 == s2Len-1 && 0 <= i+7 && i+7 < s1Len) {
        col[(i+7)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 8);
    }
    if (i+8 == s1Len-1 && 0 <= j-8 && j-8 < s2Len) {
        row[j-8] = (int16_t)_mm256_extract_epi16_rpl(vWH, 7);
    }
    if (j-8 == s2Len-1 && 0 <= i+8 && i+8 < s1Len) {
        col[(i+8)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 7);
    }
    if (i+9 == s1Len-1 && 0 <= j-9 && j-9 < s2Len) {
        row[j-9] = (int16_t)_mm256_extract_epi16_rpl(vWH, 6);
    }
    if (j-9 == s2Len-1 && 0 <= i+9 && i+9 < s1Len) {
        col[(i+9)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 6);
    }
    if (i+10 == s1Len-1 && 0 <= j-10 && j-10 < s2Len) {
        row[j-10] = (int16_t)_mm256_extract_epi16_rpl(vWH, 5);
    }
    if (j-10 == s2Len-1 && 0 <= i+10 && i+10 < s1Len) {
        col[(i+10)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 5);
    }
    if (i+11 == s1Len-1 && 0 <= j-11 && j-11 < s2Len) {
        row[j-11] = (int16_t)_mm256_extract_epi16_rpl(vWH, 4);
    }
    if (j-11 == s2Len-1 && 0 <= i+11 && i+11 < s1Len) {
        col[(i+11)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 4);
    }
    if (i+12 == s1Len-1 && 0 <= j-12 && j-12 < s2Len) {
        row[j-12] = (int16_t)_mm256_extract_epi16_rpl(vWH, 3);
    }
    if (j-12 == s2Len-1 && 0 <= i+12 && i+12 < s1Len) {
        col[(i+12)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 3);
    }
    if (i+13 == s1Len-1 && 0 <= j-13 && j-13 < s2Len) {
        row[j-13] = (int16_t)_mm256_extract_epi16_rpl(vWH, 2);
    }
    if (j-13 == s2Len-1 && 0 <= i+13 && i+13 < s1Len) {
        col[(i+13)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 2);
    }
    if (i+14 == s1Len-1 && 0 <= j-14 && j-14 < s2Len) {
        row[j-14] = (int16_t)_mm256_extract_epi16_rpl(vWH, 1);
    }
    if (j-14 == s2Len-1 && 0 <= i+14 && i+14 < s1Len) {
        col[(i+14)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 1);
    }
    if (i+15 == s1Len-1 && 0 <= j-15 && j-15 < s2Len) {
        row[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWH, 0);
    }
    if (j-15 == s2Len-1 && 0 <= i+15 && i+15 < s1Len) {
        col[(i+15)] = (int16_t)_mm256_extract_epi16_rpl(vWH, 0);
    }
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sg_flags_stats_table_diag_avx2_256_16
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sg_flags_stats_rowcol_diag_avx2_256_16
#else
#define FNAME parasail_sg_flags_stats_diag_avx2_256_16
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix,
        int s1_beg, int s1_end, int s2_beg, int s2_end)
{
    const int32_t N = 16; /* number of values in vector */
    const int32_t PAD = N-1;
    const int32_t PAD2 = PAD*2;
    const int32_t s1Len_PAD = s1Len+PAD;
    const int32_t s2Len_PAD = s2Len+PAD;
    int16_t * const restrict s1      = parasail_memalign_int16_t(32, s1Len+PAD);
    int16_t * const restrict s2B     = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _H_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _HM_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _HS_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _HL_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _F_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _FM_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _FS_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict _FL_pr = parasail_memalign_int16_t(32, s2Len+PAD2);
    int16_t * const restrict s2 = s2B+PAD; /* will allow later for negative indices */
    int16_t * const restrict H_pr = _H_pr+PAD;
    int16_t * const restrict HM_pr = _HM_pr+PAD;
    int16_t * const restrict HS_pr = _HS_pr+PAD;
    int16_t * const restrict HL_pr = _HL_pr+PAD;
    int16_t * const restrict F_pr = _F_pr+PAD;
    int16_t * const restrict FM_pr = _FM_pr+PAD;
    int16_t * const restrict FS_pr = _FS_pr+PAD;
    int16_t * const restrict FL_pr = _FL_pr+PAD;
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table3(s1Len, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol3(s1Len, s2Len);
#else
    parasail_result_t *result = parasail_result_new_stats();
#endif
#endif
    int32_t i = 0;
    int32_t j = 0;
    int32_t end_query = s1Len-1;
    int32_t end_ref = s2Len-1;
    const int16_t NEG_LIMIT = (-open < matrix->min ?
        INT16_MIN + open : INT16_MIN - matrix->min) + 1;
    const int16_t POS_LIMIT = INT16_MAX - matrix->max - 1;
    int16_t score = NEG_LIMIT;
    int16_t matches = NEG_LIMIT;
    int16_t similar = NEG_LIMIT;
    int16_t length = NEG_LIMIT;
    __m256i vNegLimit = _mm256_set1_epi16(NEG_LIMIT);
    __m256i vPosLimit = _mm256_set1_epi16(POS_LIMIT);
    __m256i vSaturationCheckMin = vPosLimit;
    __m256i vSaturationCheckMax = vNegLimit;
    __m256i vNegInf = _mm256_set1_epi16(NEG_LIMIT);
    __m256i vOpen = _mm256_set1_epi16(open);
    __m256i vGap  = _mm256_set1_epi16(gap);
    __m256i vZero = _mm256_set1_epi16(0);
    __m256i vOne = _mm256_set1_epi16(1);
    __m256i vN = _mm256_set1_epi16(N);
    __m256i vGapN = s1_beg ? _mm256_set1_epi16(0) : _mm256_set1_epi16(gap*N);
    __m256i vNegOne = _mm256_set1_epi16(-1);
    __m256i vI = _mm256_set_epi16(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
    __m256i vJreset = _mm256_set_epi16(0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15);
    __m256i vMaxHRow = vNegInf;
    __m256i vMaxMRow = vNegInf;
    __m256i vMaxSRow = vNegInf;
    __m256i vMaxLRow = vNegInf;
    __m256i vMaxHCol = vNegInf;
    __m256i vMaxMCol = vNegInf;
    __m256i vMaxSCol = vNegInf;
    __m256i vMaxLCol = vNegInf;
    __m256i vLastValH = vNegInf;
    __m256i vLastValM = vNegInf;
    __m256i vLastValS = vNegInf;
    __m256i vLastValL = vNegInf;
    __m256i vEndI = vNegInf;
    __m256i vEndJ = vNegInf;
    __m256i vILimit = _mm256_set1_epi16(s1Len);
    __m256i vILimit1 = _mm256_sub_epi16(vILimit, vOne);
    __m256i vJLimit = _mm256_set1_epi16(s2Len);
    __m256i vJLimit1 = _mm256_sub_epi16(vJLimit, vOne);
    __m256i vIBoundary = s1_beg ? _mm256_set1_epi16(0) : _mm256_set_epi16(
            -open-0*gap,
            -open-1*gap,
            -open-2*gap,
            -open-3*gap,
            -open-4*gap,
            -open-5*gap,
            -open-6*gap,
            -open-7*gap,
            -open-8*gap,
            -open-9*gap,
            -open-10*gap,
            -open-11*gap,
            -open-12*gap,
            -open-13*gap,
            -open-14*gap,
            -open-15*gap);

    /* convert _s1 from char to int in range 0-23 */
    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    /* pad back of s1 with dummy values */
    for (i=s1Len; i<s1Len_PAD; ++i) {
        s1[i] = 0; /* point to first matrix row because we don't care */
    }

    /* convert _s2 from char to int in range 0-23 */
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }
    /* pad front of s2 with dummy values */
    for (j=-PAD; j<0; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }
    /* pad back of s2 with dummy values */
    for (j=s2Len; j<s2Len_PAD; ++j) {
        s2[j] = 0; /* point to first matrix row because we don't care */
    }

    /* set initial values for stored row */
    if (s2_beg) {
        for (j=0; j<s2Len; ++j) {
            H_pr[j] = 0;
            HM_pr[j] = 0;
            HS_pr[j] = 0;
            HL_pr[j] = 0;
            F_pr[j] = NEG_INF;
            FM_pr[j] = 0;
            FS_pr[j] = 0;
            FL_pr[j] = 0;
        }
    }
    else {
        for (j=0; j<s2Len; ++j) {
            H_pr[j] = -open - j*gap;
            HM_pr[j] = 0;
            HS_pr[j] = 0;
            HL_pr[j] = 0;
            F_pr[j] = NEG_INF;
            FM_pr[j] = 0;
            FS_pr[j] = 0;
            FL_pr[j] = 0;
        }
    }
    /* pad front of stored row values */
    for (j=-PAD; j<0; ++j) {
        H_pr[j] = 0;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = 0;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    /* pad back of stored row values */
    for (j=s2Len; j<s2Len+PAD; ++j) {
        H_pr[j] = 0;
        HM_pr[j] = 0;
        HS_pr[j] = 0;
        HL_pr[j] = 0;
        F_pr[j] = 0;
        FM_pr[j] = 0;
        FS_pr[j] = 0;
        FL_pr[j] = 0;
    }
    H_pr[-1] = 0; /* upper left corner */

    /* iterate over query sequence */
    for (i=0; i<s1Len; i+=N) {
        __m256i case1 = vZero;
        __m256i case2 = vZero;
        __m256i vNH = vNegInf;
        __m256i vNM = vZero;
        __m256i vNS = vZero;
        __m256i vNL = vZero;
        __m256i vWH = vNegInf;
        __m256i vWM = vZero;
        __m256i vWS = vZero;
        __m256i vWL = vZero;
        __m256i vE = vNegInf;
        __m256i vE_opn = vNegInf;
        __m256i vE_ext = vNegInf;
        __m256i vEM = vZero;
        __m256i vES = vZero;
        __m256i vEL = vZero;
        __m256i vF = vNegInf;
        __m256i vF_opn = vNegInf;
        __m256i vF_ext = vNegInf;
        __m256i vFM = vZero;
        __m256i vFS = vZero;
        __m256i vFL = vZero;
        __m256i vJ = vJreset;
        __m256i vs1 = _mm256_set_epi16(
                s1[i+0],
                s1[i+1],
                s1[i+2],
                s1[i+3],
                s1[i+4],
                s1[i+5],
                s1[i+6],
                s1[i+7],
                s1[i+8],
                s1[i+9],
                s1[i+10],
                s1[i+11],
                s1[i+12],
                s1[i+13],
                s1[i+14],
                s1[i+15]);
        __m256i vs2 = vNegInf;
        const int * const restrict matrow0 = &matrix->matrix[matrix->size*s1[i+0]];
        const int * const restrict matrow1 = &matrix->matrix[matrix->size*s1[i+1]];
        const int * const restrict matrow2 = &matrix->matrix[matrix->size*s1[i+2]];
        const int * const restrict matrow3 = &matrix->matrix[matrix->size*s1[i+3]];
        const int * const restrict matrow4 = &matrix->matrix[matrix->size*s1[i+4]];
        const int * const restrict matrow5 = &matrix->matrix[matrix->size*s1[i+5]];
        const int * const restrict matrow6 = &matrix->matrix[matrix->size*s1[i+6]];
        const int * const restrict matrow7 = &matrix->matrix[matrix->size*s1[i+7]];
        const int * const restrict matrow8 = &matrix->matrix[matrix->size*s1[i+8]];
        const int * const restrict matrow9 = &matrix->matrix[matrix->size*s1[i+9]];
        const int * const restrict matrow10 = &matrix->matrix[matrix->size*s1[i+10]];
        const int * const restrict matrow11 = &matrix->matrix[matrix->size*s1[i+11]];
        const int * const restrict matrow12 = &matrix->matrix[matrix->size*s1[i+12]];
        const int * const restrict matrow13 = &matrix->matrix[matrix->size*s1[i+13]];
        const int * const restrict matrow14 = &matrix->matrix[matrix->size*s1[i+14]];
        const int * const restrict matrow15 = &matrix->matrix[matrix->size*s1[i+15]];
        __m256i vIltLimit = _mm256_cmplt_epi16_rpl(vI, vILimit);
        __m256i vIeqLimit1 = _mm256_cmpeq_epi16(vI, vILimit1);
        vNH = _mm256_srli_si256_rpl(vNH, 2);
        vNH = _mm256_insert_epi16_rpl(vNH, H_pr[-1], 15);
        vWH = _mm256_srli_si256_rpl(vWH, 2);
        vWH = _mm256_insert_epi16_rpl(vWH, s1_beg ? 0 : (-open - i*gap), 15);
        H_pr[-1] = -open - (i+N)*gap;
        /* iterate over database sequence */
        for (j=0; j<s2Len+PAD; ++j) {
            __m256i vMat;
            __m256i vNWH = vNH;
            __m256i vNWM = vNM;
            __m256i vNWS = vNS;
            __m256i vNWL = vNL;
            vNH = _mm256_srli_si256_rpl(vWH, 2);
            vNH = _mm256_insert_epi16_rpl(vNH, H_pr[j], 15);
            vNM = _mm256_srli_si256_rpl(vWM, 2);
            vNM = _mm256_insert_epi16_rpl(vNM, HM_pr[j], 15);
            vNS = _mm256_srli_si256_rpl(vWS, 2);
            vNS = _mm256_insert_epi16_rpl(vNS, HS_pr[j], 15);
            vNL = _mm256_srli_si256_rpl(vWL, 2);
            vNL = _mm256_insert_epi16_rpl(vNL, HL_pr[j], 15);
            vF = _mm256_srli_si256_rpl(vF, 2);
            vF = _mm256_insert_epi16_rpl(vF, F_pr[j], 15);
            vFM = _mm256_srli_si256_rpl(vFM, 2);
            vFM = _mm256_insert_epi16_rpl(vFM, FM_pr[j], 15);
            vFS = _mm256_srli_si256_rpl(vFS, 2);
            vFS = _mm256_insert_epi16_rpl(vFS, FS_pr[j], 15);
            vFL = _mm256_srli_si256_rpl(vFL, 2);
            vFL = _mm256_insert_epi16_rpl(vFL, FL_pr[j], 15);
            vF_opn = _mm256_sub_epi16(vNH, vOpen);
            vF_ext = _mm256_sub_epi16(vF, vGap);
            vF = _mm256_max_epi16(vF_opn, vF_ext);
            case1 = _mm256_cmpgt_epi16(vF_opn, vF_ext);
            vFM = _mm256_blendv_epi8(vFM, vNM, case1);
            vFS = _mm256_blendv_epi8(vFS, vNS, case1);
            vFL = _mm256_blendv_epi8(vFL, vNL, case1);
            vFL = _mm256_add_epi16(vFL, vOne);
            vE_opn = _mm256_sub_epi16(vWH, vOpen);
            vE_ext = _mm256_sub_epi16(vE, vGap);
            vE = _mm256_max_epi16(vE_opn, vE_ext);
            case1 = _mm256_cmpgt_epi16(vE_opn, vE_ext);
            vEM = _mm256_blendv_epi8(vEM, vWM, case1);
            vES = _mm256_blendv_epi8(vES, vWS, case1);
            vEL = _mm256_blendv_epi8(vEL, vWL, case1);
            vEL = _mm256_add_epi16(vEL, vOne);
            vs2 = _mm256_srli_si256_rpl(vs2, 2);
            vs2 = _mm256_insert_epi16_rpl(vs2, s2[j], 15);
            vMat = _mm256_set_epi16(
                    matrow0[s2[j-0]],
                    matrow1[s2[j-1]],
                    matrow2[s2[j-2]],
                    matrow3[s2[j-3]],
                    matrow4[s2[j-4]],
                    matrow5[s2[j-5]],
                    matrow6[s2[j-6]],
                    matrow7[s2[j-7]],
                    matrow8[s2[j-8]],
                    matrow9[s2[j-9]],
                    matrow10[s2[j-10]],
                    matrow11[s2[j-11]],
                    matrow12[s2[j-12]],
                    matrow13[s2[j-13]],
                    matrow14[s2[j-14]],
                    matrow15[s2[j-15]]
                    );
            vNWH = _mm256_add_epi16(vNWH, vMat);
            vWH = _mm256_max_epi16(vNWH, vE);
            vWH = _mm256_max_epi16(vWH, vF);
            case1 = _mm256_cmpeq_epi16(vWH, vNWH);
            case2 = _mm256_cmpeq_epi16(vWH, vF);
            vWM = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEM, vFM, case2),
                    _mm256_add_epi16(vNWM,
                        _mm256_and_si256(
                            _mm256_cmpeq_epi16(vs1,vs2),
                            vOne)),
                    case1);
            vWS = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vES, vFS, case2),
                    _mm256_add_epi16(vNWS,
                        _mm256_and_si256(
                            _mm256_cmpgt_epi16(vMat,vZero),
                            vOne)),
                    case1);
            vWL = _mm256_blendv_epi8(
                    _mm256_blendv_epi8(vEL, vFL, case2),
                    _mm256_add_epi16(vNWL, vOne), case1);
            /* as minor diagonal vector passes across the j=-1 boundary,
             * assign the appropriate boundary conditions */
            {
                __m256i cond = _mm256_cmpeq_epi16(vJ,vNegOne);
                vWH = _mm256_blendv_epi8(vWH, vIBoundary, cond);
                vWM = _mm256_andnot_si256(cond, vWM);
                vWS = _mm256_andnot_si256(cond, vWS);
                vWL = _mm256_andnot_si256(cond, vWL);
                vE = _mm256_blendv_epi8(vE, vNegInf, cond);
                vEM = _mm256_andnot_si256(cond, vEM);
                vES = _mm256_andnot_si256(cond, vES);
                vEL = _mm256_andnot_si256(cond, vEL);
            }
            vSaturationCheckMin = _mm256_min_epi16(vSaturationCheckMin, vWH);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vWH);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vWM);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vWS);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vWL);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vWL);
            vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vJ);
#ifdef PARASAIL_TABLE
            arr_store_si256(result->stats->tables->score_table, vWH, i, s1Len, j, s2Len);
            arr_store_si256(result->stats->tables->matches_table, vWM, i, s1Len, j, s2Len);
            arr_store_si256(result->stats->tables->similar_table, vWS, i, s1Len, j, s2Len);
            arr_store_si256(result->stats->tables->length_table, vWL, i, s1Len, j, s2Len);
#endif
#ifdef PARASAIL_ROWCOL
            arr_store_rowcol(result->stats->rowcols->score_row,   result->stats->rowcols->score_col, vWH, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->matches_row, result->stats->rowcols->matches_col, vWM, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->similar_row, result->stats->rowcols->similar_col, vWS, i, s1Len, j, s2Len);
            arr_store_rowcol(result->stats->rowcols->length_row,  result->stats->rowcols->length_col, vWL, i, s1Len, j, s2Len);
#endif
            H_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWH,0);
            HM_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWM,0);
            HS_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWS,0);
            HL_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vWL,0);
            F_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vF,0);
            FM_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vFM,0);
            FS_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vFS,0);
            FL_pr[j-15] = (int16_t)_mm256_extract_epi16_rpl(vFL,0);
            /* as minor diagonal vector passes across the i or j limit
             * boundary, extract the last value of the column or row */
            {
                __m256i vJeqLimit1 = _mm256_cmpeq_epi16(vJ, vJLimit1);
                __m256i vJgtNegOne = _mm256_cmpgt_epi16(vJ, vNegOne);
                __m256i vJltLimit = _mm256_cmplt_epi16_rpl(vJ, vJLimit);
                __m256i cond_j = _mm256_and_si256(vIltLimit, vJeqLimit1);
                __m256i cond_i = _mm256_and_si256(vIeqLimit1,
                        _mm256_and_si256(vJgtNegOne, vJltLimit));
                __m256i cond_max_row = _mm256_cmpgt_epi16(vWH, vMaxHRow);
                __m256i cond_max_col = _mm256_cmpgt_epi16(vWH, vMaxHCol);
                __m256i cond_last_val = _mm256_and_si256(vIeqLimit1, vJeqLimit1);
                __m256i cond_all_row = _mm256_and_si256(cond_max_row, cond_i);
                __m256i cond_all_col = _mm256_and_si256(cond_max_col, cond_j);
                vMaxHRow = _mm256_blendv_epi8(vMaxHRow, vWH, cond_all_row);
                vMaxMRow = _mm256_blendv_epi8(vMaxMRow, vWM, cond_all_row);
                vMaxSRow = _mm256_blendv_epi8(vMaxSRow, vWS, cond_all_row);
                vMaxLRow = _mm256_blendv_epi8(vMaxLRow, vWL, cond_all_row);
                vMaxHCol = _mm256_blendv_epi8(vMaxHCol, vWH, cond_all_col);
                vMaxMCol = _mm256_blendv_epi8(vMaxMCol, vWM, cond_all_col);
                vMaxSCol = _mm256_blendv_epi8(vMaxSCol, vWS, cond_all_col);
                vMaxLCol = _mm256_blendv_epi8(vMaxLCol, vWL, cond_all_col);
                vLastValH = _mm256_blendv_epi8(vLastValH, vWH, cond_last_val);
                vLastValM = _mm256_blendv_epi8(vLastValM, vWM, cond_last_val);
                vLastValS = _mm256_blendv_epi8(vLastValS, vWS, cond_last_val);
                vLastValL = _mm256_blendv_epi8(vLastValL, vWL, cond_last_val);
                vEndI = _mm256_blendv_epi8(vEndI, vI, cond_all_col);
                vEndJ = _mm256_blendv_epi8(vEndJ, vJ, cond_all_row);
            }
            vJ = _mm256_add_epi16(vJ, vOne);
        }
        vI = _mm256_add_epi16(vI, vN);
        vIBoundary = _mm256_sub_epi16(vIBoundary, vGapN);
        vSaturationCheckMax = _mm256_max_epi16(vSaturationCheckMax, vI);
    }

    /* alignment ending position */
    {
        int16_t max_rowh = NEG_INF;
        int16_t max_rowm = NEG_INF;
        int16_t max_rows = NEG_INF;
        int16_t max_rowl = NEG_INF;
        int16_t max_colh = NEG_INF;
        int16_t max_colm = NEG_INF;
        int16_t max_cols = NEG_INF;
        int16_t max_coll = NEG_INF;
        int16_t last_valh = NEG_INF;
        int16_t last_valm = NEG_INF;
        int16_t last_vals = NEG_INF;
        int16_t last_vall = NEG_INF;
        int16_t *rh = (int16_t*)&vMaxHRow;
        int16_t *rm = (int16_t*)&vMaxMRow;
        int16_t *rs = (int16_t*)&vMaxSRow;
        int16_t *rl = (int16_t*)&vMaxLRow;
        int16_t *ch = (int16_t*)&vMaxHCol;
        int16_t *cm = (int16_t*)&vMaxMCol;
        int16_t *cs = (int16_t*)&vMaxSCol;
        int16_t *cl = (int16_t*)&vMaxLCol;
        int16_t *lh = (int16_t*)&vLastValH;
        int16_t *lm = (int16_t*)&vLastValM;
        int16_t *ls = (int16_t*)&vLastValS;
        int16_t *ll = (int16_t*)&vLastValL;
        int16_t *i = (int16_t*)&vEndI;
        int16_t *j = (int16_t*)&vEndJ;
        int32_t k;
        for (k=0; k<N; ++k, ++rh, ++rm, ++rs, ++rl, ++ch, ++cm, ++cs, ++cl, ++lh, ++lm, ++ls, ++ll, ++i, ++j) {
            if (*ch > max_colh || (*ch == max_colh && *i < end_query)) {
                max_colh = *ch;
                end_query = *i;
                max_colm = *cm;
                max_cols = *cs;
                max_coll = *cl;
            }
            if (*rh > max_rowh) {
                max_rowh = *rh;
                end_ref = *j;
                max_rowm = *rm;
                max_rows = *rs;
                max_rowl = *rl;
            }
            if (*lh > last_valh) {
                last_valh = *lh;
                last_valm = *lm;
                last_vals = *ls;
                last_vall = *ll;
            }
        }
        if (s1_end && s2_end) {
            if (max_colh > max_rowh || (max_colh == max_rowh && end_ref == s2Len-1)) {
                score = max_colh;
                end_ref = s2Len-1;
                matches = max_colm;
                similar = max_cols;
                length = max_coll;
            }
            else {
                score = max_rowh;
                end_query = s1Len-1;
                matches = max_rowm;
                similar = max_rows;
                length = max_rowl;
            }
        }
        else if (s1_end) {
            score = max_colh;
            end_ref = s2Len-1;
            matches = max_colm;
            similar = max_cols;
            length = max_coll;
        }
        else if (s2_end) {
            score = max_rowh;
            end_query = s1Len-1;
            matches = max_rowm;
            similar = max_rows;
            length = max_rowl;
        }
        else {
            score = last_valh;
            end_query = s1Len-1;
            end_ref = s2Len-1;
            matches = last_valm;
            similar = last_vals;
            length = last_vall;
        }
    }

    if (_mm256_movemask_epi8(_mm256_or_si256(
            _mm256_cmplt_epi16_rpl(vSaturationCheckMin, vNegLimit),
            _mm256_cmpgt_epi16(vSaturationCheckMax, vPosLimit)))) {
        result->flag |= PARASAIL_FLAG_SATURATED;
        score = 0;
        matches = 0;
        similar = 0;
        length = 0;
        end_query = 0;
        end_ref = 0;
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;
    result->stats->matches = matches;
    result->stats->similar = similar;
    result->stats->length = length;
    result->flag |= PARASAIL_FLAG_SG | PARASAIL_FLAG_DIAG
        | PARASAIL_FLAG_STATS
        | PARASAIL_FLAG_BITS_16 | PARASAIL_FLAG_LANES_16;
    result->flag |= s1_beg ? PARASAIL_FLAG_SG_S1_BEG : 0;
    result->flag |= s1_end ? PARASAIL_FLAG_SG_S1_END : 0;
    result->flag |= s2_beg ? PARASAIL_FLAG_SG_S2_BEG : 0;
    result->flag |= s2_end ? PARASAIL_FLAG_SG_S2_END : 0;
#ifdef PARASAIL_TABLE
    result->flag |= PARASAIL_FLAG_TABLE;
#endif
#ifdef PARASAIL_ROWCOL
    result->flag |= PARASAIL_FLAG_ROWCOL;
#endif

    parasail_free(_FL_pr);
    parasail_free(_FS_pr);
    parasail_free(_FM_pr);
    parasail_free(_F_pr);
    parasail_free(_HL_pr);
    parasail_free(_HS_pr);
    parasail_free(_HM_pr);
    parasail_free(_H_pr);
    parasail_free(s2B);
    parasail_free(s1);

    return result;
}

SG_IMPL_ALL


