
/* Create macros so that the matrices are stored in column-major order */

#define A(i, j) a[(j) *lda + (i)]
#define B(i, j) b[(j) *ldb + (i)]
#define C(i, j) c[(j) *ldc + (i)]

/* Block sizes */
#define mc 256
#define kc 128

#define min(i, j) ((i) < (j) ? (i) : (j))

/* Routine for computing C = A * B + C */

void AddDot4x4(int, double *, int, double *, int, double *, int);

void InnerKernel(int m, int n, int k, double *a, int lda,
                 double *b, int ldb,
                 double *c, int ldc);

/**
 * // lms 每次只计算一个小block的C的值, block大小为mcxn，然后通过innerKernel计算固定的固定的4x4
 * @param m
 * @param n
 * @param k
 * @param a
 * @param lda
 * @param b
 * @param ldb
 * @param c
 * @param ldc
 */
void MY_MMult(int m, int n, int k, double *a, int lda,
              double *b, int ldb,
              double *c, int ldc)
{
    int i, j, p, pb, ib;

    /* This time, we compute a mc x n block of C by a call to the InnerKernel */

    // lms loop method 1
    // 方法1 才是每次计算mc x n块的值
    //    for (int i = 0; i < m; i += mc)
    //    {
    //        ib = min(m - i, mc);
    //        for (p = 0; p < k; p += kc)
    //        {
    //            pb = min(k - p, kc);
    //            InnerKernel(ib, n, pb, &A(i, p), lda, &B(p, 0), ldb, &C(i, 0), ldc);
    //        }
    //    }

    // lms loop method 2
    // 方法2相当于每次计算 mc x n 的kc/k 份，最终通过外层循环加起来
    for (p = 0; p < k; p += kc)
    {
        pb = min(k - p, kc);
        for (i = 0; i < m; i += mc)
        {
            ib = min(m - i, mc); // lms 不能被mc整除时，m-i < mc，此时处理剩余的部分
            // lms 计算一个mcxk * k*n -> mcxn 的C块值，同时又在K维度上做了kc的拆分，Inner内部
            InnerKernel(ib, n, pb, &A(i, p), lda, &B(p, 0), ldb, &C(i, 0), ldc);
        }
    }
}

void InnerKernel(int m, int n, int k, double *a, int lda,
                 double *b, int ldb,
                 double *c, int ldc)
{
    int i, j;

    for (j = 0; j < n; j += 4)
    { /* Loop over the columns of C, unrolled by 4 */
        for (i = 0; i < m; i += 4)
        { /* Loop over the rows of C */
            /* Update C( i,j ), C( i,j+1 ), C( i,j+2 ), and C( i,j+3 ) in
               one routine (four inner products) */

            AddDot4x4(k, &A(i, 0), lda, &B(0, j), ldb, &C(i, j), ldc);
        }
    }
}

#include <emmintrin.h> // SSE3
#include <mmintrin.h>
#include <pmmintrin.h> // SSE2
#include <xmmintrin.h> // SSE

typedef union
{
    __m128d v;
    double d[2];
} v2df_t;

void AddDot4x4(int k, double *a, int lda, double *b, int ldb, double *c, int ldc)
{
    /* So, this routine computes a 4x4 block of matrix A

             C( 0, 0 ), C( 0, 1 ), C( 0, 2 ), C( 0, 3 ).
             C( 1, 0 ), C( 1, 1 ), C( 1, 2 ), C( 1, 3 ).
             C( 2, 0 ), C( 2, 1 ), C( 2, 2 ), C( 2, 3 ).
             C( 3, 0 ), C( 3, 1 ), C( 3, 2 ), C( 3, 3 ).

       Notice that this routine is called with c = C( i, j ) in the
       previous routine, so these are actually the elements

             C( i  , j ), C( i  , j+1 ), C( i  , j+2 ), C( i  , j+3 )
             C( i+1, j ), C( i+1, j+1 ), C( i+1, j+2 ), C( i+1, j+3 )
             C( i+2, j ), C( i+2, j+1 ), C( i+2, j+2 ), C( i+2, j+3 )
             C( i+3, j ), C( i+3, j+1 ), C( i+3, j+2 ), C( i+3, j+3 )

       in the original matrix C

       And now we use vector registers and instructions */

    int p;
    v2df_t
        c_00_c_10_vreg,
        c_01_c_11_vreg, c_02_c_12_vreg, c_03_c_13_vreg,
        c_20_c_30_vreg, c_21_c_31_vreg, c_22_c_32_vreg, c_23_c_33_vreg,
        a_0p_a_1p_vreg,
        a_2p_a_3p_vreg,
        b_p0_vreg, b_p1_vreg, b_p2_vreg, b_p3_vreg;

    double
        /* Point to the current elements in the four columns of B */
        *b_p0_pntr,
        *b_p1_pntr, *b_p2_pntr, *b_p3_pntr;

    b_p0_pntr = &B(0, 0);
    b_p1_pntr = &B(0, 1);
    b_p2_pntr = &B(0, 2);
    b_p3_pntr = &B(0, 3);

    c_00_c_10_vreg.v = _mm_setzero_pd();
    c_01_c_11_vreg.v = _mm_setzero_pd();
    c_02_c_12_vreg.v = _mm_setzero_pd();
    c_03_c_13_vreg.v = _mm_setzero_pd();
    c_20_c_30_vreg.v = _mm_setzero_pd();
    c_21_c_31_vreg.v = _mm_setzero_pd();
    c_22_c_32_vreg.v = _mm_setzero_pd();
    c_23_c_33_vreg.v = _mm_setzero_pd();

    for (p = 0; p < k; p++) // k 列
    {
        // lms: lda 很大的话。可能造成cache miss，因为跳着o去取A的元素
        // A[0,p], A[0,p+1]中间可能隔得很远
        a_0p_a_1p_vreg.v = _mm_load_pd((double *) &A(0, p)); // lms load A[0,p] A[1,p]
        a_2p_a_3p_vreg.v = _mm_load_pd((double *) &A(2, p)); // lms load A[2,p], A[3,p]

        b_p0_vreg.v = _mm_loaddup_pd((double *) b_p0_pntr++); /* load and duplicate */
        b_p1_vreg.v = _mm_loaddup_pd((double *) b_p1_pntr++); /* load and duplicate */
        b_p2_vreg.v = _mm_loaddup_pd((double *) b_p2_pntr++); /* load and duplicate */
        b_p3_vreg.v = _mm_loaddup_pd((double *) b_p3_pntr++); /* load and duplicate */

        /* First row and second rows */
        c_00_c_10_vreg.v += a_0p_a_1p_vreg.v * b_p0_vreg.v;
        c_01_c_11_vreg.v += a_0p_a_1p_vreg.v * b_p1_vreg.v;
        c_02_c_12_vreg.v += a_0p_a_1p_vreg.v * b_p2_vreg.v;
        c_03_c_13_vreg.v += a_0p_a_1p_vreg.v * b_p3_vreg.v;

        /* Third and fourth rows */
        c_20_c_30_vreg.v += a_2p_a_3p_vreg.v * b_p0_vreg.v;
        c_21_c_31_vreg.v += a_2p_a_3p_vreg.v * b_p1_vreg.v;
        c_22_c_32_vreg.v += a_2p_a_3p_vreg.v * b_p2_vreg.v;
        c_23_c_33_vreg.v += a_2p_a_3p_vreg.v * b_p3_vreg.v;
    }

    C(0, 0) += c_00_c_10_vreg.d[0];
    C(0, 1) += c_01_c_11_vreg.d[0];
    C(0, 2) += c_02_c_12_vreg.d[0];
    C(0, 3) += c_03_c_13_vreg.d[0];

    C(1, 0) += c_00_c_10_vreg.d[1];
    C(1, 1) += c_01_c_11_vreg.d[1];
    C(1, 2) += c_02_c_12_vreg.d[1];
    C(1, 3) += c_03_c_13_vreg.d[1];

    C(2, 0) += c_20_c_30_vreg.d[0];
    C(2, 1) += c_21_c_31_vreg.d[0];
    C(2, 2) += c_22_c_32_vreg.d[0];
    C(2, 3) += c_23_c_33_vreg.d[0];

    C(3, 0) += c_20_c_30_vreg.d[1];
    C(3, 1) += c_21_c_31_vreg.d[1];
    C(3, 2) += c_22_c_32_vreg.d[1];
    C(3, 3) += c_23_c_33_vreg.d[1];
}
