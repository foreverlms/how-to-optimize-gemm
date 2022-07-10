/* Create macros so that the matrices are stored in column-major order */

#define A(i, j) a[(j) *lda + (i)]
#define B(i, j) b[(j) *ldb + (i)]
#define C(i, j) c[(j) *ldc + (i)]

/* Routine for computing C = A * B + C */

void AddDot1x4(int, double *, int, double *, int, double *, int);

/**
 * 循环展开 + 减少cache miss
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
    int i, j;

    for (j = 0; j < n; j += 4)
    { /* Loop over the columns of C, unrolled by 4 */
        for (i = 0; i < m; i += 1)
        { /* Loop over the rows of C */
            /* Update C( i,j ), C( i,j+1 ), C( i,j+2 ), and C( i,j+3 ) in
               one routine (four inner products) */

            AddDot1x4(k, &A(i, 0), lda, &B(0, j), ldb, &C(i, j), ldc);
        }
    }
}

void AddDot1x4(int k, double *a, int lda, double *b, int ldb, double *c, int ldc)
{
    /* So, this routine computes four elements of C:

             C( 0, 0 ), C( 0, 1 ), C( 0, 2 ), C( 0, 3 ).

       Notice that this routine is called with c = C( i, j ) in the
       previous routine, so these are actually the elements

             C( i, j ), C( i, j+1 ), C( i, j+2 ), C( i, j+3 )

       in the original matrix C.

       In this version, we merge the four loops, computing four inner
       products simultaneously. */

    int p;

    // 两个提升的地方：
    //  1. The index p needs only be updated once every eight floating point operations.
    //
    // 2. 只在矩阵比较大时会产生明显的效益提升：
    //    A(0,p)只需要一次从内存里取
    //      Element A( 0, p ) needs only be brought in from memory once
    //      instead of four times. (This only becomes a benefit when the
    //      matrices no longer fit in the L2 cache.)
    //

    //  AddDot( k, &A( 0, 0 ), lda, &B( 0, 0 ), &C( 0, 0 ) );
    //  AddDot( k, &A( 0, 0 ), lda, &B( 0, 1 ), &C( 0, 1 ) );
    //  AddDot( k, &A( 0, 0 ), lda, &B( 0, 2 ), &C( 0, 2 ) );
    //  AddDot( k, &A( 0, 0 ), lda, &B( 0, 3 ), &C( 0, 3 ) );
    for (p = 0; p < k; p++)
    {
        C(0, 0) += A(0, p) * B(p, 0);
        C(0, 1) += A(0, p) * B(p, 1);
        C(0, 2) += A(0, p) * B(p, 2);
        C(0, 3) += A(0, p) * B(p, 3);
    }
}
