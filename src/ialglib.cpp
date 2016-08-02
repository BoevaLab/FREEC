/********************************************************************
optimized ALGLIB subroutines.
********************************************************************/

#include "ialglib.h"
//#include "emmintrin.h"
//#include "mmintrin.h"

static const int alglib_simd_alignment = 16;
static const int alglib_r_block        = 32;
static const int alglib_c_block        = 24;
static const int alglib_half_r_block   = alglib_r_block/2;
static const int alglib_half_c_block   = alglib_c_block/2;
static const int alglib_twice_r_block  = alglib_r_block*2;
static const int alglib_twice_c_block  = alglib_c_block*2;
//#define ABLAS_PREFETCH(x) _mm_prefetch((const char*)(x),_MM_HINT_T0)
//#define ABLAS_PREFETCH(x)

static void* alglib_align(void *ptr, size_t alignment)
{
    char *result;
    if( alignment<=1 )
        return ptr;
    result = (char*)ptr;
    if( (result-(char*)0)%alignment!=0 )
        result += alignment - (result-(char*)0)%alignment;
    return result;
}

/********************************************************************
This subroutine calculates fast 32x32 real matrix-vector product:

    y := beta*y + alpha*A*x

using either generic C code or native optimizations (if available)

IMPORTANT:
* A must be stored in row-major order,
  stride is alglib_r_block,
  aligned on alglib_simd_alignment boundary
* X must be aligned on alglib_simd_alignment boundary
* Y may be non-aligned
********************************************************************/
void ialglib::mv_32(const double *a, const double *x, double *y, int stride, double alpha, double beta)
{
    int i, k;
    const double *pa0, *pa1, *pb;

    pa0 = a;
    pa1 = a+alglib_r_block;
    pb = x;
    for(i=0; i<16; i++)
    {
        double v0 = 0, v1 = 0;
        for(k=0; k<4; k++)
        {
            v0 += pa0[0]*pb[0];
            v1 += pa1[0]*pb[0];
            v0 += pa0[1]*pb[1];
            v1 += pa1[1]*pb[1];
            v0 += pa0[2]*pb[2];
            v1 += pa1[2]*pb[2];
            v0 += pa0[3]*pb[3];
            v1 += pa1[3]*pb[3];
            v0 += pa0[4]*pb[4];
            v1 += pa1[4]*pb[4];
            v0 += pa0[5]*pb[5];
            v1 += pa1[5]*pb[5];
            v0 += pa0[6]*pb[6];
            v1 += pa1[6]*pb[6];
            v0 += pa0[7]*pb[7];
            v1 += pa1[7]*pb[7];
            pa0 += 8;
            pa1 += 8;
            pb  += 8;
        }
        y[0] = beta*y[0]+alpha*v0;
        y[stride] = beta*y[stride]+alpha*v1;

        //
        // now we've processed rows I and I+1,
        // pa0 and pa1 are pointing to rows I+1 and I+2.
        // move to I+2 and I+3.
        //
        pa0 += alglib_r_block;
        pa1 += alglib_r_block;
        pb = x;
        y+=2*stride;
    }
}


/********************************************************************
This subroutine calculates fast MxN real matrix-vector product:

    y := beta*y + alpha*A*x

using either generic C code or native optimizations (if available).
It calls mv_32 if both M=32 and N=32.

IMPORTANT:
* 0<=M<=alglib_r_block, 0<=N<=alglib_r_block
* A must be stored in row-major order,
  stride is alglib_r_block

ALIGNMENT REQUIREMENTS:
1. everything may be non-aligned. in such cases we call
   generic C code.
2. for better performance, each row of A should be aligned,
   X should be aligned too. Y may be non-aligned - does not matter.
********************************************************************/
void ialglib::mv(int m, int n, const double *a, const double *x, double *y, int stride, double alpha, double beta)
{
    mv_generic(m, n, a, x, y, stride, alpha, beta);
}

/********************************************************************
This is generic C implementation of mv.
It may work with unaligned data.
********************************************************************/
void ialglib::mv_generic(int m, int n, const double *a, const double *x, double *y, int stride, double alpha, double beta)
{
    if( m==32 && n==32 )
    {
        //
        // 32x32, may be we have something better than general implementation
        //
        mv_32(a, x, y, stride, alpha, beta);
    }
    else
    {
        int i, k, m2, n8, n2, ntrail2;
        const double *pa0, *pa1, *pb;

        //
        // First M/2 rows of A are processed in pairs.
        // Highly optimized code is used.
        //
        m2 = m/2;
        n8 = n/8;
        ntrail2 = (n-8*n8)/2;
        for(i=0; i<m2; i++)
        {
            double v0 = 0, v1 = 0;

            //
            // 'a' points to the part of the matrix which
            // is not processed yet
            //
            pb = x;
            pa0 = a;
            pa1 = a+alglib_r_block;
            a += alglib_twice_r_block;

            //
            // 8 elements per iteration
            //
            for(k=0; k<n8; k++)
            {
                v0 += pa0[0]*pb[0];
                v1 += pa1[0]*pb[0];
                v0 += pa0[1]*pb[1];
                v1 += pa1[1]*pb[1];
                v0 += pa0[2]*pb[2];
                v1 += pa1[2]*pb[2];
                v0 += pa0[3]*pb[3];
                v1 += pa1[3]*pb[3];
                v0 += pa0[4]*pb[4];
                v1 += pa1[4]*pb[4];
                v0 += pa0[5]*pb[5];
                v1 += pa1[5]*pb[5];
                v0 += pa0[6]*pb[6];
                v1 += pa1[6]*pb[6];
                v0 += pa0[7]*pb[7];
                v1 += pa1[7]*pb[7];
                pa0 += 8;
                pa1 += 8;
                pb  += 8;
            }

            //
            // 2 elements per iteration
            //
            for(k=0; k<ntrail2; k++)
            {
                v0 += pa0[0]*pb[0];
                v1 += pa1[0]*pb[0];
                v0 += pa0[1]*pb[1];
                v1 += pa1[1]*pb[1];
                pa0 += 2;
                pa1 += 2;
                pb  += 2;
            }

            //
            // last element, if needed
            //
            if( n%2!=0 )
            {
                v0 += pa0[0]*pb[0];
                v1 += pa1[0]*pb[0];
            }

            //
            // final update
            //
            y[0] = beta*y[0]+alpha*v0;
            y[stride] = beta*y[stride]+alpha*v1;

            //
            // move to the next pair of elements
            //
            y+=2*stride;
        }


        //
        // Last (odd) row is processed with less optimized code.
        //
        if( m%2!=0 )
        {
            double v0 = 0;

            //
            // 'a' points to the part of the matrix which
            // is not processed yet
            //
            pb = x;
            pa0 = a;

            //
            // 2 elements per iteration
            //
            n2 = n/2;
            for(k=0; k<n2; k++)
            {
                v0 += pa0[0]*pb[0]+pa0[1]*pb[1];
                pa0 += 2;
                pb  += 2;
            }

            //
            // last element, if needed
            //
            if( n%2!=0 )
                v0 += pa0[0]*pb[0];

            //
            // final update
            //
            y[0] = beta*y[0]+alpha*v0;
        }
    }
}

/********************************************************************
This subroutine calculates fast MxN complex matrix-vector product:

    y := beta*y + alpha*A*x

using either generic C code or native optimizations (if available).

IMPORTANT:
* 0<=M<=alglib_c_block, 0<=N<=alglib_c_block
* A must be stored in row-major order,
  stride is alglib_c_block,
* Y may be referenced by cy (pointer to ap::complex) or
  dy (pointer to double) depending on what type of output you
  wish. Pass pointer to Y as one of these parameters,
  AND SET OTHER PARAMETER TO NULL.


ALIGNMENT REQUIREMENTS:
1. everything may be non-aligned. in such cases we call
   generic C code.
2. for better performance, each row of A should be aligned,
   X should be aligned too. Y may be non-aligned - does not matter.
********************************************************************/
void ialglib::mv_complex(int m, int n, const double *a, const double *x, ap::complex *cy, double *dy, int stride, ap::complex alpha, ap::complex beta)
{
    mv_complex_generic(m, n, a, x, cy, dy, stride, alpha, beta);
}

/********************************************************************
This is generic C implementation of mv_complex
It may work with unaligned data.
********************************************************************/
void ialglib::mv_complex_generic(int m, int n, const double *a, const double *x, ap::complex *cy, double *dy, int stride, ap::complex alpha, ap::complex beta)
{
    int i, j;
    const double *pa, *parow, *pb;

    parow = a;
    for(i=0; i<m; i++)
    {
        double v0 = 0, v1 = 0;
        pa = parow;
        pb = x;
        for(j=0; j<n; j++)
        {
            v0 += pa[0]*pb[0];
            v1 += pa[0]*pb[1];
            v0 -= pa[1]*pb[1];
            v1 += pa[1]*pb[0];

            pa  += 2;
            pb  += 2;
        }
        if( cy!=NULL )
        {
            double tx = (beta.x*cy->x-beta.y*cy->y)+(alpha.x*v0-alpha.y*v1);
            double ty = (beta.x*cy->y+beta.y*cy->x)+(alpha.x*v1+alpha.y*v0);
            cy->x = tx;
            cy->y = ty;
            cy+=stride;
        }
        else
        {
            double tx = (beta.x*dy[0]-beta.y*dy[1])+(alpha.x*v0-alpha.y*v1);
            double ty = (beta.x*dy[1]+beta.y*dy[0])+(alpha.x*v1+alpha.y*v0);
            dy[0] = tx;
            dy[1] = ty;
            dy += 2*stride;
        }
        parow += 2*alglib_c_block;
    }
}


/********************************************************************
This subroutine sets vector to zero
********************************************************************/
void ialglib::vzero(int n, double *p, int stride)
{
    int i;
    if( stride==1 )
    {
        for(i=0; i<n; i++,p++)
            *p = 0.0;
    }
    else
    {
        for(i=0; i<n; i++,p+=stride)
            *p = 0.0;
    }
}

/********************************************************************
This subroutine sets vector to zero
********************************************************************/
void ialglib::vzero_complex(int n, ap::complex *p, int stride)
{
    int i;
    if( stride==1 )
    {
        for(i=0; i<n; i++,p++)
        {
            p->x = 0.0;
            p->y = 0.0;
        }
    }
    else
    {
        for(i=0; i<n; i++,p+=stride)
        {
            p->x = 0.0;
            p->y = 0.0;
        }
    }
}


/********************************************************************
This subroutine copies unaligned real vector
********************************************************************/
void ialglib::vcopy(int n, const double *a, int stridea, double *b, int strideb)
{
    int i, n2;
    if( stridea==1 && strideb==1 )
    {
        n2 = n/2;
        for(i=n2; i!=0; i--, a+=2, b+=2)
        {
            b[0] = a[0];
            b[1] = a[1];
        }
        if( n%2!=0 )
            b[0] = a[0];
    }
    else
    {
        for(i=0; i<n; i++,a+=stridea,b+=strideb)
            *b = *a;
    }
}


/********************************************************************
This subroutine copies unaligned complex vector
(passed as ap::complex*)

1. strideb is stride measured in complex numbers, not doubles
2. conj may be "N" (no conj.) or "C" (conj.)
********************************************************************/
void ialglib::vcopy_complex(int n, const ap::complex *a, int stridea, double *b, int strideb, char *conj)
{
    int i;

    //
    // more general case
    //
    if( conj[0]=='N' || conj[0]=='n' )
    {
        for(i=0; i<n; i++,a+=stridea,b+=2*strideb)
        {
            b[0] = a->x;
            b[1] = a->y;
        }
    }
    else
    {
        for(i=0; i<n; i++,a+=stridea,b+=2*strideb)
        {
            b[0] = a->x;
            b[1] = -a->y;
        }
    }
}


/********************************************************************
This subroutine copies unaligned complex vector (passed as double*)

1. strideb is stride measured in complex numbers, not doubles
2. conj may be "N" (no conj.) or "C" (conj.)
********************************************************************/
void ialglib::vcopy_complex(int n, const double *a, int stridea, double *b, int strideb, char *conj)
{
    int i;

    //
    // more general case
    //
    if( conj[0]=='N' || conj[0]=='n' )
    {
        for(i=0; i<n; i++,a+=2*stridea,b+=2*strideb)
        {
            b[0] = a[0];
            b[1] = a[1];
        }
    }
    else
    {
        for(i=0; i<n; i++,a+=2*stridea,b+=2*strideb)
        {
            b[0] = a[0];
            b[1] = -a[1];
        }
    }
}


/********************************************************************
This subroutine copies matrix from  non-aligned non-contigous storage
to aligned contigous storage

A:
* MxN
* non-aligned
* non-contigous
* may be transformed during copying (as prescribed by op)

B:
* alglib_r_block*alglib_r_block (only MxN/NxM submatrix is used)
* aligned
* stride is alglib_r_block

Transformation types:
* 0 - no transform
* 1 - transposition
********************************************************************/
void ialglib::mcopyblock(int m, int n, const double *a, int op, int stride, double *b)
{
    int i, j, n2;
    const double *psrc;
    double *pdst;
    if( op==0 )
    {
        n2 = n/2;
        for(i=0,psrc=a; i<m; i++,a+=stride,b+=alglib_r_block,psrc=a)
        {
            for(j=0,pdst=b; j<n2; j++,pdst+=2,psrc+=2)
            {
                pdst[0] = psrc[0];
                pdst[1] = psrc[1];
            }
            if( n%2!=0 )
                pdst[0] = psrc[0];
        }
    }
    else
    {
        n2 = n/2;
        for(i=0,psrc=a; i<m; i++,a+=stride,b+=1,psrc=a)
        {
            for(j=0,pdst=b; j<n2; j++,pdst+=alglib_twice_r_block,psrc+=2)
            {
                pdst[0] = psrc[0];
                pdst[alglib_r_block] = psrc[1];
            }
            if( n%2!=0 )
                pdst[0] = psrc[0];
        }
    }
}


/********************************************************************
This subroutine copies matrix from  aligned contigous storage to non-
aligned non-contigous storage

A:
* MxN
* aligned
* contigous
* stride is alglib_r_block
* may be transformed during copying (as prescribed by op)

B:
* alglib_r_block*alglib_r_block (only MxN/NxM submatrix is used)
* non-aligned, non-contigous

Transformation types:
* 0 - no transform
* 1 - transposition
********************************************************************/
void ialglib::mcopyunblock(int m, int n, const double *a, int op, double *b, int stride)
{
    int i, j, n2;
    const double *psrc;
    double *pdst;
    if( op==0 )
    {
        n2 = n/2;
        for(i=0,psrc=a; i<m; i++,a+=alglib_r_block,b+=stride,psrc=a)
        {
            for(j=0,pdst=b; j<n2; j++,pdst+=2,psrc+=2)
            {
                pdst[0] = psrc[0];
                pdst[1] = psrc[1];
            }
            if( n%2!=0 )
                pdst[0] = psrc[0];
        }
    }
    else
    {
        n2 = n/2;
        for(i=0,psrc=a; i<m; i++,a++,b+=stride,psrc=a)
        {
            for(j=0,pdst=b; j<n2; j++,pdst+=2,psrc+=alglib_twice_r_block)
            {
                pdst[0] = psrc[0];
                pdst[1] = psrc[alglib_r_block];
            }
            if( n%2!=0 )
                pdst[0] = psrc[0];
        }
    }
}


/********************************************************************
This subroutine copies matrix from  non-aligned non-contigous storage
to aligned contigous storage

A:
* MxN
* non-aligned
* non-contigous
* may be transformed during copying (as prescribed by op)
* pointer to ap::complex is passed

B:
* 2*alglib_c_block*alglib_c_block doubles (only MxN/NxM submatrix is used)
* aligned
* stride is alglib_c_block
* pointer to double is passed

Transformation types:
* 0 - no transform
* 1 - transposition
* 2 - conjugate transposition
* 3 - conjugate, but no  transposition
********************************************************************/
void ialglib::mcopyblock_complex(int m, int n, const ap::complex *a, int op, int stride, double *b)
{
    int i, j;
    const ap::complex *psrc;
    double *pdst;
    if( op==0 )
    {
        for(i=0,psrc=a; i<m; i++,a+=stride,b+=alglib_twice_c_block,psrc=a)
            for(j=0,pdst=b; j<n; j++,pdst+=2,psrc++)
            {
                pdst[0] = psrc->x;
                pdst[1] = psrc->y;
            }
    }
    if( op==1 )
    {
        for(i=0,psrc=a; i<m; i++,a+=stride,b+=2,psrc=a)
            for(j=0,pdst=b; j<n; j++,pdst+=alglib_twice_c_block,psrc++)
            {
                pdst[0] = psrc->x;
                pdst[1] = psrc->y;
            }
    }
    if( op==2 )
    {
        for(i=0,psrc=a; i<m; i++,a+=stride,b+=2,psrc=a)
            for(j=0,pdst=b; j<n; j++,pdst+=alglib_twice_c_block,psrc++)
            {
                pdst[0] = psrc->x;
                pdst[1] = -psrc->y;
            }
    }
    if( op==3 )
    {
        for(i=0,psrc=a; i<m; i++,a+=stride,b+=alglib_twice_c_block,psrc=a)
            for(j=0,pdst=b; j<n; j++,pdst+=2,psrc++)
            {
                pdst[0] = psrc->x;
                pdst[1] = -psrc->y;
            }
    }
}


/********************************************************************
This subroutine copies matrix from aligned contigous storage to
non-aligned non-contigous storage

A:
* 2*alglib_c_block*alglib_c_block doubles (only MxN submatrix is used)
* aligned
* stride is alglib_c_block
* pointer to double is passed
* may be transformed during copying (as prescribed by op)

B:
* MxN
* non-aligned
* non-contigous
* pointer to ap::complex is passed

Transformation types:
* 0 - no transform
* 1 - transposition
* 2 - conjugate transposition
* 3 - conjugate, but no  transposition
********************************************************************/
void ialglib::mcopyunblock_complex(int m, int n, const double *a, int op, ap::complex* b, int stride)
{
    int i, j;
    const double *psrc;
    ap::complex *pdst;
    if( op==0 )
    {
        for(i=0,psrc=a; i<m; i++,a+=alglib_twice_c_block,b+=stride,psrc=a)
            for(j=0,pdst=b; j<n; j++,pdst++,psrc+=2)
            {
                pdst->x = psrc[0];
                pdst->y = psrc[1];
            }
    }
    if( op==1 )
    {
        for(i=0,psrc=a; i<m; i++,a+=2,b+=stride,psrc=a)
            for(j=0,pdst=b; j<n; j++,pdst++,psrc+=alglib_twice_c_block)
            {
                pdst->x = psrc[0];
                pdst->y = psrc[1];
            }
    }
    if( op==2 )
    {
        for(i=0,psrc=a; i<m; i++,a+=2,b+=stride,psrc=a)
            for(j=0,pdst=b; j<n; j++,pdst++,psrc+=alglib_twice_c_block)
            {
                pdst->x = psrc[0];
                pdst->y = -psrc[1];
            }
    }
    if( op==3 )
    {
        for(i=0,psrc=a; i<m; i++,a+=alglib_twice_c_block,b+=stride,psrc=a)
            for(j=0,pdst=b; j<n; j++,pdst++,psrc+=2)
            {
                pdst->x = psrc[0];
                pdst->y = -psrc[1];
            }
    }
}


/********************************************************************
This is real GEMM kernel
********************************************************************/
bool ialglib::_i_rmatrixgemmf(int m,
     int n,
     int k,
     double alpha,
     const ap::real_2d_array& _a,
     int ia,
     int ja,
     int optypea,
     const ap::real_2d_array& _b,
     int ib,
     int jb,
     int optypeb,
     double beta,
     ap::real_2d_array& _c,
     int ic,
     int jc)
{
    if( m>alglib_r_block || n>alglib_r_block || k>alglib_r_block )
        return false;

    int i, stride, cstride;
    double *crow;
    double __abuf[alglib_r_block+alglib_simd_alignment];
    double __b[alglib_r_block*alglib_r_block+alglib_simd_alignment];
    double * const abuf = (double * const) alglib_align(__abuf,alglib_simd_alignment);
    double * const b    = (double * const) alglib_align(__b,   alglib_simd_alignment);

    //
    // copy b
    //
    if( optypeb==0 )
        mcopyblock(k, n, &_b(ib,jb), 1, _b.getstride(), b);
    else
        mcopyblock(n, k, &_b(ib,jb), 0, _b.getstride(), b);

    //
    // multiply B by A (from the right, by rows)
    // and store result in C
    //
    crow  = &_c(ic,jc);
    stride = _a.getstride();
    cstride = _c.getstride();
    if( optypea==0 )
    {
        const double *arow = &_a(ia,ja);
        for(i=0; i<m; i++)
        {
            vcopy(k, arow, 1, abuf, 1);
            if( beta==0 )
                vzero(n, crow, 1);
            mv(n, k, b, abuf, crow, 1, alpha, beta);
            crow += cstride;
            arow += stride;
        }
    }
    else
    {
        const double *acol = &_a(ia,ja);
        for(i=0; i<m; i++)
        {
            vcopy(k, acol, stride, abuf, 1);
            if( beta==0 )
                vzero(n, crow, 1);
            mv(n, k, b, abuf, crow, 1, alpha, beta);
            crow += cstride;
            acol++;
        }
    }
    return true;
}


/********************************************************************
complex GEMM kernel
********************************************************************/
bool ialglib::_i_cmatrixgemmf(int m,
     int n,
     int k,
     ap::complex alpha,
     const ap::complex_2d_array& _a,
     int ia,
     int ja,
     int optypea,
     const ap::complex_2d_array& _b,
     int ib,
     int jb,
     int optypeb,
     ap::complex beta,
     ap::complex_2d_array& _c,
     int ic,
     int jc)
 {
    if( m>alglib_c_block || n>alglib_c_block || k>alglib_c_block )
        return false;

    const ap::complex *arow;
    ap::complex *crow;
    int i, stride, cstride;
    double __abuf[2*alglib_c_block+alglib_simd_alignment];
    double __b[2*alglib_c_block*alglib_c_block+alglib_simd_alignment];
    double * const abuf = (double * const) alglib_align(__abuf,alglib_simd_alignment);
    double * const b    = (double * const) alglib_align(__b,   alglib_simd_alignment);

    //
    // copy b
    //
    int brows = optypeb==0 ? k : n;
    int bcols = optypeb==0 ? n : k;
    if( optypeb==0 )
        mcopyblock_complex(brows, bcols, &_b(ib,jb), 1, _b.getstride(), b);
    if( optypeb==1 )
        mcopyblock_complex(brows, bcols, &_b(ib,jb), 0, _b.getstride(), b);
    if( optypeb==2 )
        mcopyblock_complex(brows, bcols, &_b(ib,jb), 3, _b.getstride(), b);

    //
    // multiply B by A (from the right, by rows)
    // and store result in C
    //
    arow  = &_a(ia,ja);
    crow  = &_c(ic,jc);
    stride = _a.getstride();
    cstride = _c.getstride();
    for(i=0; i<m; i++)
    {
        if( optypea==0 )
        {
			vcopy_complex(k, arow, 1, abuf, 1, "No conj");
            arow += stride;
        }
        else if( optypea==1 )
        {
            vcopy_complex(k, arow, stride, abuf, 1, "No conj");
            arow++;
        }
        else
        {
            vcopy_complex(k, arow, stride, abuf, 1, "Conj");
            arow++;
        }
        if( beta==0 )
            vzero_complex(n, crow, 1);
        mv_complex(n, k, b, abuf, crow, NULL, 1, alpha, beta);
        crow += cstride;
    }
    return true;
}

/********************************************************************
complex TRSM kernel
********************************************************************/
bool ialglib::_i_cmatrixrighttrsmf(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2)
{
    if( m>alglib_c_block || n>alglib_c_block )
        return false;


    //
    // local buffers
    //
    double *pdiag;
    int i;
    double __abuf[2*alglib_c_block*alglib_c_block+alglib_simd_alignment];
    double __xbuf[2*alglib_c_block*alglib_c_block+alglib_simd_alignment];
    double __tmpbuf[2*alglib_c_block+alglib_simd_alignment];
    double * const abuf   = (double * const) alglib_align(__abuf,  alglib_simd_alignment);
    double * const xbuf   = (double * const) alglib_align(__xbuf,  alglib_simd_alignment);
    double * const tmpbuf = (double * const) alglib_align(__tmpbuf,alglib_simd_alignment);

    //
    // Prepare
    //
    bool uppera;
    mcopyblock_complex(n, n, &a(i1,j1), optype, a.getstride(), abuf);
    mcopyblock_complex(m, n, &x(i2,j2), 0, x.getstride(), xbuf);
    if( isunit )
        for(i=0,pdiag=abuf; i<n; i++,pdiag+=2*(alglib_c_block+1))
        {
            pdiag[0] = 1.0;
            pdiag[1] = 0.0;
        }
    if( optype==0 )
        uppera = isupper;
    else
        uppera = !isupper;

    //
    // Solve Y*A^-1=X where A is upper or lower triangular
    //
    if( uppera )
    {
        for(i=0,pdiag=abuf; i<n; i++,pdiag+=2*(alglib_c_block+1))
        {
            ap::complex beta = 1.0/ap::complex(pdiag[0],pdiag[1]);
            ap::complex alpha;
            alpha.x = -beta.x;
            alpha.y = -beta.y;
            vcopy_complex(i, abuf+2*i, alglib_c_block, tmpbuf, 1, "No conj");
            mv_complex(m, i, xbuf, tmpbuf, NULL, xbuf+2*i, alglib_c_block, alpha, beta);
        }
        mcopyunblock_complex(m, n, xbuf, 0, &x(i2,j2), x.getstride());
    }
    else
    {
        for(i=n-1,pdiag=abuf+2*((n-1)*alglib_c_block+(n-1)); i>=0; i--,pdiag-=2*(alglib_c_block+1))
        {
            ap::complex beta = 1.0/ap::complex(pdiag[0],pdiag[1]);
            ap::complex alpha;
            alpha.x = -beta.x;
            alpha.y = -beta.y;
            vcopy_complex(n-1-i, pdiag+2*alglib_c_block, alglib_c_block, tmpbuf, 1, "No conj");
            mv_complex(m, n-1-i, xbuf+2*(i+1), tmpbuf, NULL, xbuf+2*i, alglib_c_block, alpha, beta);
        }
        mcopyunblock_complex(m, n, xbuf, 0, &x(i2,j2), x.getstride());
    }
    return true;
}

/********************************************************************
real TRSM kernel
********************************************************************/
bool ialglib::_i_rmatrixrighttrsmf(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2)
{
    if( m>alglib_r_block || n>alglib_r_block )
        return false;

    //
    // local buffers
    //
    double *pdiag;
    int i;
    double __abuf[alglib_r_block*alglib_r_block+alglib_simd_alignment];
    double __xbuf[alglib_r_block*alglib_r_block+alglib_simd_alignment];
    double __tmpbuf[alglib_r_block+alglib_simd_alignment];
    double * const abuf   = (double * const) alglib_align(__abuf,  alglib_simd_alignment);
    double * const xbuf   = (double * const) alglib_align(__xbuf,  alglib_simd_alignment);
    double * const tmpbuf = (double * const) alglib_align(__tmpbuf,alglib_simd_alignment);

    //
    // Prepare
    //
    bool uppera;
    mcopyblock(n, n, &a(i1,j1), optype, a.getstride(), abuf);
    mcopyblock(m, n, &x(i2,j2), 0, x.getstride(), xbuf);
    if( isunit )
        for(i=0,pdiag=abuf; i<n; i++,pdiag+=alglib_r_block+1)
            *pdiag = 1.0;
    if( optype==0 )
        uppera = isupper;
    else
        uppera = !isupper;

    //
    // Solve Y*A^-1=X where A is upper or lower triangular
    //
    if( uppera )
    {
        for(i=0,pdiag=abuf; i<n; i++,pdiag+=alglib_r_block+1)
        {
            double beta  = 1.0/(*pdiag);
            double alpha = -beta;
            vcopy(i, abuf+i, alglib_r_block, tmpbuf, 1);
            mv(m, i, xbuf, tmpbuf, xbuf+i, alglib_r_block, alpha, beta);
        }
        mcopyunblock(m, n, xbuf, 0, &x(i2,j2), x.getstride());
    }
    else
    {
        for(i=n-1,pdiag=abuf+(n-1)*alglib_r_block+(n-1); i>=0; i--,pdiag-=alglib_r_block+1)
        {
            double beta = 1.0/(*pdiag);
            double alpha = -beta;
            vcopy(n-1-i, pdiag+alglib_r_block, alglib_r_block, tmpbuf, 1);
            mv(m, n-1-i, xbuf+i+1, tmpbuf, xbuf+i, alglib_r_block, alpha, beta);
        }
        mcopyunblock(m, n, xbuf, 0, &x(i2,j2), x.getstride());
    }
    return true;
}

/********************************************************************
complex TRSM kernel
********************************************************************/
bool ialglib::_i_cmatrixlefttrsmf(int m,
     int n,
     const ap::complex_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::complex_2d_array& x,
     int i2,
     int j2)
{
    if( m>alglib_c_block || n>alglib_c_block )
        return false;
    
    //
    // local buffers
    //
    double *pdiag, *arow;
    int i;
    double __abuf[2*alglib_c_block*alglib_c_block+alglib_simd_alignment];
    double __xbuf[2*alglib_c_block*alglib_c_block+alglib_simd_alignment];
    double __tmpbuf[2*alglib_c_block+alglib_simd_alignment];
    double * const abuf   = (double * const) alglib_align(__abuf,  alglib_simd_alignment);
    double * const xbuf   = (double * const) alglib_align(__xbuf,  alglib_simd_alignment);
    double * const tmpbuf = (double * const) alglib_align(__tmpbuf,alglib_simd_alignment);

    //
    // Prepare
    // Transpose X (so we may use mv, which calculates A*x, but not x*A)
    //
    bool uppera;
    mcopyblock_complex(m, m, &a(i1,j1), optype, a.getstride(), abuf);
    mcopyblock_complex(m, n, &x(i2,j2), 1, x.getstride(), xbuf);
    if( isunit )
        for(i=0,pdiag=abuf; i<m; i++,pdiag+=2*(alglib_c_block+1))
        {
            pdiag[0] = 1.0;
            pdiag[1] = 0.0;
        }
    if( optype==0 )
        uppera = isupper;
    else
        uppera = !isupper;

    //
    // Solve A^-1*Y^T=X^T where A is upper or lower triangular
    //
    if( uppera )
    {
        for(i=m-1,pdiag=abuf+2*((m-1)*alglib_c_block+(m-1)); i>=0; i--,pdiag-=2*(alglib_c_block+1))
        {
            ap::complex beta = 1.0/ap::complex(pdiag[0],pdiag[1]);
            ap::complex alpha;
            alpha.x = -beta.x;
            alpha.y = -beta.y;
            vcopy_complex(m-1-i, pdiag+2, 1, tmpbuf, 1, "No conj");
            mv_complex(n, m-1-i, xbuf+2*(i+1), tmpbuf, NULL, xbuf+2*i, alglib_c_block, alpha, beta);
        }
        mcopyunblock_complex(m, n, xbuf, 1, &x(i2,j2), x.getstride());
    }
    else
    {   for(i=0,pdiag=abuf,arow=abuf; i<m; i++,pdiag+=2*(alglib_c_block+1),arow+=2*alglib_c_block)
        {
            ap::complex beta = 1.0/ap::complex(pdiag[0],pdiag[1]);
            ap::complex alpha;
            alpha.x = -beta.x;
            alpha.y = -beta.y;
            vcopy_complex(i, arow, 1, tmpbuf, 1, "No conj");
            mv_complex(n, i, xbuf, tmpbuf, NULL, xbuf+2*i, alglib_c_block, alpha, beta);
        }
        mcopyunblock_complex(m, n, xbuf, 1, &x(i2,j2), x.getstride());
    }
    return true;
}


/********************************************************************
real TRSM kernel
********************************************************************/
bool ialglib::_i_rmatrixlefttrsmf(int m,
     int n,
     const ap::real_2d_array& a,
     int i1,
     int j1,
     bool isupper,
     bool isunit,
     int optype,
     ap::real_2d_array& x,
     int i2,
     int j2)
{
    if( m>alglib_r_block || n>alglib_r_block )
        return false;
    
    //
    // local buffers
    //
    double *pdiag, *arow;
    int i;
    double __abuf[alglib_r_block*alglib_r_block+alglib_simd_alignment];
    double __xbuf[alglib_r_block*alglib_r_block+alglib_simd_alignment];
    double __tmpbuf[alglib_r_block+alglib_simd_alignment];
    double * const abuf   = (double * const) alglib_align(__abuf,  alglib_simd_alignment);
    double * const xbuf   = (double * const) alglib_align(__xbuf,  alglib_simd_alignment);
    double * const tmpbuf = (double * const) alglib_align(__tmpbuf,alglib_simd_alignment);

    //
    // Prepare
    // Transpose X (so we may use mv, which calculates A*x, but not x*A)
    //
    bool uppera;
    mcopyblock(m, m, &a(i1,j1), optype, a.getstride(), abuf);
    mcopyblock(m, n, &x(i2,j2), 1, x.getstride(), xbuf);
    if( isunit )
        for(i=0,pdiag=abuf; i<m; i++,pdiag+=alglib_r_block+1)
            *pdiag = 1.0;
    if( optype==0 )
        uppera = isupper;
    else
        uppera = !isupper;

    //
    // Solve A^-1*Y^T=X^T where A is upper or lower triangular
    //
    if( uppera )
    {
        for(i=m-1,pdiag=abuf+(m-1)*alglib_r_block+(m-1); i>=0; i--,pdiag-=alglib_r_block+1)
        {
            double beta = 1.0/(*pdiag);
            double alpha = -beta;
            vcopy(m-1-i, pdiag+1, 1, tmpbuf, 1);
            mv(n, m-1-i, xbuf+i+1, tmpbuf, xbuf+i, alglib_r_block, alpha, beta);
        }
        mcopyunblock(m, n, xbuf, 1, &x(i2,j2), x.getstride());
    }
    else
    {   for(i=0,pdiag=abuf,arow=abuf; i<m; i++,pdiag+=alglib_r_block+1,arow+=alglib_r_block)
        {
            double beta = 1.0/(*pdiag);
            double alpha = -beta;
            vcopy(i, arow, 1, tmpbuf, 1);
            mv(n, i, xbuf, tmpbuf, xbuf+i, alglib_r_block, alpha, beta);
        }
        mcopyunblock(m, n, xbuf, 1, &x(i2,j2), x.getstride());
    }
    return true;
}


/********************************************************************
complex SYRK kernel
********************************************************************/
bool ialglib::_i_cmatrixsyrkf(int n,
     int k,
     double alpha,
     const ap::complex_2d_array& a,
     int ia,
     int ja,
     int optypea,
     double beta,
     ap::complex_2d_array& c,
     int ic,
     int jc,
     bool isupper)
{
    if( n>alglib_c_block || k>alglib_c_block )
        return false;
    if( n==0 )
        return true;
    
    //
    // local buffers
    //
    double *arow, *crow;
    int i;
    double __abuf[2*alglib_c_block*alglib_c_block+alglib_simd_alignment];
    double __cbuf[2*alglib_c_block*alglib_c_block+alglib_simd_alignment];
    double __tmpbuf[2*alglib_c_block+alglib_simd_alignment];
    double * const abuf   = (double * const) alglib_align(__abuf,  alglib_simd_alignment);
    double * const cbuf   = (double * const) alglib_align(__cbuf,  alglib_simd_alignment);
    double * const tmpbuf = (double * const) alglib_align(__tmpbuf,alglib_simd_alignment);

    //
    // copy A and C, task is transformed to "A*A^H"-form.
    // if beta==0, then C is filled by zeros (and not referenced)
    //
    // alpha==0 or k==0 are correctly processed (A is not referenced)
    //
    if( alpha==0 )
        k = 0;
    if( k>0 )
    {
        if( optypea==0 )
            mcopyblock_complex(n, k, &a(ia,ja), 3, a.getstride(), abuf);
        else
            mcopyblock_complex(k, n, &a(ia,ja), 1, a.getstride(), abuf);
    }
    mcopyblock_complex(n, n, &c(ic,jc), 0, c.getstride(), cbuf);
    if( beta==0 )
    {
        for(i=0,crow=cbuf; i<n; i++,crow+=2*alglib_c_block)
            if( isupper )
                vzero(2*(n-i), crow+2*i, 1);
            else
                vzero(2*(i+1), crow, 1);
    }


    //
    // update C
    //
    if( isupper )
    {
        for(i=0,arow=abuf,crow=cbuf; i<n; i++,arow+=2*alglib_c_block,crow+=2*alglib_c_block)
        {
            vcopy_complex(k, arow, 1, tmpbuf, 1, "Conj");
            mv_complex(n-i, k, arow, tmpbuf, NULL, crow+2*i, 1, alpha, beta);
        }
    }
    else
    {
        for(i=0,arow=abuf,crow=cbuf; i<n; i++,arow+=2*alglib_c_block,crow+=2*alglib_c_block)
        {
            vcopy_complex(k, arow, 1, tmpbuf, 1, "Conj");
            mv_complex(i+1, k, abuf, tmpbuf, NULL, crow, 1, alpha, beta);
        }
    }

    //
    // copy back
    //
    mcopyunblock_complex(n, n, cbuf, 0, &c(ic,jc), c.getstride());

    return true;
}


/********************************************************************
real SYRK kernel
********************************************************************/
bool ialglib::_i_rmatrixsyrkf(int n,
     int k,
     double alpha,
     const ap::real_2d_array& a,
     int ia,
     int ja,
     int optypea,
     double beta,
     ap::real_2d_array& c,
     int ic,
     int jc,
     bool isupper)
{
    if( n>alglib_r_block || k>alglib_r_block )
        return false;
    if( n==0 )
        return true;
    
    //
    // local buffers
    //
    double *arow, *crow;
    int i;
    double __abuf[alglib_r_block*alglib_r_block+alglib_simd_alignment];
    double __cbuf[alglib_r_block*alglib_r_block+alglib_simd_alignment];
    double * const abuf   = (double * const) alglib_align(__abuf,  alglib_simd_alignment);
    double * const cbuf   = (double * const) alglib_align(__cbuf,  alglib_simd_alignment);

    //
    // copy A and C, task is transformed to "A*A^T"-form.
    // if beta==0, then C is filled by zeros (and not referenced)
    //
    // alpha==0 or k==0 are correctly processed (A is not referenced)
    //
    if( alpha==0 )
        k = 0;
    if( k>0 )
    {
        if( optypea==0 )
            mcopyblock(n, k, &a(ia,ja), 0, a.getstride(), abuf);
        else
            mcopyblock(k, n, &a(ia,ja), 1, a.getstride(), abuf);
    }
    mcopyblock(n, n, &c(ic,jc), 0, c.getstride(), cbuf);
    if( beta==0 )
    {
        for(i=0,crow=cbuf; i<n; i++,crow+=alglib_r_block)
            if( isupper )
                vzero(n-i, crow+i, 1);
            else
                vzero(i+1, crow, 1);
    }


    //
    // update C
    //
    if( isupper )
    {
        for(i=0,arow=abuf,crow=cbuf; i<n; i++,arow+=alglib_r_block,crow+=alglib_r_block)
        {
            mv(n-i, k, arow, arow, crow+i, 1, alpha, beta);
        }
    }
    else
    {
        for(i=0,arow=abuf,crow=cbuf; i<n; i++,arow+=alglib_r_block,crow+=alglib_r_block)
        {
            mv(i+1, k, abuf, arow, crow, 1, alpha, beta);
        }
    }

    //
    // copy back
    //
    mcopyunblock(n, n, cbuf, 0, &c(ic,jc), c.getstride());

    return true;
}


/********************************************************************
complex rank-1 kernel
********************************************************************/
bool ialglib::_i_cmatrixrank1f(int m,
     int n,
     ap::complex_2d_array& a,
     int ia,
     int ja,
     ap::complex_1d_array& u,
     int uoffs,
     ap::complex_1d_array& v,
     int voffs)
{
    ap::complex *arow, *pu, *pv, *vtmp, *dst;
    int n2 = n/2;
    int stride  = a.getstride();
    int i, j;

    //
    // update pairs of rows
    //
    arow  = &a(ia,ja);
    pu    = &u(uoffs);
    vtmp  = &v(voffs);
    for(i=0; i<m; i++, arow+=stride, pu++)
    {
        //
        // update by two
        //
        for(j=0,pv=vtmp, dst=arow; j<n2; j++, dst+=2, pv+=2)
        {
            double ux  = pu[0].x;
            double uy  = pu[0].y;
            double v0x = pv[0].x;
            double v0y = pv[0].y;
            double v1x = pv[1].x;
            double v1y = pv[1].y;
            dst[0].x += ux*v0x-uy*v0y;
            dst[0].y += ux*v0y+uy*v0x;
            dst[1].x += ux*v1x-uy*v1y;
            dst[1].y += ux*v1y+uy*v1x;
            //dst[0] += pu[0]*pv[0];
            //dst[1] += pu[0]*pv[1];
        }

        //
        // final update
        //
        if( n%2!=0 )
            dst[0] += pu[0]*pv[0];
    }
    return true;
}


/********************************************************************
real rank-1 kernel
********************************************************************/
bool ialglib::_i_rmatrixrank1f(int m,
     int n,
     ap::real_2d_array& a,
     int ia,
     int ja,
     ap::real_1d_array& u,
     int uoffs,
     ap::real_1d_array& v,
     int voffs)
{
    double *arow0, *arow1, *pu, *pv, *vtmp, *dst0, *dst1;
    int m2 = m/2;
    int n2 = n/2;
    int stride  = a.getstride();
    int stride2 = 2*a.getstride();
    int i, j;

    //
    // update pairs of rows
    //
    arow0 = &a(ia,ja);
    arow1 = arow0+stride;
    pu    = &u(uoffs);
    vtmp  = &v(voffs);
    for(i=0; i<m2; i++,arow0+=stride2,arow1+=stride2,pu+=2)
    {
        //
        // update by two
        //
        for(j=0,pv=vtmp, dst0=arow0, dst1=arow1; j<n2; j++, dst0+=2, dst1+=2, pv+=2)
        {
            dst0[0] += pu[0]*pv[0];
            dst0[1] += pu[0]*pv[1];
            dst1[0] += pu[1]*pv[0];
            dst1[1] += pu[1]*pv[1];
        }

        //
        // final update
        //
        if( n%2!=0 )
        {
            dst0[0] += pu[0]*pv[0];
            dst1[0] += pu[1]*pv[0];
        }
    }

    //
    // update last row
    //
    if( m%2!=0 )
    {
        //
        // update by two
        //
        for(j=0,pv=vtmp, dst0=arow0; j<n2; j++, dst0+=2, pv+=2)
        {
            dst0[0] += pu[0]*pv[0];
            dst0[1] += pu[0]*pv[1];
        }

        //
        // final update
        //
        if( n%2!=0 )
            dst0[0] += pu[0]*pv[0];
    }
    return true;
}


/*void ap::mv_32_sse2(const double *a, const double *x, double *y, double alpha, double beta)
{
    int j;
    const double *pa, *pb;

    pa = x;
    pb = a;
    for(j=0; j<16; j++)
    {
        __m128d va0, vb0, vb1, vm0, vm1, vs, vs1, vt;
        double vd, vd1, vd2;

        vs = _mm_setzero_pd();
        vs1 = _mm_setzero_pd();

        va0 = _mm_load_pd(pa+0);
        vb0 = _mm_load_pd(pb+0);
        vb1 = _mm_load_pd(pb+0+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+2);
        vb0 = _mm_load_pd(pb+2);
        vb1 = _mm_load_pd(pb+2+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+4);
        vb0 = _mm_load_pd(pb+4);
        vb1 = _mm_load_pd(pb+4+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+6);
        vb0 = _mm_load_pd(pb+6);
        vb1 = _mm_load_pd(pb+6+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+8);
        vb0 = _mm_load_pd(pb+8);
        vb1 = _mm_load_pd(pb+8+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+10);
        vb0 = _mm_load_pd(pb+10);
        vb1 = _mm_load_pd(pb+10+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+12);
        vb0 = _mm_load_pd(pb+12);
        vb1 = _mm_load_pd(pb+12+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+14);
        vb0 = _mm_load_pd(pb+14);
        vb1 = _mm_load_pd(pb+14+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+16);
        vb0 = _mm_load_pd(pb+16);
        vb1 = _mm_load_pd(pb+16+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+18);
        vb0 = _mm_load_pd(pb+18);
        vb1 = _mm_load_pd(pb+18+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+20);
        vb0 = _mm_load_pd(pb+20);
        vb1 = _mm_load_pd(pb+20+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+22);
        vb0 = _mm_load_pd(pb+22);
        vb1 = _mm_load_pd(pb+22+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+24);
        vb0 = _mm_load_pd(pb+24);
        vb1 = _mm_load_pd(pb+24+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+26);
        vb0 = _mm_load_pd(pb+26);
        vb1 = _mm_load_pd(pb+26+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+28);
        vb0 = _mm_load_pd(pb+28);
        vb1 = _mm_load_pd(pb+28+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);
        va0 = _mm_load_pd(pa+30);
        vb0 = _mm_load_pd(pb+30);
        vb1 = _mm_load_pd(pb+30+32);
        vm0 = _mm_mul_pd(va0,vb0);
        vm1 = _mm_mul_pd(va0,vb1);
        vs  = _mm_add_pd(vs,vm0);
        vs1  = _mm_add_pd(vs1,vm1);

        _mm_storeh_pd(&vd1,vs);
        _mm_storel_pd(&vd2,vs);
        *y = beta*(*y)+alpha*(vd1+vd2);
        _mm_storeh_pd(&vd1,vs1);
        _mm_storel_pd(&vd2,vs1);
        y[1] = beta*y[1]+alpha*(vd1+vd2);

        pb+= 2*32;
        y+=2;
    }
}*/

/*void ap::mv_complex_sse2(int m, int n, const double *a, const double *x, ap::complex *y, ap::complex alpha, ap::complex beta)
{
    int i, j, n2;
    const double *pa, *parow, *pb;
    double v0, v1, v2, v3;
    __m128d mmalpha, mmbeta;
        
    parow = a;
    pb = x;
    n2 = n>>1;
    mmalpha = _mm_load_sd(&alpha.x);
    mmalpha = _mm_loadh_pd(mmalpha, &alpha.y);
    mmbeta  = _mm_load_sd(&beta.x);
    mmbeta  = _mm_loadh_pd(mmbeta, &beta.y);
    for(i=0; i<m; i++)
    {
        __m128d v_ac_bd, v_ad_bc, v00, v01, v10, v11, mma0, mma1, mmb0, mmb1, mmalpha, mmbeta;
        pa = parow;
        v_ac_bd = _mm_setzero_pd();
        v_ad_bc = _mm_setzero_pd();
        for(j=0; j<n2; j++)
        {
            mma0 = _mm_load_pd(pa+0);
            mmb0 = _mm_load_pd(pb+0);
            mma1 = _mm_load_pd(pa+2);
            mmb1 = _mm_load_pd(pb+2);
            v01 = _mm_shuffle_pd(mmb0, mmb0, _MM_SHUFFLE2(0,1));
            v11 = _mm_shuffle_pd(mmb1, mmb1, _MM_SHUFFLE2(0,1));
            v00 = _mm_mul_pd(mma0, mmb0);
            v01 = _mm_mul_pd(mma0, v01);
            v10 = _mm_mul_pd(mma1, mmb1);
            v11 = _mm_mul_pd(mma1, v11);
            v_ac_bd = _mm_add_pd(v_ac_bd , v00);
            v_ad_bc = _mm_add_pd(v_ad_bc , v01);
            v_ac_bd = _mm_add_pd(v_ac_bd , v10);
            v_ad_bc = _mm_add_pd(v_ad_bc , v11);
            pa  += 4;
            pb  += 4;
        }
        if( n%2!=0 )
        {
            mma0 = _mm_load_pd(pa+0);
            mmb0 = _mm_load_pd(pb+0);
            v01 = _mm_shuffle_pd(mmb0, mmb0, _MM_SHUFFLE2(0,1));
            v00 = _mm_mul_pd(mma0, mmb0);
            v01 = _mm_mul_pd(mma0, v01);
            v_ac_bd = _mm_add_pd(v_ac_bd , v00);
            v_ad_bc = _mm_add_pd(v_ad_bc , v01);
        }
        v00 = _mm_shuffle_pd(v_ac_bd, v_ac_bd, _MM_SHUFFLE2(1,1));
        v_ac_bd = _mm_sub_sd(v_ac_bd, v00);
        v00 = _mm_shuffle_pd(v_ad_bc, v_ad_bc, _MM_SHUFFLE2(1,1));
        v_ad_bc = _mm_add_sd(v_ad_bc, v00);
        v00 = _mm_shuffle_pd(v_ac_bd, v_ad_bc, _MM_SHUFFLE2(0,0));
        _mm_storel_pd(&v0, v00);
        _mm_storeh_pd(&v1, v00);
        double tx = (beta.x*y->x-beta.y*y->y)+(alpha.x*v0-alpha.y*v1);
        double ty = (beta.x*y->y+beta.y*y->x)+(alpha.x*v1+alpha.y*v0);
        y->x = tx;
        y->y = ty;
        pb = x;
        y++;
        parow += 2*24;
    }
}*/
