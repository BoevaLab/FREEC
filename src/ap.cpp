/*************************************************************************
AP library 1.3
Copyright (c) 2003-2009 Sergey Bochkanov (ALGLIB project).

>>> LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

#include "ap.h"

const double ap::machineepsilon = 5E-16;
const double ap::maxrealnumber  = 1E300;
const double ap::minrealnumber  = 1E-300;

/********************************************************************
ap::complex operations
********************************************************************/
const bool ap::operator==(const ap::complex& lhs, const ap::complex& rhs)
{
    volatile double x1 = lhs.x;
    volatile double x2 = rhs.x;
    volatile double y1 = lhs.y;
    volatile double y2 = rhs.y;
    return x1==x2 && y1==y2;
}

const bool ap::operator!=(const ap::complex& lhs, const ap::complex& rhs)
{ return !(lhs==rhs); }

const ap::complex ap::operator+(const ap::complex& lhs)
{ return lhs; }

const ap::complex ap::operator-(const ap::complex& lhs)
{ return ap::complex(-lhs.x, -lhs.y); }

const ap::complex ap::operator+(const ap::complex& lhs, const ap::complex& rhs)
{ ap::complex r = lhs; r += rhs; return r; }

const ap::complex ap::operator+(const ap::complex& lhs, const double& rhs)
{ ap::complex r = lhs; r += rhs; return r; }

const ap::complex ap::operator+(const double& lhs, const ap::complex& rhs)
{ ap::complex r = rhs; r += lhs; return r; }

const ap::complex ap::operator-(const ap::complex& lhs, const ap::complex& rhs)
{ ap::complex r = lhs; r -= rhs; return r; }

const ap::complex ap::operator-(const ap::complex& lhs, const double& rhs)
{ ap::complex r = lhs; r -= rhs; return r; }

const ap::complex ap::operator-(const double& lhs, const ap::complex& rhs)
{ ap::complex r = lhs; r -= rhs; return r; }

const ap::complex ap::operator*(const ap::complex& lhs, const ap::complex& rhs)
{ return ap::complex(lhs.x*rhs.x - lhs.y*rhs.y,  lhs.x*rhs.y + lhs.y*rhs.x); }

const ap::complex ap::operator*(const ap::complex& lhs, const double& rhs)
{ return ap::complex(lhs.x*rhs,  lhs.y*rhs); }

const ap::complex ap::operator*(const double& lhs, const ap::complex& rhs)
{ return ap::complex(lhs*rhs.x,  lhs*rhs.y); }

const ap::complex ap::operator/(const ap::complex& lhs, const ap::complex& rhs)
{
    ap::complex result;
    double e;
    double f;
    if( fabs(rhs.y)<fabs(rhs.x) )
    {
        e = rhs.y/rhs.x;
        f = rhs.x+rhs.y*e;
        result.x = (lhs.x+lhs.y*e)/f;
        result.y = (lhs.y-lhs.x*e)/f;
    }
    else
    {
        e = rhs.x/rhs.y;
        f = rhs.y+rhs.x*e;
        result.x = (lhs.y+lhs.x*e)/f;
        result.y = (-lhs.x+lhs.y*e)/f;
    }
    return result;
}

const ap::complex ap::operator/(const double& lhs, const ap::complex& rhs)
{
    ap::complex result;
    double e;
    double f;
    if( fabs(rhs.y)<fabs(rhs.x) )
    {
        e = rhs.y/rhs.x;
        f = rhs.x+rhs.y*e;
        result.x = lhs/f;
        result.y = -lhs*e/f;
    }
    else
    {
        e = rhs.x/rhs.y;
        f = rhs.y+rhs.x*e;
        result.x = lhs*e/f;
        result.y = -lhs/f;
    }
    return result;
}

const ap::complex ap::operator/(const ap::complex& lhs, const double& rhs)
{ return ap::complex(lhs.x/rhs, lhs.y/rhs); }

const double ap::abscomplex(const ap::complex &z)
{
    double w;
    double xabs;
    double yabs;
    double v;

    xabs = fabs(z.x);
    yabs = fabs(z.y);
    w = xabs>yabs ? xabs : yabs;
    v = xabs<yabs ? xabs : yabs;
    if( v==0 )
        return w;
    else
    {
        double t = v/w;
        return w*sqrt(1+t*t);
    }
}

const ap::complex ap::conj(const ap::complex &z)
{ return ap::complex(z.x, -z.y); }

const ap::complex ap::csqr(const ap::complex &z)
{ return ap::complex(z.x*z.x-z.y*z.y, 2*z.x*z.y); }


/********************************************************************
Level 1 BLAS functions
********************************************************************/
double ap::vdotproduct(const double *v0, int stride0, const double *v1, int stride1, int n)
{
    double result = 0;
    int i;
    if( stride0!=1 || stride1!=1 )
    {
        //
        // slow general code
        //
        for(i=0; i<n; i++, v0+=stride0, v1+=stride1)
            result += (*v0)*(*v1);
    }
    else
    {
        //
        // optimized code for stride=1
        //
        int n4 = n/4;
        int nleft = n%4;
        for(i=0; i<n4; i++, v0+=4, v1+=4)
            result += v0[0]*v1[0]+v0[1]*v1[1]+v0[2]*v1[2]+v0[3]*v1[3];
        for(i=0; i<nleft; i++, v0++, v1++)
            result += v0[0]*v1[0];
    }
    return result;
}

ap::complex ap::vdotproduct(const ap::complex *v0, int stride0, const char *conj0, const ap::complex *v1, int stride1, const char *conj1, int n)
{
    double rx = 0, ry = 0;
    int i;
    bool bconj0 = !((conj0[0]=='N') || (conj0[0]=='n'));
    bool bconj1 = !((conj1[0]=='N') || (conj1[0]=='n'));
    if( bconj0 && bconj1 )
    {
        double v0x, v0y, v1x, v1y;
        for(i=0; i<n; i++, v0+=stride0, v1+=stride1)
        {
            v0x = v0->x;
            v0y = -v0->y;
            v1x = v1->x;
            v1y = -v1->y;
            rx += v0x*v1x-v0y*v1y;
            ry += v0x*v1y+v0y*v1x;
        }
    }
    if( !bconj0 && bconj1 )
    {
        double v0x, v0y, v1x, v1y;
        for(i=0; i<n; i++, v0+=stride0, v1+=stride1)
        {
            v0x = v0->x;
            v0y = v0->y;
            v1x = v1->x;
            v1y = -v1->y;
            rx += v0x*v1x-v0y*v1y;
            ry += v0x*v1y+v0y*v1x;
        }
    }
    if( bconj0 && !bconj1 )
    {
        double v0x, v0y, v1x, v1y;
        for(i=0; i<n; i++, v0+=stride0, v1+=stride1)
        {
            v0x = v0->x;
            v0y = -v0->y;
            v1x = v1->x;
            v1y = v1->y;
            rx += v0x*v1x-v0y*v1y;
            ry += v0x*v1y+v0y*v1x;
        }
    }
    if( !bconj0 && !bconj1 )
    {
        double v0x, v0y, v1x, v1y;
        for(i=0; i<n; i++, v0+=stride0, v1+=stride1)
        {
            v0x = v0->x;
            v0y = v0->y;
            v1x = v1->x;
            v1y = v1->y;
            rx += v0x*v1x-v0y*v1y;
            ry += v0x*v1y+v0y*v1x;
        }
    }
    return ap::complex(rx,ry);
}

void ap::vmove(double *vdst, int stride_dst, const double* vsrc,  int stride_src, int n)
{
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            *vdst = *vsrc;
    }
    else
    {
        //
        // highly optimized case
        //
        int n2 = n/2;
        for(i=0; i<n2; i++, vdst+=2, vsrc+=2)
        {
            vdst[0] = vsrc[0];
            vdst[1] = vsrc[1];
        }
        if( n%2!=0 )
            vdst[0] = vsrc[0];
    }
}

void ap::vmove(ap::complex *vdst, int stride_dst, const ap::complex* vsrc, int stride_src, const char *conj_src, int n)
{
    bool bconj = !((conj_src[0]=='N') || (conj_src[0]=='n'));
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x =  vsrc->x;
                vdst->y = -vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
                *vdst = *vsrc;
        }
    }
    else
    {
        //
        // highly optimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x =  vsrc->x;
                vdst->y = -vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
                *vdst = *vsrc;
        }
    }
}

void ap::vmoveneg(double *vdst,  int stride_dst, const double* vsrc,  int stride_src, int n)
{
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            *vdst = -*vsrc;
    }
    else
    {
        //
        // highly optimized case
        //
        int n2 = n/2;
        for(i=0; i<n2; i++, vdst+=2, vsrc+=2)
        {
            vdst[0] = -vsrc[0];
            vdst[1] = -vsrc[1];
        }
        if( n%2!=0 )
            vdst[0] = -vsrc[0];
    }
}

void ap::vmoveneg(ap::complex *vdst, int stride_dst, const ap::complex* vsrc, int stride_src, const char *conj_src, int n)
{
    bool bconj = !((conj_src[0]=='N') || (conj_src[0]=='n'));
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x = -vsrc->x;
                vdst->y =  vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x = -vsrc->x;
                vdst->y = -vsrc->y;
            }
        }
    }
    else
    {
        //
        // highly optimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x = -vsrc->x;
                vdst->y =  vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x = -vsrc->x;
                vdst->y = -vsrc->y;
            }
        }
    }
}

void ap::vmove(double *vdst,  int stride_dst, const double* vsrc,  int stride_src, int n, double alpha)
{
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            *vdst = alpha*(*vsrc);
    }
    else
    {
        //
        // highly optimized case
        //
        int n2 = n/2;
        for(i=0; i<n2; i++, vdst+=2, vsrc+=2)
        {
            vdst[0] = alpha*vsrc[0];
            vdst[1] = alpha*vsrc[1];
        }
        if( n%2!=0 )
            vdst[0] = alpha*vsrc[0];
    }
}

void ap::vmove(ap::complex *vdst, int stride_dst, const ap::complex* vsrc, int stride_src, const char *conj_src, int n, double alpha)
{
    bool bconj = !((conj_src[0]=='N') || (conj_src[0]=='n'));
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x =  alpha*vsrc->x;
                vdst->y = -alpha*vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x = alpha*vsrc->x;
                vdst->y = alpha*vsrc->y;
            }
        }
    }
    else
    {
        //
        // highly optimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x =  alpha*vsrc->x;
                vdst->y = -alpha*vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x = alpha*vsrc->x;
                vdst->y = alpha*vsrc->y;
            }
        }
    }
}

void ap::vmove(ap::complex *vdst, int stride_dst, const ap::complex* vsrc, int stride_src, const char *conj_src, int n, ap::complex alpha)
{
    bool bconj = !((conj_src[0]=='N') || (conj_src[0]=='n'));
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        if( bconj )
        {
            double ax = alpha.x, ay = alpha.y;
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x =  ax*vsrc->x+ay*vsrc->y;
                vdst->y = -ax*vsrc->y+ay*vsrc->x;
            }
        }
        else
        {
            double ax = alpha.x, ay = alpha.y;
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x = ax*vsrc->x-ay*vsrc->y;
                vdst->y = ax*vsrc->y+ay*vsrc->x;
            }
        }
    }
    else
    {
        //
        // highly optimized case
        //
        if( bconj )
        {
            double ax = alpha.x, ay = alpha.y;
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x =  ax*vsrc->x+ay*vsrc->y;
                vdst->y = -ax*vsrc->y+ay*vsrc->x;
            }
        }
        else
        {
            double ax = alpha.x, ay = alpha.y;
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x = ax*vsrc->x-ay*vsrc->y;
                vdst->y = ax*vsrc->y+ay*vsrc->x;
            }
        }
    }
}

void ap::vadd(double *vdst,  int stride_dst, const double *vsrc,  int stride_src, int n)
{
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            *vdst += *vsrc;
    }
    else
    {
        //
        // highly optimized case
        //
        int n2 = n/2;
        for(i=0; i<n2; i++, vdst+=2, vsrc+=2)
        {
            vdst[0] += vsrc[0];
            vdst[1] += vsrc[1];
        }
        if( n%2!=0 )
            vdst[0] += vsrc[0];
    }
}

void ap::vadd(ap::complex *vdst, int stride_dst, const ap::complex *vsrc, int stride_src, const char *conj_src, int n)
{
    bool bconj = !((conj_src[0]=='N') || (conj_src[0]=='n'));
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x += vsrc->x;
                vdst->y -= vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x += vsrc->x;
                vdst->y += vsrc->y;
            }
        }
    }
    else
    {
        //
        // highly optimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x += vsrc->x;
                vdst->y -= vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x += vsrc->x;
                vdst->y += vsrc->y;
            }
        }
    }
}

void ap::vadd(double *vdst,  int stride_dst, const double *vsrc,  int stride_src, int n, double alpha)
{
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            *vdst += alpha*(*vsrc);
    }
    else
    {
        //
        // highly optimized case
        //
        int n2 = n/2;
        for(i=0; i<n2; i++, vdst+=2, vsrc+=2)
        {
            vdst[0] += alpha*vsrc[0];
            vdst[1] += alpha*vsrc[1];
        }
        if( n%2!=0 )
            vdst[0] += alpha*vsrc[0];
    }
}

void ap::vadd(ap::complex *vdst, int stride_dst, const ap::complex *vsrc, int stride_src, const char *conj_src, int n, double alpha)
{
    bool bconj = !((conj_src[0]=='N') || (conj_src[0]=='n'));
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x += alpha*vsrc->x;
                vdst->y -= alpha*vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x += alpha*vsrc->x;
                vdst->y += alpha*vsrc->y;
            }
        }
    }
    else
    {
        //
        // highly optimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x += alpha*vsrc->x;
                vdst->y -= alpha*vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x += alpha*vsrc->x;
                vdst->y += alpha*vsrc->y;
            }
        }
    }
}

void ap::vadd(ap::complex *vdst, int stride_dst, const ap::complex *vsrc, int stride_src, const char *conj_src, int n, ap::complex alpha)
{
    bool bconj = !((conj_src[0]=='N') || (conj_src[0]=='n'));
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        double ax = alpha.x, ay = alpha.y;
        if( bconj )
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x += ax*vsrc->x+ay*vsrc->y;
                vdst->y -= ax*vsrc->y-ay*vsrc->x;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x += ax*vsrc->x-ay*vsrc->y;
                vdst->y += ax*vsrc->y+ay*vsrc->x;
            }
        }
    }
    else
    {
        //
        // highly optimized case
        //
        double ax = alpha.x, ay = alpha.y;
        if( bconj )
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x += ax*vsrc->x+ay*vsrc->y;
                vdst->y -= ax*vsrc->y-ay*vsrc->x;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x += ax*vsrc->x-ay*vsrc->y;
                vdst->y += ax*vsrc->y+ay*vsrc->x;
            }
        }
    }
}

void ap::vsub(double *vdst,  int stride_dst, const double *vsrc,  int stride_src, int n)
{
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            *vdst -= *vsrc;
    }
    else
    {
        //
        // highly optimized case
        //
        int n2 = n/2;
        for(i=0; i<n2; i++, vdst+=2, vsrc+=2)
        {
            vdst[0] -= vsrc[0];
            vdst[1] -= vsrc[1];
        }
        if( n%2!=0 )
            vdst[0] -= vsrc[0];
    }
}

void ap::vsub(ap::complex *vdst, int stride_dst, const ap::complex *vsrc, int stride_src, const char *conj_src, int n)
{
    bool bconj = !((conj_src[0]=='N') || (conj_src[0]=='n'));
    int i;
    if( stride_dst!=1 || stride_src!=1 )
    {
        //
        // general unoptimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x -= vsrc->x;
                vdst->y += vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst+=stride_dst, vsrc+=stride_src)
            {
                vdst->x -= vsrc->x;
                vdst->y -= vsrc->y;
            }
        }
    }
    else
    {
        //
        // highly optimized case
        //
        if( bconj )
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x -= vsrc->x;
                vdst->y += vsrc->y;
            }
        }
        else
        {
            for(i=0; i<n; i++, vdst++, vsrc++)
            {
                vdst->x -= vsrc->x;
                vdst->y -= vsrc->y;
            }
        }
    }
}

void ap::vsub(double *vdst,  int stride_dst, const double *vsrc,  int stride_src, int n, double alpha)
{
    vadd(vdst, stride_dst, vsrc, stride_src, n, -alpha);
}

void ap::vsub(ap::complex *vdst, int stride_dst, const ap::complex *vsrc, int stride_src, const char *conj_src, int n, double alpha)
{
    vadd(vdst, stride_dst, vsrc, stride_src, conj_src, n, -alpha);
}

void ap::vsub(ap::complex *vdst, int stride_dst, const ap::complex *vsrc, int stride_src, const char *conj_src, int n, ap::complex alpha)
{
    vadd(vdst, stride_dst, vsrc, stride_src, conj_src, n, -alpha);
}

void ap::vmul(double *vdst,  int stride_dst, int n, double alpha)
{
    int i;
    if( stride_dst!=1 )
    {
        //
        // general unoptimized case
        //
        for(i=0; i<n; i++, vdst+=stride_dst)
            *vdst *= alpha;
    }
    else
    {
        //
        // highly optimized case
        //
        for(i=0; i<n; i++, vdst++)
            *vdst *= alpha;
    }
}

void ap::vmul(ap::complex *vdst, int stride_dst, int n, double alpha)
{
    int i;
    if( stride_dst!=1 )
    {
        //
        // general unoptimized case
        //
        for(i=0; i<n; i++, vdst+=stride_dst)
        {
            vdst->x *= alpha;
            vdst->y *= alpha;
        }
    }
    else
    {
        //
        // highly optimized case
        //
        for(i=0; i<n; i++, vdst++)
        {
            vdst->x *= alpha;
            vdst->y *= alpha;
        }
    }
}

void ap::vmul(ap::complex *vdst, int stride_dst, int n, ap::complex alpha)
{
    int i;
    if( stride_dst!=1 )
    {
        //
        // general unoptimized case
        //
        double ax = alpha.x, ay = alpha.y;
        for(i=0; i<n; i++, vdst+=stride_dst)
        {
            double  dstx = vdst->x, dsty = vdst->y;
            vdst->x = ax*dstx-ay*dsty;
            vdst->y = ax*dsty+ay*dstx;
        }
    }
    else
    {
        //
        // highly optimized case
        //
        double ax = alpha.x, ay = alpha.y;
        for(i=0; i<n; i++, vdst++)
        {
            double  dstx = vdst->x, dsty = vdst->y;
            vdst->x = ax*dstx-ay*dsty;
            vdst->y = ax*dsty+ay*dstx;
        }
    }
}

/********************************************************************
Obsolete BLAS functions
********************************************************************/
double ap::vdotproduct(const double *v1, const double *v2, int N)
{
    return ap::_vdotproduct<double>(v1, v2, N);
}

ap::complex ap::vdotproduct(const ap::complex *v1, const ap::complex *v2, int N)
{
    return ap::_vdotproduct<ap::complex>(v1, v2, N);
}

void ap::vmove(double *vdst, const double* vsrc, int N)
{
    ap::_vmove<double>(vdst, vsrc, N);
}

void ap::vmove(ap::complex *vdst, const ap::complex* vsrc, int N)
{
    ap::_vmove<ap::complex>(vdst, vsrc, N);
}

void ap::vmoveneg(double *vdst, const double *vsrc, int N)
{
    ap::_vmoveneg<double>(vdst, vsrc, N);
}

void ap::vmoveneg(ap::complex *vdst, const ap::complex *vsrc, int N)
{
    ap::_vmoveneg<ap::complex>(vdst, vsrc, N);
}

void ap::vmove(double *vdst, const double *vsrc, int N, double alpha)
{
    ap::_vmove<double,double>(vdst, vsrc, N, alpha);
}

void ap::vmove(ap::complex *vdst, const ap::complex *vsrc, int N, double alpha)
{
    ap::_vmove<ap::complex,double>(vdst, vsrc, N, alpha);
}

void ap::vmove(ap::complex *vdst, const ap::complex *vsrc, int N, ap::complex alpha)
{
    ap::_vmove<ap::complex,ap::complex>(vdst, vsrc, N, alpha);
}

void ap::vadd(double *vdst, const double *vsrc, int N)
{
    ap::_vadd<double>(vdst, vsrc, N);
}

void ap::vadd(ap::complex *vdst, const ap::complex *vsrc, int N)
{
    ap::_vadd<ap::complex>(vdst, vsrc, N);
}

void ap::vadd(double *vdst, const double *vsrc, int N, double alpha)
{
    ap::_vadd<double,double>(vdst, vsrc, N, alpha);
}

void ap::vadd(ap::complex *vdst, const ap::complex *vsrc, int N, double alpha)
{
    ap::_vadd<ap::complex,double>(vdst, vsrc, N, alpha);
}

void ap::vadd(ap::complex *vdst, const ap::complex *vsrc, int N, ap::complex alpha)
{
    ap::_vadd<ap::complex,ap::complex>(vdst, vsrc, N, alpha);
}

void ap::vsub(double *vdst, const double *vsrc, int N)
{
    ap::_vsub<double>(vdst, vsrc, N);
}

void ap::vsub(ap::complex *vdst, const ap::complex *vsrc, int N)
{
    ap::_vsub<ap::complex>(vdst, vsrc, N);
}

void ap::vsub(double *vdst, const double *vsrc, int N, double alpha)
{
    ap::_vsub<double,double>(vdst, vsrc, N, alpha);
}

void ap::vsub(ap::complex *vdst, const ap::complex *vsrc, int N, double alpha)
{
    ap::_vsub<ap::complex,double>(vdst, vsrc, N, alpha);
}

void ap::vsub(ap::complex *vdst, const ap::complex *vsrc, int N, ap::complex alpha)
{
    ap::_vsub<ap::complex,ap::complex>(vdst, vsrc, N, alpha);
}

void ap::vmul(double *vdst, int N, double alpha)
{
    ap::_vmul<double,double>(vdst, N, alpha);
}

void ap::vmul(ap::complex *vdst, int N, double alpha)
{
    ap::_vmul<ap::complex,double>(vdst, N, alpha);
}

void ap::vmul(ap::complex *vdst, int N, ap::complex alpha)
{
    ap::_vmul<ap::complex,ap::complex>(vdst, N, alpha);
}

/********************************************************************
standard functions
********************************************************************/
int ap::sign(double x)
{
    if( x>0 ) return  1;
    if( x<0 ) return -1;
    return 0;
}

double ap::randomreal()
{
    int i1 = rand();
    int i2 = rand();
    while(i1==RAND_MAX)
        i1 =rand();
    while(i2==RAND_MAX)
        i2 =rand();
    double mx = RAND_MAX;
    return (i1+i2/mx)/mx;
}

int ap::randominteger(int maxv)
{  return rand()%maxv; }

int ap::round_f(double x)
{ return int(floor(x+0.5)); }

int ap::trunc(double x)
{ return int(x>0 ? floor(x) : ceil(x)); }

int ap::ifloor(double x)
{ return int(floor(x)); }

int ap::iceil(double x)
{ return int(ceil(x)); }

double ap::pi()
{ return 3.14159265358979323846; }

double ap::sqr(double x)
{ return x*x; }

int ap::maxint(int m1, int m2)
{
    return m1>m2 ? m1 : m2;
}

int ap::minint(int m1, int m2)
{
    return m1>m2 ? m2 : m1;
}

double ap::maxreal(double m1, double m2)
{
    return m1>m2 ? m1 : m2;
}

double ap::minreal(double m1, double m2)
{
    return m1>m2 ? m2 : m1;
}

bool ap::fp_eq(double v1, double v2)
{
    // IEEE-strict floating point comparison
    volatile double x = v1;
    volatile double y = v2;
    return x==y;
}

bool ap::fp_neq(double v1, double v2)
{
    // IEEE-strict floating point comparison
    return !fp_eq(v1,v2);
}

bool ap::fp_less(double v1, double v2)
{
    // IEEE-strict floating point comparison
    volatile double x = v1;
    volatile double y = v2;
    return x<y;
}

bool ap::fp_less_eq(double v1, double v2)
{
    // IEEE-strict floating point comparison
    volatile double x = v1;
    volatile double y = v2;
    return x<=y;
}

bool ap::fp_greater(double v1, double v2)
{
    // IEEE-strict floating point comparison
    volatile double x = v1;
    volatile double y = v2;
    return x>y;
}

bool ap::fp_greater_eq(double v1, double v2)
{
    // IEEE-strict floating point comparison
    volatile double x = v1;
    volatile double y = v2;
    return x>=y;
}

/********************************************************************
Dataset functions
********************************************************************/
/*bool ap::readstrings(std::string file, std::list<std::string> *pOutput)
{
    return readstrings(file, pOutput, "");
}

bool ap::readstrings(std::string file, std::list<std::string> *pOutput, std::string comment)
{
    std::string cmd, s;
    FILE *f;
    char buf[32768];
    char *str;

    f = fopen(file.c_str(), "rb");
    if( !f )
        return false;
    s = "";
    pOutput->clear();
    while(str=fgets(buf, sizeof(buf), f))
    {
        // TODO: read file by small chunks, combine in one large string
        if( strlen(str)==0 )
            continue;

        //
        // trim trailing newline chars
        //
        char *eos = str+strlen(str)-1;
        if( *eos=='\n' )
        {
            *eos = 0;
            eos--;
        }
        if( *eos=='\r' )
        {
            *eos = 0;
            eos--;
        }
        s = str;

        //
        // skip comments
        //
        if( comment.length()>0 )
            if( strncmp(s.c_str(), comment.c_str(), comment.length())==0 )
            {
                s = "";
                continue;
            }

        //
        // read data
        //
        if( s.length()<1 )
        {
            fclose(f);
            throw ap::ap_error("internal error in read_strings");
        }
        pOutput->push_back(s);
    }
    fclose(f);
    return true;
}

void ap::explodestring(std::string s, char sep, std::vector<std::string> *pOutput)
{
    std::string tmp;
    int i;
    tmp = "";
    pOutput->clear();
    for(i=0; i<s.length(); i++)
    {
        if( s[i]!=sep )
        {
            tmp += s[i];
            continue;
        }
        //if( tmp.length()!=0 )
        pOutput->push_back(tmp);
        tmp = "";
    }
    if( tmp.length()!=0 )
        pOutput->push_back(tmp);
}

std::string ap::strtolower(const std::string &s)
{
    std::string r = s;
    for(int i=0; i<r.length(); i++)
        r[i] = tolower(r[i]);
    return r;
}

std::string ap::xtrim(std::string s)
{
    char *pstr = (char*)malloc(s.length()+1);
    char *p2 = pstr;
    if( pstr==NULL )
        throw "xalloc in xtrim()";
    try
    {
        bool bws;
        int i;

        //
        // special cases:
        // * zero length string
        // * string includes only spaces
        //
        if( s.length()==0 )
        {
            free(pstr);
            return "";
        }
        bws = true;
        for(i=0; i<s.length(); i++)
            if( s[i]!=' ' )
                bws = false;
        if( bws )
        {
            free(pstr);
            return "";
        }

        //
        // merge internal spaces
        //
        bws = false;
        for(i=0; i<s.length(); i++)
        {
            if( s[i]==' ' && bws )
                continue;
            if( s[i]==' ' )
            {
                *p2 = ' ';
                p2++;
                bws = true;
                continue;
            }
            *p2 = s[i];
            bws = false;
            p2++;
        }
        *p2 = 0;

        //
        // trim leading/trailing spaces.
        // we expect at least one non-space character in the string
        //
        p2--;
        while(*p2==' ')
        {
            *p2 = 0;
            p2--;
        }
        p2 = pstr;
        while((*p2)==' ')
            p2++;

        //
        // result
        //
        std::string r = p2;
        free(pstr);
        return r;
    }
    catch(...)
    {
        free(pstr);
        throw "unknown exception in xtrim()";
    }
}

bool ap::opendataset(std::string file, dataset *pdataset)
{
    std::list<std::string> Lines;
    std::vector<std::string> Values, RowsArr, ColsArr, VarsArr, HeadArr;
    std::list<std::string>::iterator i;
    std::string s;
    int TrnFirst, TrnLast, ValFirst, ValLast, TstFirst, TstLast, LinesRead, j;

    //
    // Read data
    //
    if( pdataset==NULL )
        return false;
    if( !readstrings(file, &Lines, "//") )
        return false;
    i = Lines.begin();
    *pdataset = dataset();

    //
    // Read header
    //
    if( i==Lines.end() )
        return false;
    s = ap::xtrim(*i);
    ap::explodestring(s, '#', &HeadArr);
    if( HeadArr.size()!=2 )
        return false;

    //
    // Rows info
    //
    ap::explodestring(ap::xtrim(HeadArr[0]), ' ', &RowsArr);
    if( RowsArr.size()==0 || RowsArr.size()>3 )
        return false;
    if( RowsArr.size()==1 )
    {
        pdataset->totalsize = atol(RowsArr[0].c_str());
        pdataset->trnsize = pdataset->totalsize;
    }
    if( RowsArr.size()==2 )
    {
        pdataset->trnsize = atol(RowsArr[0].c_str());
        pdataset->tstsize = atol(RowsArr[1].c_str());
        pdataset->totalsize = pdataset->trnsize + pdataset->tstsize;
    }
    if( RowsArr.size()==3 )
    {
        pdataset->trnsize = atol(RowsArr[0].c_str());
        pdataset->valsize = atol(RowsArr[1].c_str());
        pdataset->tstsize = atol(RowsArr[2].c_str());
        pdataset->totalsize = pdataset->trnsize + pdataset->valsize + pdataset->tstsize;
    }
    if( pdataset->totalsize<=0 || pdataset->trnsize<0 || pdataset->valsize<0 || pdataset->tstsize<0 )
        return false;
    TrnFirst = 0;
    TrnLast = TrnFirst + pdataset->trnsize;
    ValFirst = TrnLast;
    ValLast = ValFirst + pdataset->valsize;
    TstFirst = ValLast;
    TstLast = TstFirst + pdataset->tstsize;

    //
    // columns
    //
    ap::explodestring(ap::xtrim(HeadArr[1]), ' ', &ColsArr);
    if( ColsArr.size()!=1 && ColsArr.size()!=4 )
        return false;
    if( ColsArr.size()==1 )
    {
        pdataset->nin = atoi(ColsArr[0].c_str());
        if( pdataset->nin<=0 )
            return false;
    }
    if( ColsArr.size()==4 )
    {
        if( ap::strtolower(ColsArr[0])!="reg" && ap::strtolower(ColsArr[0])!="cls" )
            return false;
        if( ColsArr[2]!="=>" )
            return false;
        pdataset->nin = atol(ColsArr[1].c_str());
        if( pdataset->nin<1 )
            return false;
        if( ap::strtolower(ColsArr[0])=="reg" )
        {
            pdataset->nclasses = 0;
            pdataset->nout = atol(ColsArr[3].c_str());
            if( pdataset->nout<1 )
                return false;
        }
        else
        {
            pdataset->nclasses = atol(ColsArr[3].c_str());
            pdataset->nout = 1;
            if( pdataset->nclasses<2 )
                return false;
        }
    }

    //
    // initialize arrays
    //
    pdataset->all.setlength(pdataset->totalsize, pdataset->nin+pdataset->nout);
    if( pdataset->trnsize>0 ) pdataset->trn.setlength(pdataset->trnsize, pdataset->nin+pdataset->nout);
    if( pdataset->valsize>0 ) pdataset->val.setlength(pdataset->valsize, pdataset->nin+pdataset->nout);
    if( pdataset->tstsize>0 ) pdataset->tst.setlength(pdataset->tstsize, pdataset->nin+pdataset->nout);

    //
    // read data
    //
    for(LinesRead=0, i++; i!=Lines.end() && LinesRead<pdataset->totalsize; i++, LinesRead++)
    {
        std::string sss = *i;
        ap::explodestring(ap::xtrim(*i), ' ', &VarsArr);
        if( VarsArr.size()!=pdataset->nin+pdataset->nout )
            return false;
        int tmpc = ap::round(atof(VarsArr[pdataset->nin+pdataset->nout-1].c_str()));
        if( pdataset->nclasses>0 && (tmpc<0 || tmpc>=pdataset->nclasses) )
            return false;
        for(j=0; j<pdataset->nin+pdataset->nout; j++)
        {
            pdataset->all(LinesRead,j) = atof(VarsArr[j].c_str());
            if( LinesRead>=TrnFirst && LinesRead<TrnLast )
                pdataset->trn(LinesRead-TrnFirst,j) = atof(VarsArr[j].c_str());
            if( LinesRead>=ValFirst && LinesRead<ValLast )
                pdataset->val(LinesRead-ValFirst,j) = atof(VarsArr[j].c_str());
            if( LinesRead>=TstFirst && LinesRead<TstLast )
                pdataset->tst(LinesRead-TstFirst,j) = atof(VarsArr[j].c_str());
        }
    }
    if( LinesRead!=pdataset->totalsize )
        return false;
    return true;
}*/

/*
previous variant
bool ap::opendataset(std::string file, dataset *pdataset)
{
    std::list<std::string> Lines;
    std::vector<std::string> Values;
    std::list<std::string>::iterator i;
    int nCol, nRow, nSplitted;
    int nColumns, nRows;

    //
    // Read data
    //
    if( pdataset==NULL )
        return false;
    if( !readstrings(file, &Lines, "//") )
        return false;
    i = Lines.begin();
    *pdataset = dataset();

    //
    // Read columns info
    //
    if( i==Lines.end() )
        return false;
    if( sscanf(i->c_str(), " columns = %d %d ", &pdataset->nin, &pdataset->nout)!=2 )
        return false;
    if( pdataset->nin<=0 || pdataset->nout==0 || pdataset->nout==-1)
        return false;
    if( pdataset->nout<0 )
    {
        pdataset->nclasses = -pdataset->nout;
        pdataset->nout = 1;
        pdataset->iscls = true;
    }
    else
    {
        pdataset->isreg = true;
    }
    nColumns = pdataset->nin+pdataset->nout;
    i++;

    //
    // Read rows info
    //
    if( i==Lines.end() )
        return false;
    if( sscanf(i->c_str(), " rows = %d %d %d ", &pdataset->trnsize, &pdataset->valsize, &pdataset->tstsize)!=3 )
        return false;
    if( (pdataset->trnsize<0) || (pdataset->valsize<0) || (pdataset->tstsize<0) )
        return false;
    if( (pdataset->trnsize==0) && (pdataset->valsize==0) && (pdataset->tstsize==0) )
        return false;
    nRows = pdataset->trnsize+pdataset->valsize+pdataset->tstsize;
    pdataset->size = nRows;
    if( Lines.size()!=nRows+2 )
        return false;
    i++;

    //
    // Read all cases
    //
    ap::real_2d_array &arr = pdataset->all;
    arr.setbounds(0, nRows-1, 0, nColumns-1);
    for(nRow=0; nRow<nRows; nRow++)
    {
        ap::ap_error::make_assertion(i!=Lines.end());
        explodestring(*i, '\t', &Values);
        if( Values.size()!=nColumns )
            return false;
        for(nCol=0; nCol<nColumns; nCol++)
        {
            double v;
            if( sscanf(Values[nCol].c_str(), "%lg", &v)!=1 )
                return false;
            if( (nCol==nColumns-1) && pdataset->iscls && ((round(v)<0) || (round(v)>=pdataset->nclasses)) )
                return false;
            if( (nCol==nColumns-1) && pdataset->iscls )
                arr(nRow, nCol) = round(v);
            else
                arr(nRow, nCol) = v;
        }
        i++;
    }

    //
    // Split to training, validation and test sets
    //
    if( pdataset->trnsize>0 )
        pdataset->trn.setbounds(0, pdataset->trnsize-1, 0, nColumns-1);
    if( pdataset->valsize>0 )
        pdataset->val.setbounds(0, pdataset->valsize-1, 0, nColumns-1);
    if( pdataset->tstsize>0 )
        pdataset->tst.setbounds(0, pdataset->tstsize-1, 0, nColumns-1);
    nSplitted=0;
    for(nRow=0; nRow<=pdataset->trnsize-1; nRow++, nSplitted++)
        for(nCol=0; nCol<=nColumns-1; nCol++)
            pdataset->trn(nRow,nCol) = arr(nSplitted,nCol);
    for(nRow=0; nRow<=pdataset->valsize-1; nRow++, nSplitted++)
        for(nCol=0; nCol<=nColumns-1; nCol++)
            pdataset->val(nRow,nCol) = arr(nSplitted,nCol);
    for(nRow=0; nRow<=pdataset->tstsize-1; nRow++, nSplitted++)
        for(nCol=0; nCol<=nColumns-1; nCol++)
            pdataset->tst(nRow,nCol) = arr(nSplitted,nCol);
    return true;
}*/

/********************************************************************
Service routines:
********************************************************************/
void* ap::amalloc(size_t size, size_t alignment)
{
    if( alignment<=1 )
    {
        //
        // no alignment, just call malloc
        //
        void *block = malloc(sizeof(void*)+size);
        void **p = (void**)block;
        *p = block;
        return (void*)((char*)block+sizeof(void*));
    }
    else
    {
        //
        // align.
        //
        void *block = malloc(alignment-1+sizeof(void*)+size);
        char *result = (char*)block+sizeof(void*);
        //if( ((unsigned int)(result))%alignment!=0 )
        //    result += alignment - ((unsigned int)(result))%alignment;
        if( (result-(char*)0)%alignment!=0 )
            result += alignment - (result-(char*)0)%alignment;
        *((void**)(result-sizeof(void*))) = block;
        return result;
    }
}

void ap::afree(void *block)
{
    void *p = *((void**)((char*)block-sizeof(void*)));
    free(p);
}

int ap::vlen(int n1, int n2)
{
    return n2-n1+1;
}

