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

#ifndef AP_H
#define AP_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#ifdef __BORLANDC__
#include <list.h>
#include <vector.h>
#else
#include <list>
#include <vector>
#endif

/********************************************************************
Array bounds check
********************************************************************/
#define NO_AP_ASSERT

#ifndef AP_ASSERT     //
#define NO_AP_ASSERT  // This code avoids definition of the
#endif                // both AP_ASSERT and NO_AP_ASSERT symbols
#ifdef NO_AP_ASSERT   //
#ifdef AP_ASSERT      //
#undef NO_AP_ASSERT   //
#endif                //
#endif                //


/********************************************************************
Current environment.
********************************************************************/
#ifndef AP_WIN32
#ifndef AP_UNKNOWN
#define AP_UNKNOWN
#endif
#endif
#ifdef AP_WIN32
#ifdef AP_UNKNOWN
#error Multiple environments are declared!
#endif
#endif

/********************************************************************
This symbol is used for debugging. Do not define it and do not remove
comments.
********************************************************************/
//#define UNSAFE_MEM_COPY


/********************************************************************
Namespace of a standard library AlgoPascal.
********************************************************************/
namespace ap
{

/********************************************************************
Service routines:
    amalloc - allocates an aligned block of size bytes
    afree - frees block allocated by amalloc
    vlen - just alias for n2-n1+1
********************************************************************/
void* amalloc(size_t size, size_t alignment);
void afree(void *block);
int vlen(int n1, int n2);

/********************************************************************
Exception class.
********************************************************************/
class ap_error
{
public:
    ap_error(){};
    ap_error(const char *s){ msg = s; };

    std::string msg;

    static void make_assertion(bool bClause)
        { if(!bClause) throw ap_error(); };
    static void make_assertion(bool bClause, const char *msg)
        { if(!bClause) throw ap_error(msg); };
private:
};

/********************************************************************
Class defining a complex number with double precision.
********************************************************************/
class complex;

class complex
{
public:
    complex():x(0.0),y(0.0){};
    complex(const double &_x):x(_x),y(0.0){};
    complex(const double &_x, const double &_y):x(_x),y(_y){};
    complex(const complex &z):x(z.x),y(z.y){};

    complex& operator= (const double& v){ x  = v; y = 0.0; return *this; };
    complex& operator+=(const double& v){ x += v;          return *this; };
    complex& operator-=(const double& v){ x -= v;          return *this; };
    complex& operator*=(const double& v){ x *= v; y *= v;  return *this; };
    complex& operator/=(const double& v){ x /= v; y /= v;  return *this; };

    complex& operator= (const complex& z){ x  = z.x; y  = z.y; return *this; };
    complex& operator+=(const complex& z){ x += z.x; y += z.y; return *this; };
    complex& operator-=(const complex& z){ x -= z.x; y -= z.y; return *this; };
    complex& operator*=(const complex& z){ double t = x*z.x-y*z.y; y = x*z.y+y*z.x; x = t; return *this; };
    complex& operator/=(const complex& z)
    {
        ap::complex result;
        double e;
        double f;
        if( fabs(z.y)<fabs(z.x) )
        {
            e = z.y/z.x;
            f = z.x+z.y*e;
            result.x = (z.x+z.y*e)/f;
            result.y = (z.y-z.x*e)/f;
        }
        else
        {
            e = z.x/z.y;
            f = z.y+z.x*e;
            result.x = (z.y+z.x*e)/f;
            result.y = (-z.x+z.y*e)/f;
        }
        *this = result;
        return *this;
    };

    double x, y;
};

const complex operator/(const complex& lhs, const complex& rhs);
const bool operator==(const complex& lhs, const complex& rhs);
const bool operator!=(const complex& lhs, const complex& rhs);
const complex operator+(const complex& lhs);
const complex operator-(const complex& lhs);
const complex operator+(const complex& lhs, const complex& rhs);
const complex operator+(const complex& lhs, const double& rhs);
const complex operator+(const double& lhs, const complex& rhs);
const complex operator-(const complex& lhs, const complex& rhs);
const complex operator-(const complex& lhs, const double& rhs);
const complex operator-(const double& lhs, const complex& rhs);
const complex operator*(const complex& lhs, const complex& rhs);
const complex operator*(const complex& lhs, const double& rhs);
const complex operator*(const double& lhs, const complex& rhs);
const complex operator/(const complex& lhs, const complex& rhs);
const complex operator/(const double& lhs, const complex& rhs);
const complex operator/(const complex& lhs, const double& rhs);
const double abscomplex(const complex &z);
const complex conj(const complex &z);
const complex csqr(const complex &z);


/********************************************************************
Templates for vector operations
********************************************************************/
#include "apvt.h"

/********************************************************************
Level 1 BLAS functions
********************************************************************/
double vdotproduct(const double *v0, int stride0, const double *v1, int stride1, int n);
complex vdotproduct(const complex *v0, int stride0, const char *conj0, const complex *v1, int stride1, const char *conj1, int n);

void vmove(double *vdst,  int stride_dst, const double* vsrc,  int stride_src, int n);
void vmove(complex *vdst, int stride_dst, const complex* vsrc, int stride_src, const char *conj_src, int n);

void vmoveneg(double *vdst,  int stride_dst, const double* vsrc,  int stride_src, int n);
void vmoveneg(complex *vdst, int stride_dst, const complex* vsrc, int stride_src, const char *conj_src, int n);

void vmove(double *vdst,  int stride_dst, const double* vsrc,  int stride_src, int n, double alpha);
void vmove(complex *vdst, int stride_dst, const complex* vsrc, int stride_src, const char *conj_src, int n, double alpha);
void vmove(complex *vdst, int stride_dst, const complex* vsrc, int stride_src, const char *conj_src, int n, complex alpha);

void vadd(double *vdst,  int stride_dst, const double *vsrc,  int stride_src, int n);
void vadd(complex *vdst, int stride_dst, const complex *vsrc, int stride_src, const char *conj_src, int n);

void vadd(double *vdst,  int stride_dst, const double *vsrc,  int stride_src, int n, double alpha);
void vadd(complex *vdst, int stride_dst, const complex *vsrc, int stride_src, const char *conj_src, int n, double alpha);
void vadd(complex *vdst, int stride_dst, const complex *vsrc, int stride_src, const char *conj_src, int n, complex alpha);

void vsub(double *vdst,  int stride_dst, const double *vsrc,  int stride_src, int n);
void vsub(complex *vdst, int stride_dst, const complex *vsrc, int stride_src, const char *conj_src, int n);

void vsub(double *vdst,  int stride_dst, const double *vsrc,  int stride_src, int n, double alpha);
void vsub(complex *vdst, int stride_dst, const complex *vsrc, int stride_src, const char *conj_src, int n, double alpha);
void vsub(complex *vdst, int stride_dst, const complex *vsrc, int stride_src, const char *conj_src, int n, complex alpha);

void vmul(double *vdst,  int stride_dst, int n, double alpha);
void vmul(complex *vdst, int stride_dst, int n, double alpha);
void vmul(complex *vdst, int stride_dst, int n, complex alpha);

/********************************************************************
Obsolete BLAS functions
********************************************************************/
double vdotproduct(const double *v1, const double *v2, int N);
complex vdotproduct(const complex *v1, const complex *v2, int N);

void vmove(double *vdst, const double* vsrc, int N);
void vmove(complex *vdst, const complex* vsrc, int N);

void vmoveneg(double *vdst, const double *vsrc, int N);
void vmoveneg(complex *vdst, const complex *vsrc, int N);

void vmove(double *vdst, const double *vsrc, int N, double alpha);
void vmove(complex *vdst, const complex *vsrc, int N, double alpha);
void vmove(complex *vdst, const complex *vsrc, int N, complex alpha);

void vadd(double *vdst, const double *vsrc, int N);
void vadd(complex *vdst, const complex *vsrc, int N);

void vadd(double *vdst, const double *vsrc, int N, double alpha);
void vadd(complex *vdst, const complex *vsrc, int N, double alpha);
void vadd(complex *vdst, const complex *vsrc, int N, complex alpha);

void vsub(double *vdst, const double *vsrc, int N);
void vsub(complex *vdst, const complex *vsrc, int N);

void vsub(double *vdst, const double *vsrc, int N, double alpha);
void vsub(complex *vdst, const complex *vsrc, int N, double alpha);
void vsub(complex *vdst, const complex *vsrc, int N, complex alpha);

void vmul(double *vdst, int N, double alpha);
void vmul(complex *vdst, int N, double alpha);
void vmul(complex *vdst, int N, complex alpha);


/********************************************************************
Template of a dynamical one-dimensional array
********************************************************************/
template<class T, bool Aligned = false>
class template_1d_array
{
public:
    template_1d_array()
    {
        m_Vec=0;
        m_iVecSize = 0;
        m_iLow = 0;
        m_iHigh = -1;
    };

    ~template_1d_array()
    {
        if(m_Vec)
        {
            if( Aligned )
                ap::afree(m_Vec);
            else
                delete[] m_Vec;
        }
    };

    template_1d_array(const template_1d_array &rhs)
    {
        m_Vec=0;
        m_iVecSize = 0;
        m_iLow = 0;
        m_iHigh = -1;
        if( rhs.m_iVecSize!=0 )
            setcontent(rhs.m_iLow, rhs.m_iHigh, rhs.getcontent());
    };


    const template_1d_array& operator=(const template_1d_array &rhs)
    {
        if( this==&rhs )
            return *this;

        if( rhs.m_iVecSize!=0 )
            setcontent(rhs.m_iLow, rhs.m_iHigh, rhs.getcontent());
        else
        {
            m_Vec=0;
            m_iVecSize = 0;
            m_iLow = 0;
            m_iHigh = -1;
        }
        return *this;
    };


    const T& operator()(int i) const
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i>=m_iLow && i<=m_iHigh);
        #endif
        return m_Vec[ i-m_iLow ];
    };


    T& operator()(int i)
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i>=m_iLow && i<=m_iHigh);
        #endif
        return m_Vec[ i-m_iLow ];
    };


    void setbounds( int iLow, int iHigh )
    {
        if(m_Vec)
        {
            if( Aligned )
                ap::afree(m_Vec);
            else
                delete[] m_Vec;
        }
        m_iLow = iLow;
        m_iHigh = iHigh;
        m_iVecSize = iHigh-iLow+1;
        if( Aligned )
            m_Vec = (T*)ap::amalloc((size_t)(m_iVecSize*sizeof(T)), 16);
        else
            m_Vec = new T[(size_t)m_iVecSize];
    };


    void setlength(int iLen)
    {
        setbounds(0, iLen-1);
    }


    void setcontent( int iLow, int iHigh, const T *pContent )
    {
        setbounds(iLow, iHigh);
        for(int i=0; i<m_iVecSize; i++)
            m_Vec[i] = pContent[i];
    };


    T* getcontent()
    {
        return m_Vec;
    };

    const T* getcontent() const
    {
        return m_Vec;
    };


    int getlowbound(int iBoundNum = 0) const
    {
        return m_iLow;
    };


    int gethighbound(int iBoundNum = 0) const
    {
        return m_iHigh;
    };

    raw_vector<T> getvector(int iStart, int iEnd)
    {
        if( iStart>iEnd || wrongIdx(iStart) || wrongIdx(iEnd) )
            return raw_vector<T>(0, 0, 1);
        else
            return raw_vector<T>(m_Vec+iStart-m_iLow, iEnd-iStart+1, 1);
    };


    const_raw_vector<T> getvector(int iStart, int iEnd) const
    {
        if( iStart>iEnd || wrongIdx(iStart) || wrongIdx(iEnd) )
            return const_raw_vector<T>(0, 0, 1);
        else
            return const_raw_vector<T>(m_Vec+iStart-m_iLow, iEnd-iStart+1, 1);
    };
private:
    bool wrongIdx(int i) const { return i<m_iLow || i>m_iHigh; };

    T         *m_Vec;
    long      m_iVecSize;
    long      m_iLow, m_iHigh;
};



/********************************************************************
Template of a dynamical two-dimensional array
********************************************************************/
template<class T, bool Aligned = false>
class template_2d_array
{
public:
    template_2d_array()
    {
        m_Vec=0;
        m_iVecSize=0;
        m_iLow1 = 0;
        m_iHigh1 = -1;
        m_iLow2 = 0;
        m_iHigh2 = -1;
    };

    ~template_2d_array()
    {
        if(m_Vec)
        {
            if( Aligned )
                ap::afree(m_Vec);
            else
                delete[] m_Vec;
        }
    };

    template_2d_array(const template_2d_array &rhs)
    {
        m_Vec=0;
        m_iVecSize=0;
        m_iLow1 = 0;
        m_iHigh1 = -1;
        m_iLow2 = 0;
        m_iHigh2 = -1;
        if( rhs.m_iVecSize!=0 )
        {
            setbounds(rhs.m_iLow1, rhs.m_iHigh1, rhs.m_iLow2, rhs.m_iHigh2);
            for(int i=m_iLow1; i<=m_iHigh1; i++)
                for(int j=m_iLow2; j<=m_iHigh2; j++)
                    operator()(i,j) = rhs(i,j);
                //vmove(&(operator()(i,m_iLow2)), &(rhs(i,m_iLow2)), m_iHigh2-m_iLow2+1);
        }
    };
    const template_2d_array& operator=(const template_2d_array &rhs)
    {
        if( this==&rhs )
            return *this;

        if( rhs.m_iVecSize!=0 )
        {
            setbounds(rhs.m_iLow1, rhs.m_iHigh1, rhs.m_iLow2, rhs.m_iHigh2);
            for(int i=m_iLow1; i<=m_iHigh1; i++)
                for(int j=m_iLow2; j<=m_iHigh2; j++)
                    operator()(i,j) = rhs(i,j);
                //vmove(&(operator()(i,m_iLow2)), &(rhs(i,m_iLow2)), m_iHigh2-m_iLow2+1);
        }
        else
        {
            if(m_Vec)
            {
                if( Aligned )
                    ap::afree(m_Vec);
                else
                    delete[] m_Vec;
            }
            m_Vec=0;
            m_iVecSize=0;
            m_iLow1 = 0;
            m_iHigh1 = -1;
            m_iLow2 = 0;
            m_iHigh2 = -1;
        }
        return *this;
    };

    const T& operator()(int i1, int i2) const
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i1>=m_iLow1 && i1<=m_iHigh1);
        ap_error::make_assertion(i2>=m_iLow2 && i2<=m_iHigh2);
        #endif
        return m_Vec[ m_iConstOffset + i2 +i1*m_iLinearMember];
    };

    T& operator()(int i1, int i2)
    {
        #ifndef NO_AP_ASSERT
        ap_error::make_assertion(i1>=m_iLow1 && i1<=m_iHigh1);
        ap_error::make_assertion(i2>=m_iLow2 && i2<=m_iHigh2);
        #endif
        return m_Vec[ m_iConstOffset + i2 +i1*m_iLinearMember];
    };

    void setbounds( int iLow1, int iHigh1, int iLow2, int iHigh2 )
    {
        if(m_Vec)
        {
            if( Aligned )
                ap::afree(m_Vec);
            else
                delete[] m_Vec;
        }
        int n1 = iHigh1-iLow1+1;
        int n2 = iHigh2-iLow2+1;
        m_iVecSize = n1*n2;
        if( Aligned )
        {
            //if( n2%2!=0 )
            while( (n2*sizeof(T))%16!=0 )
            {
                n2++;
                m_iVecSize += n1;
            }
            m_Vec = (T*)ap::amalloc((size_t)(m_iVecSize*sizeof(T)), 16);
        }
        else
            m_Vec = new T[(size_t)m_iVecSize];
        m_iLow1  = iLow1;
        m_iHigh1 = iHigh1;
        m_iLow2  = iLow2;
        m_iHigh2 = iHigh2;
        m_iConstOffset = -m_iLow2-m_iLow1*n2;
        m_iLinearMember = n2;
    };

    void setlength(int iLen1, int iLen2)
    {
        setbounds(0, iLen1-1, 0, iLen2-1);
    }

    void setcontent( int iLow1, int iHigh1, int iLow2, int iHigh2, const T *pContent )
    {
        setbounds(iLow1, iHigh1, iLow2, iHigh2);
        for(int i=m_iLow1; i<=m_iHigh1; i++, pContent += m_iHigh2-m_iLow2+1)
            for(int j=m_iLow2; j<=m_iHigh2; j++)
                operator()(i,j) = pContent[j-m_iLow2];
            //vmove(&(operator()(i,m_iLow2)), pContent, m_iHigh2-m_iLow2+1);
    };

    int getlowbound(int iBoundNum) const
    {
        return iBoundNum==1 ? m_iLow1 : m_iLow2;
    };

    int gethighbound(int iBoundNum) const
    {
        return iBoundNum==1 ? m_iHigh1 : m_iHigh2;
    };

    raw_vector<T> getcolumn(int iColumn, int iRowStart, int iRowEnd)
    {
        if( (iRowStart>iRowEnd) || wrongColumn(iColumn) || wrongRow(iRowStart) ||wrongRow(iRowEnd) )
            return raw_vector<T>(0, 0, 1);
        else
            return raw_vector<T>(&((*this)(iRowStart, iColumn)), iRowEnd-iRowStart+1, m_iLinearMember);
    };

    raw_vector<T> getrow(int iRow, int iColumnStart, int iColumnEnd)
    {
        if( (iColumnStart>iColumnEnd) || wrongRow(iRow) || wrongColumn(iColumnStart) || wrongColumn(iColumnEnd))
            return raw_vector<T>(0, 0, 1);
        else
            return raw_vector<T>(&((*this)(iRow, iColumnStart)), iColumnEnd-iColumnStart+1, 1);
    };

    const_raw_vector<T> getcolumn(int iColumn, int iRowStart, int iRowEnd) const
    {
        if( (iRowStart>iRowEnd) || wrongColumn(iColumn) || wrongRow(iRowStart) ||wrongRow(iRowEnd) )
            return const_raw_vector<T>(0, 0, 1);
        else
            return const_raw_vector<T>(&((*this)(iRowStart, iColumn)), iRowEnd-iRowStart+1, m_iLinearMember);
    };

    const_raw_vector<T> getrow(int iRow, int iColumnStart, int iColumnEnd) const
    {
        if( (iColumnStart>iColumnEnd) || wrongRow(iRow) || wrongColumn(iColumnStart) || wrongColumn(iColumnEnd))
            return const_raw_vector<T>(0, 0, 1);
        else
            return const_raw_vector<T>(&((*this)(iRow, iColumnStart)), iColumnEnd-iColumnStart+1, 1);
    };

	int getstride() const
    {
        return m_iLinearMember;
    };
private:
    bool wrongRow(int i) const { return i<m_iLow1 || i>m_iHigh1; };
    bool wrongColumn(int j) const { return j<m_iLow2 || j>m_iHigh2; };

    T           *m_Vec;
    long        m_iVecSize;
    long        m_iLow1, m_iLow2, m_iHigh1, m_iHigh2;
    long        m_iConstOffset, m_iLinearMember;
};


typedef template_1d_array<int>          integer_1d_array;
typedef template_1d_array<double,true>  real_1d_array;
typedef template_1d_array<complex>      complex_1d_array;
typedef template_1d_array<bool>         boolean_1d_array;

typedef template_2d_array<int>          integer_2d_array;
typedef template_2d_array<double,true>  real_2d_array;
typedef template_2d_array<complex>      complex_2d_array;
typedef template_2d_array<bool>         boolean_2d_array;


/********************************************************************
dataset information.

can store regression dataset, classification dataset, or non-labeled
task:
* nout==0 means non-labeled task (clustering, for example)
* nout>0 && nclasses==0 means regression task
* nout>0 && nclasses>0 means classification task
********************************************************************/
/*class dataset
{
public:
    dataset():nin(0), nout(0), nclasses(0), trnsize(0), valsize(0), tstsize(0), totalsize(0){};

    int nin, nout, nclasses;

    int trnsize;
    int valsize;
    int tstsize;
    int totalsize;

    ap::real_2d_array trn;
    ap::real_2d_array val;
    ap::real_2d_array tst;
    ap::real_2d_array all;
};

bool opendataset(std::string file, dataset *pdataset);

//
// internal functions
//
std::string strtolower(const std::string &s);
bool readstrings(std::string file, std::list<std::string> *pOutput);
bool readstrings(std::string file, std::list<std::string> *pOutput, std::string comment);
void explodestring(std::string s, char sep, std::vector<std::string> *pOutput);
std::string xtrim(std::string s);*/

/********************************************************************
reverse communication state
********************************************************************/
struct rcommstate
{
    int stage;
    ap::integer_1d_array ia;
    ap::boolean_1d_array ba;
    ap::real_1d_array ra;
    ap::complex_1d_array ca;
};


/********************************************************************
Constants and functions introduced for compatibility with AlgoPascal
********************************************************************/
extern const double machineepsilon;
extern const double maxrealnumber;
extern const double minrealnumber;

int sign(double x);
double randomreal();
int randominteger(int maxv);
int round_f(double x);
int trunc(double x);
int ifloor(double x);
int iceil(double x);
double pi();
double sqr(double x);
int maxint(int m1, int m2);
int minint(int m1, int m2);
double maxreal(double m1, double m2);
double minreal(double m1, double m2);
bool fp_eq(double v1, double v2);
bool fp_neq(double v1, double v2);
bool fp_less(double v1, double v2);
bool fp_less_eq(double v1, double v2);
bool fp_greater(double v1, double v2);
bool fp_greater_eq(double v1, double v2);

}//namespace ap


#endif
