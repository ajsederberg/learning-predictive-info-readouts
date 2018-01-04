#define GETVALUE(vble)  ( (vble) = getvalue(#vble) )
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "mex.h"
#include "mt19937ar.h"

#define eps 10e-17

int ndata;
int n_yvals;
int n_xvals;
int nfracreps;
int nfracs;

long dum=1;
unsigned long seed;

void calc_info(double *px,double *py,double *pf,double *pixy,double *pixyr,double *ppxy, double *ppxyr)
{
    
    float ran2(long *);
    void shuffle(int *,int);
    long seed2 = -1;
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    int i,j,k,f,ndata_frac,indkf,indij;
    int ipxy,ipxyr,ip,index,ixy_val,ixy_valR,index_i,index_j,index_jR;
    int np = n_xvals*n_yvals;
    double pstep;
    int *isamples;
    int *isamplesR;
    double *pPX;
    double *pPY;
    double *pPYR;
    double *pPXY;
    double *pPXYR;
    double *pPXPY;
    double *pPXPYR;
    
    init_genrand(seed);
    
    pPX    = mxCalloc(n_xvals,sizeof(double));
    pPY    = mxCalloc(n_yvals,sizeof(double));
    pPYR   = mxCalloc(n_yvals,sizeof(double));
    pPXY   = mxCalloc(np,sizeof(double));
    pPXYR  = mxCalloc(np,sizeof(double));
    pPXPY  = mxCalloc(np,sizeof(double));
    pPXPYR = mxCalloc(np,sizeof(double));
    
    isamples = mxCalloc(ndata,sizeof(int));
    isamplesR = mxCalloc(ndata,sizeof(int));
    
    for (i=0;i<ndata;++i){
        isamples[i]=i+1;
        isamplesR[i]=i+1;
    }
    for (k=0;k<nfracreps;++k){
        for (f=0;f<nfracs;++f){
            indkf = k+nfracreps*f;
            pixy[indkf]  = 0.0;
            pixyr[indkf] = 0.0;
        }
    }
    for (i=0;i<np;++i){
        ppxy[i]   = 0.0;
        ppxyr[i]  = 0.0;
        pPXY[i]   = 0.0;
        pPXYR[i]  = 0.0;
        pPXPY[i]  = 0.0;
        pPXPYR[i] = 0.0;
    }
    for (i=0;i<n_xvals;++i){
        pPX[i]  = 0.0;
    }
    for (i=0;i<n_yvals;++i){
        pPY[i]  = 0.0;
        pPYR[i] = 0.0;
    }
    
    for (f=0;f<nfracs;++f){
        ndata_frac = floor(pf[f]*ndata);
        pstep = 1.0/(ndata_frac);
        for (k=0;k<nfracreps;++k){
            shuffle(isamples,ndata);
            shuffle(isamplesR,ndata);
            indkf = k+nfracreps*f;
            if (f==0 && k>0) {
                ip = 0;
                pixy[indkf]  = pixy[ip];
                pixyr[indkf] = pixyr[ip];
            }
            else{
                for (i=0;i<n_xvals;++i){
                    for (j=0;j<n_yvals;++j){
                        indij = i+n_xvals*j;
                        pPXY[indij]  = 0.0;
                        pPXYR[indij] = 0.0;
                        pPXPY[indij] = 0.0;
                    }
                    pPX[i] = 0.0;
                }
                for (i=0;i<n_yvals;++i){
                    pPY[i]  = 0.0;
                    pPYR[i] = 0.0;
                }
                for (index=0;index<ndata_frac;++index){
                    ixy_val  = isamples[index]-1;
                    ixy_valR = isamplesR[index]-1;
                    index_i  = (int)(px[ixy_val])-1;
                    index_j  = (int)(py[ixy_val])-1;
                    index_jR = (int)(py[ixy_valR])-1;
                    ipxy    = index_i+n_xvals*index_j;
                    ipxyr   = index_i+n_xvals*index_jR;
                    pPXY[ipxy]     += pstep;
                    pPXYR[ipxyr]   += pstep;
                    pPX[index_i]   += pstep;
                    pPY[index_j]   += pstep;
                    pPYR[index_jR] += pstep;
                }
                for (i=0;i<n_xvals;++i){
                    for (j=0;j<n_yvals;++j){
                        indij = i+n_xvals*j;
                        pPXPY[indij]  = pPX[i]*pPY[j];
                        pPXPYR[indij] = pPX[i]*pPYR[j];
                    }
                }
                for (i=0;i<np;++i){
                    pixy[indkf]  += pPXY[i]*log((pPXY[i]/(pPXPY[i]+eps))+eps)/log(2);
                    pixyr[indkf] += pPXYR[i]*log((pPXYR[i]/(pPXPYR[i]+eps))+eps)/log(2);
                }
                if (f==0 && k==0){
                    for (i=0;i<np;++i){
                    ppxy[i] = pPXY[i];
                    ppxyr[i] = pPXYR[i];
                    }
                }
                
            }  /* if-else for 1.0 data frac */
        }  /* k frac reps */
    }  /* f fracs */
    
    mxFree(pPX);
    mxFree(pPY);
    mxFree(pPYR);
    mxFree(pPXY);
    mxFree(pPXYR);
    mxFree(pPXPY);
    mxFree(pPXPYR);
    mxFree(isamples);
    mxFree(isamplesR);
    
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    double getvalue();
    void calc_info(double *,double *,double *,double *,double *,double *,double *);
    
    double *pIxyr;
    double *pIxy;
    double *pPxy;
    double *pPxyR;
    double *px;
    double *py;
    double *pfracs;
    
    GETVALUE(seed);
    GETVALUE(ndata);
    GETVALUE(n_yvals);
    GETVALUE(n_xvals);
    GETVALUE(nfracreps);
    GETVALUE(nfracs);
    
    
    px=mxGetPr(prhs[0]);
    py=mxGetPr(prhs[1]);
    pfracs=mxGetPr(prhs[2]);
    
    plhs[0]=mxCreateDoubleMatrix(nfracreps,nfracs,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(nfracreps,nfracs,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(n_xvals,n_yvals,mxREAL);
    plhs[3]=mxCreateDoubleMatrix(n_xvals,n_yvals,mxREAL);
   
    pIxy=mxGetPr(plhs[0]);
    pIxyr=mxGetPr(plhs[1]);
    pPxy=mxGetPr(plhs[2]);
    pPxyR=mxGetPr(plhs[3]);
    
    calc_info(px,py,pfracs,pIxy,pIxyr,pPxy,pPxyR);
    
}

/* function adapted from Anton Krukowski's code */
double getvalue(char *vble)

{
    const mxArray *pscalar;
    double val;
    
    if ( (pscalar = mexGetVariablePtr("caller",vble)) == NULL) {
        mexPrintf("Couldn't find MATLAB matrix %s.\n",vble);
    }
    else {
        val = mxGetScalar(pscalar);
        return(val);
    }
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

void shuffle(int *array, int n)
{
    int i,hold,swap;
    if (n > 1) {
        for (i=0;i<n-1;++i){
            swap = (int)(genrand_real2()*(double)(n-i)+(double)i);
            hold = array[swap];
            array[swap] = array[i];
            array[i] = hold;
        }
    }
}

float ran2(long *idum)
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
#undef MAXINT
/* (C) Copr. 1986-92 Numerical Recipes Software #.3. */


/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.
 
   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).
 
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.
   Copyright (C) 2005, Mutsuo Saito,
   All rights reserved.
 
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
 
     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
 
     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
 
     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.
 
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 
   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
 */

#include <stdio.h>
#include "mt19937ar.h"

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
        (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
        + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
        - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }
    
    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= N) { /* generate N words at one time */
        int kk;
        
        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */
        
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        
        mti = 0;
    }
    
    y = mt[mti++];
    
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    
    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void)
{
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */
