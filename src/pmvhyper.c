#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include "mvhyper.h"
void C_pmvhyper(int *x, int *nL, int *L, int *n, double *p, int *lower, int *logp){
/*
x:     number of elements overlap between all subsets
nL:    number of subsets
L:     subset sizes
n:     background size
p:     output probability
lower: 1, lower tail probability Pr(overlap <= x); 0, upper tail probability Pr(overlap > x)
logp:  return log probability
*/
	int i;
	int i0=0;
	double p0=0.0;
	double logVal[*n];
	int minL=min(L,*nL);

	for(i=1; i<= *n ; i++){
		logVal[i-1]=log((double)i);
	}
	*p=0.0;
	for(i=minL;i >= *x+1;i--){ //iteration from largest to smallest index, more accurate?
		C_dmvhyper_logVal(&i, nL, L, n, &p0, &i0, logVal);
		*p += p0;
	}
	if(*lower>0) *p=1.0-*p;
	if(*logp>0) *p=log(*p);
	return;
}
void C_pmvhyper0(int *x, int *nL, int *L, int *n, double *p, int *lower, int *logp){
// much slower than C_pmvhyper
/*
x:     number of elements overlap between all subsets
nL:    number of subsets
L:     subset sizes
n:     background size
p:     output probability
lower: 1, lower tail probability Pr(overlap <= x); 0, upper tail probability Pr(overlap > x)
logp:  return log probability
*/
	int i;
	int i0=0;
	double p0=0.0;
	int minL=min(L,*nL);

	*p=0.0;
	for(i=minL;i >= *x+1;i--){ //iteration from largest to smallest index, more accurate?
		C_dmvhyper(&i, nL, L, n, &p0, &i0);
		*p += p0;
	}
	if(*lower>0) *p=1.0-*p;
	if(*logp>0) *p=log(*p);
	return;
}
