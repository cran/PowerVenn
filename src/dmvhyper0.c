#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include "mvhyper.h"
void C_dmvhyper0(int *x, int *nL, int *L, int *n, double *p, int *logp){
// the differences between C_dmvhyper and C_dmvhyper0 are in the max loop index and dhyper
/*
nL:   number of subsets
L:    subset sizes
n:    background size
x:    number of elements overlap between all subsets
p:    output probability
*/
	int i, j, k, l, iter;
	int aSize=max(L,*nL)-*x+1;
	double f1[aSize], f0[aSize];

	//from inner-most to outer-most
	for(i=1; i <= *nL-1; i++){
		if(i==1){
			for(l=*x; l <= L[*nL-1]; l++){ //dhyper is R function
				f1[l-*x]=dhyper((double)*x,(double)l,(double)(*n-l),(double)L[*nL-1],(int)0);
			}
			continue;
		}
		for(iter=0;iter<aSize;iter++){
			f0[iter]=f1[iter];
		}
		if(*nL-i>2){
			for(k=*x;k<=L[*nL-i-1];k++){ //calculate f_l(k)
				f1[k-*x]=0;
				for(l=max2(*x,k+L[*nL-i]-*n);l <= min2(k,L[*nL-i]); l++){ //sum over l for each k
					f1[k-*x] += dhyper((double)l,(double)L[*nL-i],(double)(*n-L[*nL-i]),(double)k,(int)0) * f0[l-*x];
				}
			}
			continue;
		}
		//second inner-most
		if(*nL-i==2){
			for(j=*x; j <= min2(L[0],L[1]); j++){ //calculate f_k(j)
				f1[j-*x]=0;
				for(k=max2(*x,j+L[*nL-i]-*n);k <= min2(j,L[*nL-i]); k++){ //sum over k for each j
					f1[j-*x] += dhyper((double)k,(double)L[*nL-i],(double)(*n-L[*nL-i]),(double)j,(int)0) * f0[k-*x];
				}
			}
			continue;
		}
		//final integration
		*p=0;
		for(j=*x;j <= min2(L[0],L[1]);j++){
			*p += dhyper((double)j,(double)L[1],(double)(*n-L[1]),(double)L[0],(int)0) * f1[j-*x];
		}
	}
	if(*logp>0) *p=log(*p);
	return;
}
