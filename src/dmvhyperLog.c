#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mvhyper.h"
void C_dmvhyper_logVal(int *x, int *nL, int *L, int *n, double *p, int *logp, double *logVal){
/*
x:     number of elements overlap between all subsets
nL:    number of subsets
L:     subset sizes
n:     background size
p:     output probability
logp:  return log probability
*/
	int i, j, k, l;
	int i0=0;
	int aSize=max(L,*nL)-*x+1;
	double f1[aSize], f0[aSize];
	double temp;
	int minL=min(L,*nL);
	if(*nL==2){
		*p=C_dhyper_logVal(*x,L[0],*n-L[0],L[1],*logp,logVal);
		return;
	}
/*
Recursive update dhyper(x,a,N-a,b):
 dhyper(x,a+1,N-a-1,b)=dhyper(x,a,N-a,b) * ((a+1)/(a+1-x)) * ((N-a-b+x) /(N-a))
 dhyper(x,a,N-a,b)=dhyper(x,a-1,N-a+1,b) * (a/(a-x)) * ((N-a+1-b+x) /(N-a+1))
 dhyper(x+1,a,N-a,b)=dhyper(x,a,N-a,b) * ((a-x)/(x+1)) * ((b-x) /(N-a-b+x+1))
 dhyper(x,a,N-a,b)=dhyper(x-1,a,N-a,b) * ((a-x+1)/x) * ((b-x+1) /(N-a-b+x))
*/
	//from inner-most to outer-most
	for(i=1; i <= *nL-1; i++){
		if(i==1){
			l=*x;
			f1[0]=C_dhyper_logVal(*x,l,*n-l,L[*nL-1],i0,logVal);
			for(l=*x+1; l <= minL; l++){
				f1[l-*x] = f1[l-*x-1] * ((double)(*n-l+1-L[*nL-1]+*x)/(double)(l-*x))  * ((double)l/(double)(*n-l+1));
			}
			continue;
		}
		memcpy ( f0, f1, aSize * sizeof((double) 0) );
		if(*nL-i>=2){
			for(k=*x;k <= minL;k++){ //calculate f_l(k)
				f1[k-*x]=0;
				l=max2(*x,k+L[*nL-i]-*n);
				temp = C_dhyper_logVal(l,L[*nL-i],*n-L[*nL-i],k,i0,logVal);
				f1[k-*x] += temp * f0[l-*x];
				for(l=max2(*x,k+L[*nL-i]-*n)+1;l <= k; l++){ //sum over l for each k
					temp = temp * ((double)(L[*nL-i]-l+1)/(double)l) * ((double)(k-l+1) /(double)(*n-L[*nL-i]-k+l));
					f1[k-*x] += temp * f0[l-*x];
				}
			}
			continue;
		}
		//final integration
		*p=0;
		j=*x;
		temp=C_dhyper_logVal(j,L[1],*n-L[1],L[0],i0,logVal);
		*p += temp * f1[j-*x];
		for(j=*x+1;j <= minL;j++){
			temp=temp * ((double)(L[1]-j+1)/(double)j) * ((double)(L[0]-j+1) /(double)(*n-L[1]-L[0]+j));
			*p += temp * f1[j-*x];
		}
	}
	if(*logp>0) *p=log(*p);
	return;
}
double C_dhyper_logVal(int x, int w, int b, int n, int logp, double *logVal){
//probability of getting x white balls out of n draws from an urn with w white balls and b black balls
	double result;
	result=C_logChoose_logVal(w,x,logVal)+C_logChoose_logVal(b,n-x,logVal);
	result=result-C_logChoose_logVal(w+b,n,logVal);
	if(logp==0) result=exp(result);
	return(result);
}
