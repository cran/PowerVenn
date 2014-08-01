#include <R.h>
#include <Rmath.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mvhyper.h"
void C_dmvhyper(int *x, int *nL, int *L, int *n, double *p, int *logp){
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
		*p=C_dhyper(*x,L[0],*n-L[0],L[1],*logp);
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
//			for(l=*x; l <= minL; l++){ //calculate f_m(l), x is fixed
//				//f1[l-*x]=dhyper((double)*x,(double)l,(double)(*n-l),(double)L[*nL-1],(int)0);
//				f1[l-*x]=C_dhyper(*x,l,*n-l,L[*nL-1],i0);
//			}
			l=*x;
			f1[0]=C_dhyper(*x,l,*n-l,L[*nL-1],i0);
			for(l=*x+1; l <= minL; l++){
				f1[l-*x] = f1[l-*x-1] * ((double)(*n-l+1-L[*nL-1]+*x)/(double)(l-*x))  * ((double)l/(double)(*n-l+1));
			}
			continue;
		}
		memcpy ( f0, f1, aSize * sizeof((double) 0) );
		if(*nL-i>=2){
			for(k=*x;k <= minL;k++){ //calculate f_l(k)
				f1[k-*x]=0;
//				for(l=max2(*x,k+L[*nL-i]-*n);l <= k; l++){ //sum over l for each k
//					//f1[k-*x] += dhyper((double)l,(double)L[*nL-i],(double)(*n-L[*nL-i]),(double)k,(int)0) * f0[l-*x];
//					f1[k-*x] += C_dhyper(l,L[*nL-i],*n-L[*nL-i],k,i0) * f0[l-*x];
//				}
				l=max2(*x,k+L[*nL-i]-*n);
				temp = C_dhyper(l,L[*nL-i],*n-L[*nL-i],k,i0);
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
//		for(j=*x;j <= minL;j++){
//			//*p += dhyper((double)j,(double)L[1],(double)(*n-L[1]),(double)L[0],(int)0) * f1[j-*x];
//			*p += C_dhyper(j,L[1],*n-L[1],L[0],i0) * f1[j-*x];
//		}
		j=*x;
		temp=C_dhyper(j,L[1],*n-L[1],L[0],i0);
		*p += temp * f1[j-*x];
		for(j=*x+1;j <= minL;j++){
			temp=temp * ((double)(L[1]-j+1)/(double)j) * ((double)(L[0]-j+1) /(double)(*n-L[1]-L[0]+j));
			*p += temp * f1[j-*x];
		}
	}
	if(*logp>0) *p=log(*p);
	return;
}
double C_dhyper(int x, int w, int b, int n, int logp){
//probability of getting x white balls out of n draws from an urn with w white balls and b black balls
	double result;
	result=C_logChoose(w,x)+C_logChoose(b,n-x);
	result=result-C_logChoose(w+b,n);
	if(logp==0) result=exp(result);
	return(result);
}
