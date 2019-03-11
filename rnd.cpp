#include<stdio.h>
#include<stdlib.h>
#include<random>
#include<algorithm>

int main(){
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> RndR(-1.0,1.0);

	int n=200;
	int i;

	int *chn1=(int *)malloc(sizeof(int)*n);
	int *chn2=(int *)malloc(sizeof(int)*n);
	std::uniform_int_distribution<int> RndI(0,n-1);
	for(i=0;i<n;i++){
		chn1[i]=i;
		chn2[i]=i;
	}

	int tmp,i1,i2;
	for(i=0;i<1000;i++){
		i1=RndI(mt);
		i2=RndI(mt);
		tmp=chn1[i1];
		chn1[i1]=chn1[i2];
		chn1[i2]=tmp;

		i1=RndI(mt);
		i2=RndI(mt);
		tmp=chn2[i1];
		chn2[i1]=chn2[i2];
		chn2[i2]=tmp;

	}
	std::sort(chn1,chn1+n);
	for(i=0;i<n;i++){
		printf("%d %d\n",chn1[i],chn2[i]);
	}
	return(0);
};
