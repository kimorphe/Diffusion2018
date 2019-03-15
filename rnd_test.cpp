#define DB 1
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
//https://cpprefjp.github.io/reference/random.html
#include<random>
using namespace std;

double my_rnd(std::mt19937_64 *engine){
	//std::random_device seed_gen;
	//std::mt19937_64 engine(-1);
	std::uniform_real_distribution<double>RndR(-1.0,1.0);
	std::normal_distribution<double>RndG(0.0,0.2);

	return(RndG(*engine));
};

#if DB ==0
int main(){
	std::random_device seed_gen;
	std::mt19937_64 engine(seed_gen());
	std::uniform_real_distribution<double>RndR(-1.0,1.0);
	std::normal_distribution<double>RndG(0.0,0.2);

	int nbin=100;
	int *count=(int *)malloc(sizeof(int)*nbin);
	int i;
	for(i=0;i<nbin;i++) count[i]=0;
	double dx=(1.0-(-1.0))/nbin;
	int ntry=10000;
	for(i=0;i<ntry;i++){
		//count[int((RndR(engine)+1)/dx)]++;
		count[int((RndG(engine)+1)/dx)]++;
	};
	for(i=0;i<nbin;i++) printf("%lf %d\n",-1.0+dx*i,count[i]);

	return(0);
};
#endif

#if DB==1
int main(){
	std::mt19937_64 engine(-1);
	for(int i=0;i<10;i++){
		printf("%lf\n",my_rnd(&engine));
	}
	return(0);
};
#endif
