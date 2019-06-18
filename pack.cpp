#define DB 1
//#define DB 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include "set2d.h"
#ifndef __TCNTRL__
	#define __TCNTRL__
	#include "tcntrl.h"
#endif

#ifndef __GRIDL__
	#define __GRID__
	#include "grid.h"
#endif

using namespace std;

class Grain{
	public:
		double Drng[2];	// diameter
		double Arng[2];	// aspect ratio
		double Trng[2];	// angle 
		double Wd[2];
		std::mt19937_64 mt;	// random number generator
		std::uniform_real_distribution<> Urnd;
		std::normal_distribution<> Grnd;
		void init_rand(int seed);
//		void set_Wd(double w1,double w2);
		void gen();
		double ra,rb,alph,th;
		double x,y;
		void print();
		double area();
	private:
};
void Grain::init_rand(int seed){
	mt=std::mt19937_64(seed);
	Urnd=std::uniform_real_distribution<double>(0.0,1.0); // Uniform distribution(min,max)
	Grnd=std::normal_distribution<double>(0.0,1.0); // Normal distribution(mean,stdev)	
};
//void Grain::set_Wd(double w1,double w2){
//	Wd[0]=w1; Wd[1]=w2;
//}
void Grain::gen(){
	ra=0.5*(Urnd(mt)*(Drng[1]-Drng[0])+Drng[0]);
	alph=Urnd(mt)*(Arng[1]-Arng[0])+Arng[0];
	rb=ra*alph;
	th=Urnd(mt)*(Trng[1]-Trng[0])+Trng[0];
	x=Urnd(mt)*Wd[0];
	y=Urnd(mt)*Wd[1];
};
double Grain::area(){
	double PI=4.0*atan(1.0);
	return(PI*ra*rb);
};
void Grain::print(){
	double PI=4.0*atan(1.0);
	printf("(x,y)=(%lf, %lf), th=%lf[deg]\n",x,y,th/PI*180.);
	printf("(ra,rb)=(%lf, %lf)\n",ra,rb);
	printf("alph=%lf (aspect ratio)\n",alph);
};
#if DB==1
int main(){

	FILE *finp=fopen("pack.inp","r");
	FILE *fl=fopen("pack.erg","w");
	char fn[128]="geom.dat";
	char cbff[128],fout[128];
	double Xa[2],Wd[2];
	double poro;	// design porosity 
	int Lev;	// Quad-tree height
	double T1,T2;	// temperatures (start,end)
	int nstep;	// Number of MC steps
	int seed=-1;
	double PI=4.0*atan(1.0);
	Grain gr;

//	--------------------------------------------------

	fgets(cbff,128,finp);
	fscanf(finp,"%s\n",fout);
	fgets(cbff,128,finp);
	fscanf(finp,"%le,%le\n",Xa,Xa+1);
	fgets(cbff,128,finp);
	fscanf(finp,"%le,%le\n",Wd,Wd+1);
		gr.Wd[0]=Wd[0]; gr.Wd[1]=Wd[1];
	fgets(cbff,128,finp);
	fscanf(finp,"%d\n",&Lev);
	fgets(cbff,128,finp);
	fscanf(finp,"%lf\n",&poro);

	fgets(cbff,128,finp);
	fscanf(finp,"%lf,%lf\n",gr.Drng, gr.Drng+1);
	fgets(cbff,128,finp);
	fscanf(finp,"%lf,%lf\n",gr.Arng, gr.Arng+1);
	fgets(cbff,128,finp);
	fscanf(finp,"%lf,%lf\n",gr.Trng, gr.Trng+1);
		gr.Trng[0]=gr.Trng[0]/180.*PI;
		gr.Trng[1]=gr.Trng[1]/180.*PI;

	fgets(cbff,128,finp);
	fscanf(finp,"%le,%le\n",&T1,&T2);
	fgets(cbff,128,finp);
	fscanf(finp,"%d\n",&nstep);
	printf("nstep=%d\n",nstep);
	fclose(finp);
//	--------------------------------------------------

	Temp_Hist TH(T1,T2,nstep);
	gr.init_rand(seed);
	double S0=0.0,pr=1.0;
	int np=0;
	while(pr>poro){
		gr.gen();
		S0+=gr.area();
		pr=1.0-S0/Wd[0]*Wd[1];
		np++;
	};
	printf("np=%d\n",np);
	printf("minimum achievable porosity =%lf\n",pr);
	gr.init_rand(seed);

	Solid sld(np);
	sld.set_domain(Xa,Wd);
	for(int i=0;i<np;i++){
		gr.gen();
		sld.els[i].set_xc(gr.x,gr.y);
		sld.els[i].phi=gr.th;
		sld.els[i].set_radi(gr.ra,gr.rb);
		sld.els[i].set_bbox();
		sld.isect[i]=false;
	};


	sld.area(Lev);
	double dE_tot=0.0,dE;
	while(TH.cont_iteration){
		TH.inc_Temp_exp();
		dE=sld.MC(TH);
		dE_tot+=dE;
		printf("tau=%lf, dE=%le dE_sum=%le\n",TH.tau(),dE,dE_tot);
		fprintf(fl,"%ld %le %le %le %le\n",TH.istep,TH.tau(),TH.Temp,dE,dE_tot);
		fflush(stdout);
	};
	printf("dE_tot=%lf\n",dE_tot);
	sld.draw(fn,50);
	sld.area(Lev);
	sld.write(fout);

	Tree4 tr4;
	tr4.setup(sld.els,sld.nelp,false,Lev,sld.bbox);
	tr4.draw();
	tr4.clean();
	return(0);
};
#endif

#if DB==0
int main(){
	Temp_Hist TH;
	char fnt[128]="temp_hist.inp";	// Annealing parameter (temperature control)
	TH.load(fnt);

	int np=400;	// number of particles
	int Lev=9;	// Quad-tree height
	double Wd[2]={1.0,1.0}; // Unit Cell Size
	char fn[128]="geom0.dat";
	double dE_tot=0.0,dE;

	Solid sld(np,Wd);
	sld.write(fn);
	puts("initial geometry written\n");
	int tmp;
//	scanf("%d",&tmp);
	//sld.draw(fn,50);

	double Xa[2]={0.0,0.0};
	char fns[128]="solid1.dat";
	/*
	Solid sld;
	sld.load(fns);
	sld.bbox.setup(Xa,Wd);
	*/

	sld.area(Lev);
	FILE *fl=fopen("log.dat","w");
	int i=0;
	while(TH.cont_iteration){
		dE=sld.MC(TH);
		dE_tot+=dE;
		TH.inc_Temp_exp();
		printf("tau=%lf, dE=%le dE=%le\n",TH.tau(),dE,dE_tot);
		fprintf(fl,"%ld %le %le %le %le\n",TH.istep,TH.tau(),TH.Temp,dE,dE_tot);
		fflush(stdout);
	};
	printf("dE_tot=%lf\n",dE_tot);
	sprintf(fn,"geom1.dat");
	sld.draw(fn,50);
	sld.area(Lev);

	sprintf(fn,"solid.dat");
	sld.write(fn);

	Tree4 tr4;
	tr4.setup(sld.els,sld.nelp,false,Lev,sld.bbox);
	tr4.draw();
	tr4.clean();
	return(0);
};
#endif



// Testing Grid and Node classes on which random walks will be taken
/*
int main(){
	int np=100;	// number of particles
	int Lev=9;	// Quad-tree height
	double Wd[2]={1.0,1.0}; // Unit Cell Size
	//Solid sld(np,Wd);
	char fn[128]="solid.dat";
	Solid sld;
	
	double Xa[2]={0.0,0.0};
	sld.load(fn);
	sld.bbox.setup(Xa,Wd);

	Tree4 tr4;
	tr4.setup(sld.els,sld.nelp,false,Lev,sld.bbox);
	tr4.write();
	printf("Tree data written.");
	int tmp;
	scanf("%d",&tmp);

	double xf[2]={0.5,0.384};
	printf("Is interior point? --> %d\n",QtreeFind(&(tr4.qp0),xf));
	tr4.count();
	tr4.set_grid_params();
	printf("Xa=%lf %lf\n",tr4.Xa[0],tr4.Xa[1]);
	printf("Xb=%lf %lf\n",tr4.Xb[0],tr4.Xb[1]);
	printf("Wd=%lf %lf\n",tr4.Wd[0],tr4.Wd[1]);
	printf("dx=%lf %lf\n",tr4.dx[0],tr4.dx[1]);

	Grid gd;
	gd.setup(tr4);


	std::mt19937_64 engine(-2);
	std::uniform_real_distribution<double>MT01(0.0,1.0);
	double xcod,ycod;
	int next,now;
	printf("Start grid=%d\n",now);
	
	int nwk=10000,Nt=20000,inc=100; 
	
	Node **nd0=(Node **)malloc(sizeof(Node*)*nwk);
	int i,j;
	FILE *fp=fopen("rw.out","w");
	for(i=0;i<nwk;i++){
		now=int(MT01(engine)*gd.Ng);
		nd0[i]=&gd.NDs[now];
	};

	fprintf(fp,"%d,%d,%d\n",nwk,Nt,inc);
	fprintf(fp,"%le,%le\n",tr4.dx[0],tr4.dx[1]);
	double x0,y0;
	double tolx=tr4.dx[0]*1.001;
	double toly=tr4.dx[1]*1.001;
	double *ofx=(double *)malloc(sizeof(double)*nwk);
	double *ofy=(double *)malloc(sizeof(double)*nwk);
	printf("tols=%lf %lf\n",tolx,toly);
	for(j=0;j<Nt;j++){
		for(i=0;i<nwk;i++){
			gd.grid_cod(nd0[i]->iad,&x0,&y0);

			next=int(MT01(engine)*nd0[i]->nc);
			nd0[i]=nd0[i]->cnds[next];
			gd.grid_cod(nd0[i]->iad,&xcod,&ycod);

			if(xcod-x0>tolx) ofx[i]-=Wd[0];
			if(x0-xcod>tolx) ofx[i]+=Wd[0];
			if(ycod-y0>toly) ofy[i]-=Wd[1];
			if(y0-ycod>toly) ofy[i]+=Wd[1];
			if(j%inc==0) fprintf(fp,"%lf,%lf\n",xcod+ofx[i],ycod+ofy[i]);
		}
	}
	tr4.draw();
	return(0);
};
*/
