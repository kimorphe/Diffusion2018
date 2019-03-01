#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "domain.h"

//------------------ Pixel Class ------------------
class Pixel{
	public:
		double Xa[2], Xb[2];
		double Wd[2],Xc[2];
		double xs[4],ys[4];
		double Le[4];
		double rmax,rmin;
		void print();
		void draw(char fn[128],char mode[3]);
		void draw();
		void set_Xa(double xa,double ya);
		void set_Xb(double xb,double yb);
		double dist2pt(double px, double py);
		bool ready;
		bool is_in;
		void setup();
		Pixel();
	private:
	protected:
};
Pixel::Pixel(){
	Xa[0]=0.0; Xa[1]=0.0;
	Xb[0]=1.0; Xb[1]=1.0;
	ready=false;
};
void Pixel::draw(char fn[128],char mode[3]){
	FILE *fp=fopen(fn,mode);
	if(~ready) Pixel::setup();
	for(int i=0;i<5;i++){
		fprintf(fp,"%lf %lf\n",xs[i%4],ys[i%4]);
	};
	fprintf(fp,"\n");
	fclose(fp);
};
void Pixel::draw(){
	if(~ready) Pixel::setup();
	for(int i=0;i<5;i++){
		printf("%lf %lf\n",xs[i%4],ys[i%4]);
	};
	printf("\n");
};
void Pixel::print(){
	printf("Xa=%lf %lf\n",Xa[0],Xa[1]);
	printf("Xb=%lf %lf\n",Xb[0],Xb[1]);
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
	Xc[0]=(Xa[0]+Xb[0])*0.5;
	Xc[1]=(Xa[1]+Xb[1])*0.5;
	printf("Wd=%lf %lf\n",Wd[0],Wd[1]);
	printf("Xc=%lf %lf\n",Xc[0],Xc[1]);
};
void Pixel::setup(){
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
	Xc[0]=(Xa[0]+Xb[0])*0.5;
	Xc[1]=(Xa[1]+Xb[1])*0.5;

	xs[0]=Xa[0];
	xs[1]=Xb[0];
	xs[2]=Xb[0];
	xs[3]=Xa[0];

	ys[0]=Xa[1];
	ys[1]=Xa[1];
	ys[2]=Xb[1];
	ys[3]=Xb[1];

	Le[0]=Wd[0];
	Le[1]=Wd[1];
	Le[2]=Wd[0];
	Le[3]=Wd[1];

	ready=true;
};
void Pixel::set_Xa(double xa, double ya){
	Xa[0]=xa;
	Xa[1]=ya;
};
void Pixel::set_Xb(double xb, double yb){
	Xb[0]=xb;
	Xb[1]=yb;
};
double Pixel::dist2pt(double px, double py){
	double x1[2],x2[2];
	int i;
	double r1,r2,r;
	double xf[2],wt[4],hn[4];
	double dx,dy;

	if( ~ready) Pixel::setup();

	wt[0]= px-xs[0];
	wt[2]=-(px-xs[2]);
	wt[1]= py-ys[1];
	wt[3]=-(py-ys[3]);

	hn[0]=-(py-ys[0]);
	hn[2]= (py-ys[2]);
	hn[1]= (px-xs[1]);
	hn[3]=-(px-xs[3]);

	dx=px-xs[0];
	dy=py-ys[0];
	rmin=sqrt(dx*dx+dy*dy);
	rmax=0.0;
	for(i=0;i<4;i++){
		dx=px-xs[i];
		dy=py-ys[i];
		r=sqrt(dx*dx+dy*dy);
		if(r > rmax) rmax=r;
		if(r < rmin) rmin=r;

		if(wt[i]<0.0) continue;
		if(wt[i]>Le[i]) continue;
		if(fabs(hn[i]) < rmin) rmin=fabs(hn[i]);
	};

	is_in=true;
	for(i=0;i<4;i++){
		if(hn[i]>0){
			is_in=false;
			break;
		}
	}
	return(rmin);
};
//---------------------------------------------------
int is_cross(Pixel px, Circ cr){
	double rmin,rmax;

	rmin=px.dist2pt(cr.xc[0],cr.xc[1]);
	rmax=px.rmax;
	if(rmax < cr.radi) return(1);	// P < C
	if(rmin > cr.radi){
		if(px.is_in){
			return(2); // C < P
		}else{
			return(0); // C ^ P=\phi
		}
	}
	return(3); // C^P != \phi
};
//--------------------------------------------------------
class QPatch{
	public:
		Pixel px;	// pixel data
		QPatch *par;	// parent
		QPatch *chld[4]; // children
		QPatch();	// default constructor
		void print();	// 
		void draw(char fn[128],char mode[3]);
		void draw();
		int lev;
		void set_lim(double xa[2], double xb[2]);
	private:
	protected:
};
QPatch::QPatch(){
	par=NULL;
	lev=0;
	for(int i=0;i<4;i++) chld[i]=NULL;
};
void QPatch::print(){
	px.print();
};
void QPatch::draw(char fn[128],char mode[3]){
	px.draw(fn,mode);
};
void QPatch::draw(){
	px.draw();
};
void QPatch::set_lim(double xa[2], double xb[2]){
	px.set_Xa(xa[0],xa[1]);
	px.set_Xb(xb[0],xb[1]);
	px.setup();
};	
QPatch *new_QPatch(double Xa[2], double Wd[2]){
	QPatch *qp=(QPatch *)malloc(sizeof(QPatch));
	double Xb[2];
	Xb[0]=Xa[0]+Wd[0];
	Xb[1]=Xa[1]+Wd[1];
	qp->set_lim(Xa,Xb);
	for(int i=0;i<4;i++) qp->chld[i]=NULL;
	return(qp);
};
//--------------------------------------------------------
void translate_crs(int icrs){

	if(icrs==0) puts("A ^ B = phi");
	if(icrs==1) puts("B < A ");
	if(icrs==2) puts("A < B ");
	if(icrs==3) puts("A ^ B != phi");
};
void gather_leaves(QPatch *qp_par, int *count, QPatch *qp_leaves){
	if(qp_par->chld[0]==NULL){
		//qp_par->draw();
		//printf("idx=%d\n",*count);
		qp_leaves[*count]=(*qp_par);
		(*count)++;
	}else{
		for(int k=0;k<4;k++) gather_leaves(qp_par->chld[k],count,qp_leaves);
	};
};
int Qtree(QPatch *qp, Circ cr,int *count){

	int i,j,k;
	double Xa[2],Ya[2],Wd[2];

	int lev=qp->lev;

	if(lev > 6){
		//qp->draw();
		(*count)++;
		return(lev);
	}

	int icrs=is_cross(qp->px, cr);
	//printf("Relation of G to A is %d\n",icrs);
	//translate_crs(icrs);
	Wd[0]=qp->px.Wd[0]*0.5;
	Wd[1]=qp->px.Wd[1]*0.5;
	Xa[0]=qp->px.Xa[0];
	Xa[1]=qp->px.Xa[1];
	if(icrs>1){
		k=0;
		for(j=0; j<2; j++){
			Ya[1]=Xa[1]+Wd[1]*j;
		for(i=0; i<2; i++){
			Ya[0]=Xa[0]+Wd[0]*i;
			qp->chld[k]=new_QPatch(Ya,Wd);
			qp->chld[k]->par=qp;
			qp->chld[k]->lev=lev+1;
			Qtree(qp->chld[k], cr, count);
			k++;
		}
		}
	}else{
		(*count)++;
		//qp->draw();
	};
	return(qp->lev);
};
int main(){
	int count;
	char fname[128]="qtree.out";
	char mode[3]="w";

	Circ cr;
	cr.xc[0]=1.0;
	cr.xc[1]=0.0;
	cr.radi=1.0;

	double xa[2],xb[2];
	xa[0]=0.0; xa[1]=0.0;
	xb[0]=2.0; xb[1]=1.0;

	QPatch qp0;	// root node
	qp0.set_lim(xa,xb);
	cr.draw(fname,180,mode);
	count=0;
	Qtree(&qp0, cr,&count);
	printf("number of leaves=%d\n",count);
	QPatch *qp_leaves=(QPatch *)malloc(sizeof(QPatch)*count);
	count=0;
	gather_leaves(&qp0,&count,qp_leaves);

	sprintf(mode,"a");
	for(int i=0;i<count;i++) qp_leaves[i].draw(fname,mode);

	return(0);
};
