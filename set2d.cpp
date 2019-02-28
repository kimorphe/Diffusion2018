#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "domain.h"


class Pixel{
	public:
		double Xa[2], Xb[2];
		double Wd[2],Xc[2];
		double xs[4],ys[4];
		double Le[4];
		double rmax,rmin;
		void print();
		void draw(char fn[128],char mode[3]);
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
//---------------------------------------------------
int main(){
	
	double rmin,x,y,dx,dy;
	Pixel px;
	px.set_Xa(0.0,0.0);
	px.set_Xb(2.0,1.0);
	//px.print();


	double x1=2.0,x2=2.0;
	double y1= -1.0,y2=2.0;
	int i,ndivx=1,ndivy=100;
	dx=0.0; dy=0.0;
	dx=(x2-x1)/ndivx;
	dy=(y2-y1)/ndivy;
	//for(i=0;i<=ndivx;i++){
	for(i=0;i<=ndivy;i++){
		x=x1+dx*i;
		y=y1+dy*i;
		rmin=px.dist2pt(x,y);
		printf("%lf %lf %lf\n",x,y,rmin);
	};

	Circ cr;
	cr.xc[0]=-1.00;
	cr.xc[1]=-0.5;
	cr.radi=1.0;

	char fnc[128]="circ.dat";
	char fnp[128]="pix.dat";
	char mode[3]="w";

	cr.draw(fnc,100,mode);
	px.draw(fnp,mode);
	int icrs;

	icrs=is_cross(px,cr);
	printf("icrs=%d\n",icrs);


	return(0);
};
