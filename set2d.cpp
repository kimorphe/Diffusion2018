#define DB 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec2.h"
#include "set2d.h"

using namespace std;
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
//---------------------------------------------------------------
class Ellip{
	public:
		double xc[2];	// center
		double radi[2]; // radii
		double phi;	// angle [rad]
		void draw(char fn[128], int np, char mode[3]);
		void draw(int np);
		void set_xc(double x,double y);
		void set_radi(double r1,double r2);
		void set_phi(double ang);
		bool is_in(double xf[2]);
		Ellip();
	private:
	protected:
};

Ellip::Ellip(){
	xc[0]=0.0; xc[1]=0.0;
	radi[0]=1.0; radi[1]=1.0;
	phi=0.0;
};
void Ellip::draw(int np){
	double dth=8.0*atan(1.0)/np;
	double th;
	double x,y,X,Y;
	double cosp=cos(phi),sinp=sin(phi);
	for(int i=0;i<=np;i++){
		th=i*dth;
		x=radi[0]*cos(th);
		y=radi[1]*sin(th);
		X= cosp*x+sinp*y+xc[0];
		Y=-sinp*x+cosp*y+xc[1];
	       	printf("%lf %lf\n",X,Y);
	}
};
void Ellip::draw(char fn[128], int np, char mode[3]){
	FILE *fp=fopen(fn,mode);
	double dth=8.0*atan(1.0)/np;
	double th;
	double x,y,X,Y;
	double cosp=cos(phi),sinp=sin(phi);
	for(int i=0;i<=np;i++){
		th=i*dth;
		x=radi[0]*cos(th);
		y=radi[1]*sin(th);
		X= cosp*x+sinp*y+xc[0];
		Y=-sinp*x+cosp*y+xc[1];
	       	fprintf(fp,"%lf %lf\n",X,Y);
	}

	fclose(fp);
};
void Ellip::set_xc(double x, double y){
	xc[0]=x;
	xc[1]=y;
};
void Ellip::set_radi(double r1, double r2){
	radi[0]=r1;
	radi[0]=r2;
};
void Ellip::set_phi(double ang){
	double pi=4.0*atan(1.0);
	phi=ang/180.*pi;
};
bool Ellip::is_in(double xf[2]){
	double yf[2],Y[2];

	yf[0]=xf[0]-xc[0];
	yf[1]=xf[1]-xc[1];

	double cosp=cos(phi),sinp=sin(phi);
	Y[0]=cosp*yf[0]-sinp*yf[1];
	Y[1]=sinp*yf[0]+cosp*yf[1];
	Y[0]/=radi[0];
	Y[1]/=radi[1];

	bool iin=true;

	if(Y[0]*Y[0]+Y[1]*Y[1]>1.0) iin=false;
	return(iin);
};
class Poly{
	public:
		int np;	
		double *xs,*ys;
		double xg[2];	// centroide
		void draw();
		void draw(char fn[128],char mode[3]);
		void mem_alloc();
		Poly();
		Poly(int n);
		bool is_in(double xf[2]);
		void slide(double ux, double uy);
		void rotate(double xc[2], double th);
		void rotate(int node_num, double th);
		void set_center();
	private:
	protected:
};
Poly::Poly(){};
Poly::Poly(int n){
	np=n;
	Poly::mem_alloc();
	xg[0]=0.0;
	xg[1]=0.0;
};
void Poly::mem_alloc(){
	xs=(double *)std::malloc(sizeof(double)*np);
	ys=(double *)std::malloc(sizeof(double)*np);
	xg[0]=0.0;
	xg[1]=0.0;
};
void Poly::set_center(){
	for(int i=0;i<np;i++){
		xg[0]+=xs[i];
		xg[1]+=ys[i];
	}
	xg[0]/=np;
	xg[1]/=np;
}
void Poly::draw(){
	for(int i=0;i<=np;i++){
		printf("%lf %lf\n",xs[i%np],ys[i%np]);
	};
	printf("\n");
};
void Poly::draw(char fn[128],char mode[3]){
	FILE *fp=fopen(fn,mode);
	for(int i=0;i<=np;i++){
		fprintf(fp,"%lf %lf\n",xs[i%np],ys[i%np]);
	};
	fprintf(fp,"\n");
	fclose(fp);
};
bool Poly::is_in(double xf[2]){
	int i;
	Vec2 x1,x2,yf,tx,nx;
	yf.set(xf[0],xf[1]);
	for(i=0;i<np;i++){
		x1.set(xs[i],ys[i]);	
		x2.set(xs[(i+1)%np],ys[(i+1)%np]);	
		tx=vdiff(x2,x1);
		nx.set(tx.x[1],-tx.x[0]);
		if(iprod(nx, vdiff(yf,x1))>0.0){
			return(false);
		}
	};
	return(true);
};
void Poly::slide(double ux, double uy){
	for(int i=0;i<np;i++){
		xs[i]+=ux;
		ys[i]+=uy;
	}
};
void Poly::rotate(int node_num, double th){

	Vec2 xc,ri;
	if(node_num<0){
		xc.set(xg);	// centroide
	}else{
		xc.set(xs[node_num%np],ys[node_num%np]);
	};
		

	for(int i=0; i<np; i++){
		ri.set(xs[i],ys[i]);
		ri=vdiff(ri,xc);
		ri.rotate(th);
		ri=vsum(ri,xc);
		xs[i]=ri.x[0];
		ys[i]=ri.x[1];
	}
};
void Poly::rotate(double xr0[2], double th){
	Vec2 xc(xr0);
	Vec2 ri;

	for(int i=0; i<np; i++){
		ri.set(xs[i],ys[i]);
		ri=vdiff(ri,xc);
		ri.rotate(th);
		ri=vsum(ri,xc);
		xs[i]=ri.x[0];
		ys[i]=ri.x[1];
	}
};

int poly_cross(Poly A, Poly B){
	int i,j;
	double xf[2];

	int intr=0;
	int extr=0;
	for(i=0;i<A.np;i++){
		xf[0]=A.xs[i];
		xf[1]=A.ys[i];
		if(B.is_in(xf)) intr++;
	};

	if(intr==A.np) return(2); // A < B

	int jntr=0;
	for(i=0;i<B.np;i++){
		xf[0]=B.xs[i];
		xf[1]=B.ys[i];
		if(A.is_in(xf)) jntr++;
	};

	if(jntr==B.np) return(1); // B < A
	if((intr+jntr)==0) return(0); // A^B=\phi

	return(3); // A^B != \phi
};
#if DB==2
int main(){
	Poly plA(5);
	Poly plB(4);

	plA.xs[0]=0.0; plA.ys[0]=0.0;
	plA.xs[1]=1.0; plA.ys[1]=0.0;
	plA.xs[2]=1.0; plA.ys[2]=1.0;
	plA.xs[3]=0.5; plA.ys[3]=1.5;
	plA.xs[4]=0.0; plA.ys[4]=0.8;
	plA.set_center();


	plB.xs[0]=0.5; plB.ys[0]=0.0;
	plB.xs[1]=0.0; plB.ys[1]=0.5;
	plB.xs[2]=-0.5; plB.ys[2]=0.0;
	plB.xs[3]=-0.5; plB.ys[3]=-0.5;
	plB.set_center();

	double xf[2];

	char fname[128]="temp.dat";
	char md[3]="w";
	plA.draw(fname,md);
	plB.slide(-0.3,-0.3);

	sprintf(md,"a");
	double pi=4.*atan(1.0);
	double dth=10.0/180.*pi;
	xf[0]=0.0;
	xf[1]=0.0;
	for(int i=0;i<10;i++){
		//plB.slide(0.1,0.1);

		plA.draw(fname,md);
		plB.draw(fname,md);

		printf("icrs=%d\n",poly_cross(plA,plB));
		//plA.rotate(xf,dth);
		plA.rotate(1,dth);
	}
	return(0);
};
#endif
#if DB==1
int main(){
	Ellip el;
	el.set_xc(1.0,0.5);
	el.set_radi(2.0,0.5);
	el.set_phi(130.0);
	char fn[128]="temp.dat";
	char md[3]="w";
	el.draw(fn,100,md);

	double xf[2];
	double dx=0.1;
	double pi=4.0*atan(1.0);
	double th=15.0/180.*pi;
	for(int i=0;i<50;i++){
		xf[0]=(dx*i)*cos(th);	
		xf[1]=(dx*i)*sin(th);	
		printf("%lf %lf %d\n",xf[0],xf[1],el.is_in(xf));
	};	
	return(0);
};
#endif
