#define DB 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "set2d.h"

/*
class Solid{
	public:
		int nelp;
		Ellip *els;
		Solid();
		Solid(int n);
		void draw(char fn[128],int ndat);
	private:
	protected:
};
*/
Solid::Solid(){};
Solid::Solid(int n){
	nelp=n;
	els=(Ellip *)malloc(sizeof(Ellip)*nelp);
};
void Solid::draw(char fn[128],int ndat){
	char md[3];
	for(int i=0;i<nelp;i++){
		sprintf(md,"a");
		if(i==0) sprintf(md,"w");
		els[i].draw(fn,ndat,md);
	}
};

int main(){
	Solid SLD(3);

	SLD.els[0].xc[0]=0.0;
	SLD.els[0].xc[0]=0.0;
	SLD.els[0].phi=10.0;
	SLD.els[0].radi[0]=1.0;
	SLD.els[0].radi[1]=0.6;

	SLD.els[1]=SLD.els[0];
	SLD.els[2]=SLD.els[0];

	SLD.els[1].xc[1]=2.0;
	SLD.els[1].set_phi(85.0);
	SLD.els[1].radi[0]=2.0;

	SLD.els[2].xc[1]=4.0;

	char fn[128]="log.txt";

	double xa[2],xb[2];
	xa[0]=-1.0; xa[1]=-1.0;
	xb[0]= 1.0; xb[1]=5.0;
	QPatch qp0;	// root node
	qp0.set_lim(xa,xb);

	int count;
	count=0;
	Qtree(&qp0, SLD,&count);
	printf("number of leaves=%d\n",count);

	char fname[128]="qtree.out";
	char mode[3]="w";
	SLD.draw(fname,100);

	QPatch *qp_leaves=(QPatch *)malloc(sizeof(QPatch)*count);
	count=0;
	gather_leaves(&qp0,&count,qp_leaves);
	sprintf(mode,"a");
	for(int i=0;i<count;i++) qp_leaves[i].draw(fname,mode);
	return(0);
};
#if DB==1
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
#endif
