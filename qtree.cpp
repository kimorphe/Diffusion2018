#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "set2d.h"

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
