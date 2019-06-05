#include<stdio.h>
#include<stdlib.h>
#include<random>
#include"set2d.h"
#include"vec2.h"

int main(){

	printf("start program\n");
	Solid sld;

	//int np=100;	// number of particles
	int Lev=9;	// Quad-tree height
	double Wd[2]={1.0,1.0}; // Unit Cell Size
	double Xa[2]={0.0,0.0}; // Unit Cell position (lowerleft vertex)
	char fn[128]="solid.dat";	// solid phase data file

	puts(fn);
	sld.load(fn);	// import solid phase data

	sld.bbox.setup(Xa,Wd); // set bounding box

	Tree4 tr4;	// Tree structure 
	tr4.qp0.refine[0]=true;
	tr4.qp0.show_prms();
	tr4.setup(sld.els,sld.nelp,false,Lev,sld.bbox);
	//tr4.write();
	tr4.draw();
	printf("Tree data written.\n");
	return(0);
};
