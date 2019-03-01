void mem_alloc2D(int nx, int ny, double ***pt);
void show_msg(char *fname);

//	********************************************
//
//		DOMAIN CLASS
//
//	********************************************

//-------------Dom2D Class ------------------
class Dom2D{
	public:
		double Xa[2],Xb[2],dx[2];// Computational Domain
		int Nx[2],Ndiv[2];
		int iprd[2];	// periodic B.C. (1:yes,2:No)
		int **kcell;
		int perfo(char *fn);
		int perfo_rec(char *fn);
		void perfo0(Circ cdat, int iphs);
		int mscs(char *fn);
		void out_kcell(char *fname);
		Dom2D();
		Dom2D(int ndiv1, int ndiv2);
		Dom2D(char *fname);
		void set_dx();
		int count(int i);
		double time;	// time (ps) in DEM simulation
		int draw_line(double x1[2], double x2[2],int iphs, int lw);
	private:
		void mem_alloc();
};

double dist2d(double *x1, double *x2);

