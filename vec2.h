class Vec2{
	public:
	double x[2];
	void set(double,double);	
	void set(double*);	
	Vec2 times(double);	
	Vec2 div(double);	
	double len();
	void print();
	void rotate(double th);
	Vec2();
	Vec2(double xf[2]);
	private:
};
double iprod(Vec2 a,Vec2 b);
double iprod(double a[2], double b[2]);
Vec2 vsum(Vec2 a, Vec2 b);
Vec2 vdiff(Vec2 a, Vec2 b);
double vdist(Vec2 a, Vec2 b);
double distP2L(
	double xpt[2],
	double x1[2],
	double x2[2],
	double *rmin,
	double *rmax
);
