//---------------------------------------------------------------------------------
class Material{
	public:
		double gmm[3][3];	// interfacial energy (0:gas, 1:fluid, 2: solid)
		double thE;		// contact angle
		void load(double th); 	// set gmm values 
		double gmm_max;		// maximum value in gmm
		void  normalize();	// normalize gmm
		void  print_gmm(); // print interfacial energy
		double erg;	// cell energy
};
class Cell{
	public:
		QPatch *qp0;
		int cnct[8];//connected ?
		Cell *cncl[8]; //pointer to connected cells
		int nc;	// number  of connected cells
		Cell();	// constructuor
		bool bnd; // boundary cell (T/F)
		int ID;		// linear index for 2D Grid
		int iad;	// address in PoreCells[ncell]
		int phs;	// 0=gas, 1=fluid, 2=solid 
		int phs_bff;
		double erg;
		double erg_bff;
	private:
};
class PoreCells:public Tree4{
	public:
		int ncell; 	// number of cells
		Cell *cells;	// pointer(array) to cell class instances
		PoreCells();	// constructor
		void setup(Ellip *els,int nelp,bool set_opr, int Lev_Max,Bbox bx); // generate cells 
		int find(int id);	// find cell having a given linear grid number(id).
		void connect();	// establish neghboring cell connection 
		void l2ij(int l, int *i, int *j); // index transform ( linear to 2D index)
		double Sr;	// degree of water saturation
		int n_void,n_water; // number of gas and fluid cells, resp.
		int init(double sr); // initialize phase distribution
		int *indx_w;	// index set of fluid phase
		int *indx_v;	// index set of gas phase
		double cell_energy(int iad); // evaluate interfacial energy/cell
		double total_energy(); // evaluate total interfacial energy
		Material mtrl;	// material constants (interfaceial energy)
		void load_gmm(double th); // load gmm data (th = fluid/solid contact angle)
		double Etot;	// total interfacial energy
		double swap(int id, int jd); // swap cell id & jd tempralily
		void reject_swap(int id, int jd); // apply swap 
		double MC_stepping(Temp_Hist TH);	// Monte Carlo stepping 
		void write_phs();
		int count_grids();
		int grid_type(int i,int j);
		int grid_loc(int i,int j);
		void fwrite_cells(char fn[128]);
		void load_cell_data(char fn[128]);
		void write_leaves();
	private:
};
//---------------------------------------------------------------------------------
