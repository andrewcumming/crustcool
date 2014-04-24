class Crust;

class Data {
public:
	double *t, *TT, *Te;
	int n;
	int luminosity;
	
	void read_in_data(const char *fname);
	void calculate_chisq(Crust &crust);
};



