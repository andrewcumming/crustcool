#define PI 3.14159265358979323846
// Numerical recipes vectors and matrices
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
static double minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
void nrerror(const char error_text[]);
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
