#include <Rcpp.h>

// fortran access

#ifdef __cplusplus
extern"C"{
#endif
  void chfm_(double *zre, double *zim,
             double *are, double *aim,
             double *bre, double *bim,
             double *cre, double *cim,
             int *n, int *lnchf, int *ip);
#ifdef __cplusplus
}
#endif

// constants

#define Q_LIMIT 64
#define INPUT_VECTORISED 1
#define INPUT_SINGLE 2

// functions

bool validKummerParameters(double a, double b, bool warn = true);
bool isInteger(double x, bool warn = true);
bool validProbability(double p, bool warn = false);
bool isInadmissible(double x, bool warn = false);
bool validMpbParameters(double alpha, double beta, double c, bool warn = false);
void reportGslError(int status);

// macros
#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
