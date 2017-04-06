#include <Rcpp.h>

// functions

bool validKummerParameters(double a, double b, bool warn = true);
bool isInteger(double x, bool warn = true);
bool validProbability(double p);
bool isInadmissible(double x, bool warn = false);

// macros
#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
