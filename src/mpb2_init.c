#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP mpb2_chf_1F1(SEXP, SEXP, SEXP);
extern SEXP mpb2_cpp_dmpb(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpb2_cpp_pmpb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpb2_cpp_qmpb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mpb2_cpp_rmpb(SEXP, SEXP, SEXP, SEXP);
extern SEXP mpb2_gmRNA_switch(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"mpb2_chf_1F1", (DL_FUNC) &mpb2_chf_1F1, 3},
  {"mpb2_cpp_dmpb",    (DL_FUNC) &mpb2_cpp_dmpb,    5},
  {"mpb2_cpp_pmpb",    (DL_FUNC) &mpb2_cpp_pmpb,    6},
  {"mpb2_cpp_qmpb",    (DL_FUNC) &mpb2_cpp_qmpb,    6},
  {"mpb2_cpp_rmpb",    (DL_FUNC) &mpb2_cpp_rmpb,    4},
  {"mpb2_gmRNA_switch",       (DL_FUNC) &mpb2_gmRNA_switch,       5},
  {NULL, NULL, 0}
};

void R_init_mpb2(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
