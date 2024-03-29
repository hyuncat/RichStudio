// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RcppKappaCluster2
List RcppKappaCluster2(CharacterVector myTerms, CharacterVector myGenes, NumericVector myPvalues, double kappaCutoff);
RcppExport SEXP _RichStudio_RcppKappaCluster2(SEXP myTermsSEXP, SEXP myGenesSEXP, SEXP myPvaluesSEXP, SEXP kappaCutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type myTerms(myTermsSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type myGenes(myGenesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type myPvalues(myPvaluesSEXP);
    Rcpp::traits::input_parameter< double >::type kappaCutoff(kappaCutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppKappaCluster2(myTerms, myGenes, myPvalues, kappaCutoff));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello
List rcpp_hello();
RcppExport SEXP _RichStudio_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RichStudio_RcppKappaCluster2", (DL_FUNC) &_RichStudio_RcppKappaCluster2, 4},
    {"_RichStudio_rcpp_hello", (DL_FUNC) &_RichStudio_rcpp_hello, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_RichStudio(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
