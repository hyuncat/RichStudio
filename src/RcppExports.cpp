// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calKappa
double calKappa(const std::vector<std::string> t1_genes, const std::vector<std::string> t2_genes, double totalGeneCount);
RcppExport SEXP _RichStudio_calKappa(SEXP t1_genesSEXP, SEXP t2_genesSEXP, SEXP totalGeneCountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string> >::type t1_genes(t1_genesSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string> >::type t2_genes(t2_genesSEXP);
    Rcpp::traits::input_parameter< double >::type totalGeneCount(totalGeneCountSEXP);
    rcpp_result_gen = Rcpp::wrap(calKappa(t1_genes, t2_genes, totalGeneCount));
    return rcpp_result_gen;
END_RCPP
}
// RcppKappaCluster
List RcppKappaCluster(CharacterVector myTerms, CharacterVector myGenes, double kappaCutoff);
RcppExport SEXP _RichStudio_RcppKappaCluster(SEXP myTermsSEXP, SEXP myGenesSEXP, SEXP kappaCutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type myTerms(myTermsSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type myGenes(myGenesSEXP);
    Rcpp::traits::input_parameter< double >::type kappaCutoff(kappaCutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(RcppKappaCluster(myTerms, myGenes, kappaCutoff));
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
    {"_RichStudio_calKappa", (DL_FUNC) &_RichStudio_calKappa, 3},
    {"_RichStudio_RcppKappaCluster", (DL_FUNC) &_RichStudio_RcppKappaCluster, 3},
    {"_RichStudio_rcpp_hello", (DL_FUNC) &_RichStudio_rcpp_hello, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_RichStudio(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}