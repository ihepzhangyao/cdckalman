


#include "KalFitAlg/KalFitPar.h"
using namespace std;

int     KalFitPar::muls_=1;
int     KalFitPar::loss_=1;
int     KalFitPar::wsag_=99;
int     KalFitPar::back_=1;
double  KalFitPar::pT_=0.0;
int     KalFitPar::lead_=2;
int     KalFitPar::lr_=1;
double  KalFitPar::matrixg_=100.0;
int     KalFitPar::debug_=0;
int     KalFitPar::debug_kft_=0;
int     KalFitPar::ntuple_=51;
int     KalFitPar::activeonly_=0;
string  KalFitPar::matfile_;
string  KalFitPar::cylfile_;
double  KalFitPar::deltachi2_cutf_=10.0;
double  KalFitPar::deltachi2_cuts_=10.0;
double  KalFitPar::pt_=0.0;
double  KalFitPar::pe_cut_=2.0;
double  KalFitPar::pmu_cut_=2.0;
double  KalFitPar::ppi_cut_=2.0;
double  KalFitPar::pk_cut_=2.0;
double  KalFitPar::pp_cut_=2.0;
int     KalFitPar::i_back_=1;
int     KalFitPar::tofflag_=1;
int     KalFitPar::tof_hyp_=1;
int     KalFitPar::resolution_=0;
double  KalFitPar::fstrag_=0.9;
int     KalFitPar::enhance_=0;
double  KalFitPar::fac_h1_=1.0;
double  KalFitPar::fac_h2_=1.0;
double  KalFitPar::fac_h3_=1.0;
double  KalFitPar::fac_h4_=1.0;
double  KalFitPar::fac_h5_=1.0;
int     KalFitPar::numf_=0;
double  KalFitPar::Bznom_=10.0;
int     KalFitPar::numfcor_=1;
int     KalFitPar::steplev_=0;
int     KalFitPar::Tof_correc_=1;
int     KalFitPar::strag_=1;
double  KalFitPar::factor_strag_=0.4;
int     KalFitPar::tofall_=1;
int     KalFitPar::nmdc_hit2_=500;
double  KalFitPar::chi2mdc_hit2_=0;
int     KalFitPar::LR_=1;
int     KalFitPar::resolflag_=0;

void KalFitPar::lead(int i){ lead_ = i;}
int  KalFitPar::lead(void){ return lead_;}
                                                                          
void KalFitPar::back(int i){ back_ = i;}
int  KalFitPar::back(void){ return back_;}
                                                                          
void KalFitPar::resol(int i){ resolflag_ = i;}
int  KalFitPar::resol(void){ return resolflag_;}
                                                                          
void KalFitPar::LR(int x){ LR_ = x;}

void KalFitPar::loss(int i){ loss_ = i;}
                                                                          
int  KalFitPar::muls(void){ return muls_;}
                                                                          
int  KalFitPar::loss(void){  return loss_;}

void KalFitPar::muls(int i){muls_ = i;}
