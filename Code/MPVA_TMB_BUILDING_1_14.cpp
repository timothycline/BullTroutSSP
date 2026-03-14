#include <TMB.hpp>

template<class Type>
bool isNA(Type x){ return R_IsNA(asDouble(x)); }

template<class Type>
Type objective_function<Type>::operator() () {
  using namespace density;

  // ==== DATA ====
  DATA_INTEGER(nCore);
  DATA_INTEGER(nPop);
  DATA_INTEGER(nYears);
  DATA_INTEGER(nSites);

  DATA_VECTOR(extent);        // length nPop
  DATA_VECTOR(pS);            // length nSites
  DATA_MATRIX(ReddCounts);    // nSites x nYears

  DATA_IVECTOR(site_pop);     // length nSites, 0-based
  DATA_IVECTOR(pop_core);     // length nPop, 0-based

  // ==== COVARIATE DATA ====
  // Static population-level
  DATA_VECTOR(Adfluvial);     // length nPop
  DATA_VECTOR(VB);            // length nPop, unconfined valley bottom
  DATA_VECTOR(Flood);         // length nPop, log W95 (avg days exceeding 95th pct flow)
  DATA_VECTOR(RDVB);          // length nPop, log road density in valley bottom
  // Time-varying population-level
  DATA_MATRIX(Flow);          // nPop x nYears
  DATA_MATRIX(Temp);          // nPop x nYears
  DATA_MATRIX(INV);           // nPop x nYears
  // Time-varying core-level (demeaned temporal deviations)
  DATA_MATRIX(LKT);           // nCore x nYears
  // Cross-system LKT mean per core area (informs b0r_pop prior)
  DATA_VECTOR(LKT_mean);      // length nCore

  // ==== COVARIATE FLAGS (0 = off, 1 = on) ====
  DATA_INTEGER(use_Adfluvial);
  DATA_INTEGER(use_VB);
  DATA_INTEGER(use_Flood);
  DATA_INTEGER(use_RD);
  DATA_INTEGER(use_Flow);
  DATA_INTEGER(use_Temp);
  DATA_INTEGER(use_INV);
  DATA_INTEGER(use_LKT);

  // ==== PARAMETERS ====
  PARAMETER(log_sigmaN);

  // Growth rate hierarchy
  PARAMETER(b0r_M);
  PARAMETER_VECTOR(b0r_pop);   // length nPop
  PARAMETER_VECTOR(b0r_core);  // length nCore
  PARAMETER(log_sigma_b0r_pop);
  PARAMETER(log_sigma_b0r_core);

  // Density dependence hierarchy
  PARAMETER(log_b0phi_M);
  PARAMETER_VECTOR(log_b0phi_pop);   // length nPop
  PARAMETER_VECTOR(log_b0phi_core);  // length nCore
  PARAMETER(log_sigma_b0phi_pop);
  PARAMETER(log_sigma_b0phi_core);

  // ==== STATIC COVARIATE PARAMETERS ====
  // Effect on b0r_pop prior mean (hyper-means for Adfluvial and VB remain global)
  PARAMETER(bAdf);
  PARAMETER(bVB);

  // Flood and road density: global scalars
  // (pop-level random slopes on static covariates are not identifiable
  //  when growth rate intercepts are also random at the same level)
  PARAMETER(bFlood);
  PARAMETER(bRD);

  // ==== HIERARCHICAL COVARIATE SLOPE PARAMETERS ====
  // Flow: population-level random slopes
  PARAMETER(bFlow_M);
  PARAMETER_VECTOR(bFlow_pop);   // length nPop
  PARAMETER(log_sigma_bFlow);

  // Temperature linear: population-level random slopes
  // Quadratic curvature (log_bTemp2) remains global
  PARAMETER(bTemp_M);
  PARAMETER_VECTOR(bTemp_pop);   // length nPop
  PARAMETER(log_sigma_bTemp);
  PARAMETER(log_bTemp2);

  // Invasives: population-level random slopes
  PARAMETER(bINV_M);
  PARAMETER_VECTOR(bINV_pop);    // length nPop
  PARAMETER(log_sigma_bINV);

  // Lake trout: core-level random slopes
  // bLKT_M also used for cross-system mean effect on b0r_pop prior
  PARAMETER(bLKT_M);
  PARAMETER_VECTOR(bLKT_core);   // length nCore
  PARAMETER(log_sigma_bLKT);

  // Latent states
  PARAMETER_VECTOR(logN0);    // length nPop
  PARAMETER_MATRIX(epsN);     // nPop x (nYears-1)

  // ==== TRANSFORMS ====
  Type sigmaN           = exp(log_sigmaN);
  Type sigma_b0r_pop    = exp(log_sigma_b0r_pop);
  Type sigma_b0r_core   = exp(log_sigma_b0r_core);
  Type sigma_b0phi_pop  = exp(log_sigma_b0phi_pop);
  Type sigma_b0phi_core = exp(log_sigma_b0phi_core);

  Type sigma_bFlow  = exp(log_sigma_bFlow);
  Type sigma_bTemp  = exp(log_sigma_bTemp);
  Type sigma_bINV   = exp(log_sigma_bINV);
  Type sigma_bLKT   = exp(log_sigma_bLKT);

  Type bTemp2 = -exp(log_bTemp2);

  // ==== OBJECTIVE ====
  Type nll = 0.0;

  // ==== PRIORS ====
  // Process error
  nll -= dnorm(log_sigmaN, Type(log(0.2)), Type(1.0), true);

  // Growth rate hierarchy
  nll -= dnorm(b0r_M,              Type(0.1),     Type(1.5), true);
  nll -= dnorm(log_sigma_b0r_pop,  Type(log(0.5)), Type(1.0), true);
  nll -= dnorm(log_sigma_b0r_core, Type(log(0.5)), Type(1.0), true);

  // Density dependence
  nll -= dnorm(log_b0phi_M,          Type(log(2.0)), Type(10.0), true);
  nll -= dnorm(log_sigma_b0phi_pop,  Type(log(0.3)), Type(1.0),  true);
  nll -= dnorm(log_sigma_b0phi_core, Type(log(0.3)), Type(1.0),  true);

  // Static covariate priors (global scalars)
  nll -= dnorm(bAdf,   Type(0.0), Type(1.0), true);
  nll -= dnorm(bVB,    Type(0.0), Type(1.0), true);
  nll -= dnorm(bFlood, Type(0.0), Type(1.0), true);
  nll -= dnorm(bRD,    Type(0.0), Type(1.0), true);

  // Covariate hyper-mean priors (weakly informative on z-scored covariates)
  nll -= dnorm(bFlow_M, Type(0.0), Type(1.0), true);
  nll -= dnorm(bTemp_M, Type(0.0), Type(1.0), true);
  nll -= dnorm(bINV_M,  Type(0.0), Type(1.0), true);
  nll -= dnorm(bLKT_M,  Type(0.0), Type(1.0), true);

  // Covariate slope SD priors
  nll -= dnorm(log_sigma_bFlow, Type(log(0.3)), Type(1.0), true);
  nll -= dnorm(log_sigma_bTemp, Type(log(0.3)), Type(1.0), true);
  nll -= dnorm(log_sigma_bINV,  Type(log(0.3)), Type(1.0), true);
  nll -= dnorm(log_sigma_bLKT,  Type(log(0.3)), Type(1.0), true);

  // Temperature quadratic (global)
  nll -= dnorm(log_bTemp2, Type(log(0.2)), Type(1.0), true);

  // Initial abundance
  for(int p = 0; p < nPop; p++){
    nll -= dnorm(logN0[p], Type(0.0), Type(5.0), true);
  }

  // ==== HIERARCHICAL RANDOM EFFECTS: GROWTH RATE ====

  // b0r_core: global -> core
  for(int c = 0; c < nCore; c++){
    nll -= dnorm(b0r_core[c], b0r_M, sigma_b0r_core, true);
  }

  // b0r_pop: core -> pop, with static covariates and bLKT_M * LKT_mean informing mean
  for(int p = 0; p < nPop; p++){
    int c = pop_core[p];
    Type mu_r = b0r_core[c];
    if(use_Adfluvial) mu_r += bAdf          * Adfluvial[p];
    if(use_VB)        mu_r += bVB           * VB[p];
    if(use_Flood)     mu_r += bFlood * Flood[p];
    if(use_RD)        mu_r += bRD    * RDVB[p];
    if(use_LKT)       mu_r += bLKT_M        * LKT_mean[c];
    nll -= dnorm(b0r_pop[p], mu_r, sigma_b0r_pop, true);
  }

  // ==== HIERARCHICAL RANDOM EFFECTS: DENSITY DEPENDENCE ====

  for(int c = 0; c < nCore; c++){
    nll -= dnorm(log_b0phi_core[c], log_b0phi_M, sigma_b0phi_core, true);
  }
  for(int p = 0; p < nPop; p++){
    int c = pop_core[p];
    nll -= dnorm(log_b0phi_pop[p], log_b0phi_core[c], sigma_b0phi_pop, true);
  }

  // ==== HIERARCHICAL RANDOM EFFECTS: COVARIATE SLOPES ====

  if(use_Flow){
    for(int p = 0; p < nPop; p++)
      nll -= dnorm(bFlow_pop[p], bFlow_M, sigma_bFlow, true);
  }
  if(use_Temp){
    for(int p = 0; p < nPop; p++)
      nll -= dnorm(bTemp_pop[p], bTemp_M, sigma_bTemp, true);
  }
  if(use_INV){
    for(int p = 0; p < nPop; p++)
      nll -= dnorm(bINV_pop[p], bINV_M, sigma_bINV, true);
  }
  if(use_LKT){
    for(int c = 0; c < nCore; c++)
      nll -= dnorm(bLKT_core[c], bLKT_M, sigma_bLKT, true);
  }
  // ==== LATENT STATE MODEL ====
  matrix<Type> logN(nPop, nYears);

  for(int p = 0; p < nPop; p++){
    int c = pop_core[p];
    logN(p, 0) = logN0[p];

    for(int t = 1; t < nYears; t++){
      Type R = b0r_pop[p];

      // Time-varying covariates — pop-level random slopes
      if(use_Flow) R += bFlow_pop[p] * Flow(p, t-1);
      if(use_Temp) R += bTemp_pop[p] * Temp(p, t-1) + bTemp2 * pow(Temp(p, t-1), Type(2.0));
      if(use_INV)  R += bINV_pop[p]  * INV(p,  t-1);

      // Core-level random slope for LKT (demeaned temporal deviations)
      if(use_LKT)  R += bLKT_core[c] * LKT(c, t-1);

      // Density dependence
      Type logdens = logN(p, t-1) - log(extent[p]);
      logdens = CppAD::CondExpGt(logdens, Type(10.0), Type(10.0), logdens);
      Type dens = exp(logdens) / Type(1000.0);
      R += -exp(log_b0phi_pop[p]) * dens;

      // Process error
      nll -= dnorm(epsN(p, t-1), Type(0.0), sigmaN, true);

      // State update
      logN(p, t) = logN(p, t-1) + R + epsN(p, t-1);

      if(!isfinite(logN(p, t))) nll += Type(1e10);
    }
  }

  // ==== OBSERVATION MODEL ====
  for(int s = 0; s < nSites; s++){
    int p = site_pop[s];
    for(int y = 0; y < nYears; y++){
      if(!isNA(ReddCounts(s, y))){
        Type N_pred = exp(logN(p, y));
        Type mu     = Type(0.5) * N_pred * pS[s];
        mu = CppAD::CondExpLt(mu, Type(1e-12), Type(1e-12), mu);
        nll -= dpois(ReddCounts(s, y), mu, true);
      }
    }
  }

  // ==== REPORTING ====
  ADREPORT(sigmaN);
  ADREPORT(b0r_M);
  ADREPORT(b0r_pop);
  ADREPORT(b0r_core);
  ADREPORT(sigma_b0r_pop);
  ADREPORT(sigma_b0r_core);
  ADREPORT(log_b0phi_M);
  ADREPORT(log_b0phi_pop);
  ADREPORT(log_b0phi_core);
  ADREPORT(sigma_b0phi_pop);
  ADREPORT(sigma_b0phi_core);
  ADREPORT(logN);
  ADREPORT(bAdf);
  ADREPORT(bVB);
  ADREPORT(bFlood);
  ADREPORT(bRD);
  ADREPORT(bFlow_M);   ADREPORT(bFlow_pop);   ADREPORT(sigma_bFlow);
  ADREPORT(bTemp_M);   ADREPORT(bTemp_pop);   ADREPORT(sigma_bTemp);
  ADREPORT(bTemp2);
  ADREPORT(bINV_M);    ADREPORT(bINV_pop);    ADREPORT(sigma_bINV);
  ADREPORT(bLKT_M);    ADREPORT(bLKT_core);   ADREPORT(sigma_bLKT);

  return nll;
}
