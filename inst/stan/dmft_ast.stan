// dmft_ast.stan — Random intercept mixed-effects model for DMFT/dmft
//
// Mirrors the lme4 formula:
//   y ~ 1 + [covariates] + (1|region) + (1|year)
// with inverse-variance weighting via known standard errors.
//
// AST smoothing is applied post-hoc in R (Stage 2).

data {
  int<lower=1> N;              // number of observations
  int<lower=1> N_regions;      // number of regions (provinces)
  int<lower=1> N_years;        // number of years
  int<lower=0> K;              // number of fixed-effect covariates (0 = intercept only)

  vector[N] y;                 // response: mean DMFT/dmft
  vector<lower=0>[N] se;       // known standard errors (for weighting)
  array[N] int<lower=1,upper=N_regions> region_idx;
  array[N] int<lower=1,upper=N_years>   year_idx;

  // Optional covariate matrix (N x K); ignored when K = 0
  matrix[N, K] X;
}

parameters {
  real alpha;                           // intercept
  vector[K] beta;                       // fixed-effect coefficients
  vector[N_regions] z_region;           // non-centered region effects
  vector[N_years]   z_year;             // non-centered year effects
  real<lower=0> sigma_region;           // SD of region random effects
  real<lower=0> sigma_year;             // SD of year random effects
  real<lower=0> sigma;                  // residual SD
}

transformed parameters {
  vector[N_regions] re_region = z_region * sigma_region;
  vector[N_years]   re_year   = z_year   * sigma_year;

  vector[N] mu;
  for (i in 1:N) {
    mu[i] = alpha + re_region[region_idx[i]] + re_year[year_idx[i]];
  }
  if (K > 0) {
    mu += X * beta;
  }
}

model {
  // Priors
  alpha ~ normal(5, 10);
  beta  ~ normal(0, 5);

  // Non-centered parameterization for random effects
  z_region ~ std_normal();
  z_year   ~ std_normal();

  // Weakly informative half-normal priors on SDs
  sigma_region ~ normal(0, 5);
  sigma_year   ~ normal(0, 5);
  sigma        ~ normal(0, 5);

  // Likelihood with observation-level variance
  // Total SD per observation = sqrt(sigma^2 + se_i^2)
  for (i in 1:N) {
    y[i] ~ normal(mu[i], sqrt(sigma^2 + se[i]^2));
  }
}

generated quantities {
  vector[N] y_rep;   // posterior predictive draws
  vector[N] log_lik; // pointwise log-likelihood (for LOO-CV)

  for (i in 1:N) {
    real total_sd = sqrt(sigma^2 + se[i]^2);
    y_rep[i]   = normal_rng(mu[i], total_sd);
    log_lik[i] = normal_lpdf(y[i] | mu[i], total_sd);
  }
}
