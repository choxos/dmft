// GBD-Style DMFT/dmft Bayesian Hierarchical Model
// Equivalent to INLA BYM2 + RW2 + RW1 specification
// For dental caries burden estimation in Canada

functions {
  // Soft sum-to-zero constraint for random effects
  real soft_zero_sum_lpdf(vector x, real sigma) {
    return normal_lpdf(sum(x) | 0, sigma * 0.001 * rows(x));
  }
}

data {
  int<lower=1> N;                          // Number of observations
  int<lower=1> N_prov;                     // Number of provinces (13)
  int<lower=1> N_year;                     // Number of years
  int<lower=1> N_age;                      // Number of age groups
  int<lower=1> N_edges;                    // Number of adjacency edges

  // Response (can be continuous mean or count)
  int<lower=0> use_negbin;                 // 1 = negative binomial, 0 = Gaussian
  array[N] real y_continuous;              // For Gaussian likelihood
  array[N] int<lower=0> y_count;           // For NegBin likelihood

  // Sample sizes (for weighting or as offset)
  array[N] real<lower=0> sample_size;

  // Standard errors (for Gaussian model)
  array[N] real<lower=0> se;

  // Index variables
  array[N] int<lower=1, upper=N_prov> prov_idx;
  array[N] int<lower=1, upper=N_year> year_idx;
  array[N] int<lower=1, upper=N_age> age_idx;
  array[N] int<lower=1, upper=2> sex_idx;  // 1 = Male, 2 = Female

  // Spatial adjacency structure (for ICAR)
  array[N_edges] int<lower=1, upper=N_prov> node1;
  array[N_edges] int<lower=1, upper=N_prov> node2;

  // Covariates (optional)
  int<lower=0> K;                          // Number of covariates
  matrix[N, K] X;                          // Covariate matrix

  // Scaling factor for ICAR (for variance interpretation)
  real<lower=0> scaling_factor;
}

transformed data {
  // Pre-compute for efficiency
  vector[N] log_n = log(to_vector(sample_size));
}

parameters {
  // Intercept
  real alpha;

  // Covariate effects
  vector[K] beta;

  // BYM2 spatial effects (Riebler et al. 2016)
  vector[N_prov] phi_unstructured;         // Unstructured (IID) component
  vector[N_prov] psi_structured;           // Spatially structured (ICAR) component
  real<lower=0, upper=1> rho_spatial;      // Mixing parameter (proportion spatial)
  real<lower=0> sigma_spatial;             // Marginal SD for spatial effect

  // Temporal effects (RW2)
  vector[N_year] gamma_raw;                // Raw temporal effects
  real<lower=0> sigma_temporal;            // SD for temporal innovations

  // Age effects (RW1)
  vector[N_age] delta_raw;                 // Raw age effects
  real<lower=0> sigma_age;                 // SD for age innovations

  // Sex effects (IID)
  vector[2] eta;
  real<lower=0> sigma_sex;

  // Interaction effects (optional - can be turned off)
  vector[N_prov * N_year] st_effect;       // Space-time
  real<lower=0> sigma_st;

  // Overdispersion (for NegBin)
  real<lower=0> reciprocal_phi;            // 1/phi for NegBin
}

transformed parameters {
  // Combined spatial effect (BYM2)
  vector[N_prov] S;

  // Temporal effect with RW2 structure (cumulative)
  vector[N_year] gamma;

  // Age effect with RW1 structure (cumulative)
  vector[N_age] delta;

  // Linear predictor
  vector[N] mu;

  // BYM2: convex combination of structured and unstructured
  // S = sqrt(rho) * psi_scaled + sqrt(1-rho) * phi_scaled
  S = sqrt(rho_spatial / scaling_factor) * psi_structured +
      sqrt(1 - rho_spatial) * phi_unstructured;
  S = sigma_spatial * S;

  // Construct RW2 temporal effects (consistent scaling throughout)
  gamma[1] = sigma_temporal * gamma_raw[1];
  gamma[2] = sigma_temporal * gamma_raw[2];
  for (t in 3:N_year) {
    gamma[t] = 2 * gamma[t-1] - gamma[t-2] + sigma_temporal * gamma_raw[t];
  }

  // Construct RW1 age effects (consistent scaling throughout)
  delta[1] = sigma_age * delta_raw[1];
  for (a in 2:N_age) {
    delta[a] = delta[a-1] + sigma_age * delta_raw[a];
  }

  // Build linear predictor
  mu = rep_vector(alpha, N);

  // Add covariates
  if (K > 0) {
    mu += X * beta;
  }

  // Add random effects
  for (i in 1:N) {
    mu[i] += S[prov_idx[i]] +
             gamma[year_idx[i]] +
             delta[age_idx[i]] +
             eta[sex_idx[i]];

    // Add space-time interaction
    int st_index = (prov_idx[i] - 1) * N_year + year_idx[i];
    mu[i] += st_effect[st_index];
  }
}

model {
  // Priors

  // Intercept
  alpha ~ normal(0, 5);

  // Covariate effects
  beta ~ normal(0, 2);

  // BYM2 spatial priors
  // Unstructured component
  phi_unstructured ~ normal(0, 1);

  // ICAR prior for structured component
  target += -0.5 * dot_self(psi_structured[node1] - psi_structured[node2]);
  // Soft sum-to-zero constraint
  sum(psi_structured) ~ normal(0, 0.001 * N_prov);

  // PC priors for spatial parameters (Riebler et al. 2016)
  // P(sigma > 1) = 0.01 => exponential with rate -log(0.01)/1
  sigma_spatial ~ exponential(4.6);        // -log(0.01) ≈ 4.6

  // rho_spatial ~ beta prior centered at 0.5
  rho_spatial ~ beta(1, 1);

  // Temporal RW2 priors
  gamma_raw ~ normal(0, 1);
  sigma_temporal ~ exponential(4.6);

  // Age RW1 priors
  delta_raw ~ normal(0, 1);
  sigma_age ~ exponential(4.6);

  // Sex effects
  eta ~ normal(0, sigma_sex);
  sigma_sex ~ exponential(1);

  // Space-time interaction
  st_effect ~ normal(0, sigma_st);
  sigma_st ~ exponential(1);

  // Overdispersion
  reciprocal_phi ~ exponential(1);

  // Likelihood
  if (use_negbin == 1) {
    // Negative binomial for count data
    y_count ~ neg_binomial_2_log(mu + log_n, 1.0 / reciprocal_phi);
  } else {
    // Gaussian for continuous means (weighted by SE)
    for (i in 1:N) {
      y_continuous[i] ~ normal(mu[i], se[i]);
    }
  }
}

generated quantities {
  // Posterior predictive
  array[N] real y_rep;

  // Log-likelihood for LOO-CV
  array[N] real log_lik;

  // Fitted values (on response scale)
  array[N] real fitted;

  // Generate predictions and log-likelihood
  for (i in 1:N) {
    if (use_negbin == 1) {
      fitted[i] = exp(mu[i] + log_n[i]);
      y_rep[i] = neg_binomial_2_log_rng(mu[i] + log_n[i], 1.0 / reciprocal_phi);
      log_lik[i] = neg_binomial_2_log_lpmf(y_count[i] | mu[i] + log_n[i],
                                           1.0 / reciprocal_phi);
    } else {
      // Gaussian meta-regression: mu[i] is directly the predicted mean
      // (identity link, no transformation needed)
      fitted[i] = fmax(0.0, mu[i]);
      y_rep[i] = normal_rng(mu[i], se[i]);
      log_lik[i] = normal_lpdf(y_continuous[i] | mu[i], se[i]);
    }
  }

  // Province-specific estimates (for a reference year/age/sex)
  array[N_prov] real province_effect = to_array_1d(S);

  // Year-specific estimates
  array[N_year] real year_effect = to_array_1d(gamma);

  // Age-specific estimates
  array[N_age] real age_effect = to_array_1d(delta);
}
