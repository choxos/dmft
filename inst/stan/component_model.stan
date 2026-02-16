// Component Model for D/M/F (Decayed, Missing, Filled) Analysis
// Uses Isometric Log-Ratio (ILR) transformation for compositional data
// Equivalent to INLA model for component proportions

functions {
  // ILR transformation for 3-part composition
  // Transforms (D, M, F) to 2D ILR coordinates
  vector ilr_transform(vector comp) {
    vector[2] ilr;
    // ilr1 = sqrt(1/2) * log(D/M)
    ilr[1] = sqrt(0.5) * log(comp[1] / comp[2]);
    // ilr2 = sqrt(2/3) * log(geometric_mean(D,M) / F)
    ilr[2] = sqrt(2.0/3.0) * log(sqrt(comp[1] * comp[2]) / comp[3]);
    return ilr;
  }

  // Inverse ILR transformation
  // Transforms ILR coordinates back to composition
  vector ilr_inverse(vector ilr) {
    vector[3] comp;
    real x1 = ilr[1];
    real x2 = ilr[2];

    // Back-transform
    real log_D = x1 / sqrt(2.0) + x2 / sqrt(6.0);
    real log_M = -x1 / sqrt(2.0) + x2 / sqrt(6.0);
    real log_F = -2.0 * x2 / sqrt(6.0);

    // Normalize to sum to 1
    real log_sum = log_sum_exp(log_sum_exp(log_D, log_M), log_F);
    comp[1] = exp(log_D - log_sum);
    comp[2] = exp(log_M - log_sum);
    comp[3] = exp(log_F - log_sum);

    return comp;
  }
}

data {
  int<lower=1> N;                          // Number of observations
  int<lower=1> N_prov;                     // Number of provinces
  int<lower=1> N_year;                     // Number of years
  int<lower=1> N_age;                      // Number of age groups
  int<lower=1> N_edges;                    // Number of adjacency edges

  // ILR-transformed component data (2 coordinates per observation)
  array[N] vector[2] y_ilr;

  // Standard errors for ILR coordinates
  array[N] vector<lower=0>[2] se_ilr;

  // Index variables
  array[N] int<lower=1, upper=N_prov> prov_idx;
  array[N] int<lower=1, upper=N_year> year_idx;
  array[N] int<lower=1, upper=N_age> age_idx;
  array[N] int<lower=1, upper=2> sex_idx;

  // Spatial adjacency
  array[N_edges] int<lower=1, upper=N_prov> node1;
  array[N_edges] int<lower=1, upper=N_prov> node2;

  // Scaling factor for ICAR
  real<lower=0> scaling_factor;
}

parameters {
  // Intercepts for each ILR coordinate
  vector[2] alpha;

  // BYM2 spatial effects (for each ILR coordinate)
  array[2] vector[N_prov] phi_unstructured;
  array[2] vector[N_prov] psi_structured;
  vector<lower=0, upper=1>[2] rho_spatial;
  vector<lower=0>[2] sigma_spatial;

  // Temporal effects RW2
  array[2] vector[N_year] gamma_raw;
  vector<lower=0>[2] sigma_temporal;

  // Age effects RW1
  array[2] vector[N_age] delta_raw;
  vector<lower=0>[2] sigma_age;

  // Sex effects
  array[2] vector[2] eta;
  vector<lower=0>[2] sigma_sex;
}

transformed parameters {
  // Combined spatial effects
  array[2] vector[N_prov] S;

  // Temporal effects
  array[2] vector[N_year] gamma;

  // Age effects
  array[2] vector[N_age] delta;

  // Linear predictors (2D for ILR)
  array[N] vector[2] mu;

  // Build effects for each ILR coordinate
  for (c in 1:2) {
    // BYM2 spatial
    S[c] = sqrt(rho_spatial[c] / scaling_factor) * psi_structured[c] +
           sqrt(1 - rho_spatial[c]) * phi_unstructured[c];
    S[c] = sigma_spatial[c] * S[c];

    // RW2 temporal (consistent scaling throughout)
    gamma[c][1] = sigma_temporal[c] * gamma_raw[c][1];
    gamma[c][2] = sigma_temporal[c] * gamma_raw[c][2];
    for (t in 3:N_year) {
      gamma[c][t] = 2 * gamma[c][t-1] - gamma[c][t-2] +
                    sigma_temporal[c] * gamma_raw[c][t];
    }

    // RW1 age (consistent scaling throughout)
    delta[c][1] = sigma_age[c] * delta_raw[c][1];
    for (a in 2:N_age) {
      delta[c][a] = delta[c][a-1] + sigma_age[c] * delta_raw[c][a];
    }
  }

  // Build linear predictors
  for (i in 1:N) {
    for (c in 1:2) {
      mu[i][c] = alpha[c] +
                 S[c][prov_idx[i]] +
                 gamma[c][year_idx[i]] +
                 delta[c][age_idx[i]] +
                 eta[c][sex_idx[i]];
    }
  }
}

model {
  // Priors
  alpha ~ normal(0, 2);

  // BYM2 priors for each coordinate
  for (c in 1:2) {
    phi_unstructured[c] ~ normal(0, 1);
    target += -0.5 * dot_self(psi_structured[c][node1] - psi_structured[c][node2]);
    sum(psi_structured[c]) ~ normal(0, 0.001 * N_prov);

    sigma_spatial[c] ~ exponential(4.6);
    rho_spatial[c] ~ beta(1, 1);

    gamma_raw[c] ~ normal(0, 1);
    sigma_temporal[c] ~ exponential(4.6);

    delta_raw[c] ~ normal(0, 1);
    sigma_age[c] ~ exponential(4.6);

    eta[c] ~ normal(0, sigma_sex[c]);
    sigma_sex[c] ~ exponential(1);
  }

  // Likelihood (multivariate normal for ILR coordinates)
  for (i in 1:N) {
    // Independent Gaussians for each ILR coordinate
    for (c in 1:2) {
      y_ilr[i][c] ~ normal(mu[i][c], se_ilr[i][c]);
    }
  }
}

generated quantities {
  // Posterior predictive in ILR space
  array[N] vector[2] y_ilr_rep;

  // Back-transformed compositions
  array[N] vector[3] composition_fitted;
  array[N] vector[3] composition_rep;

  // Log-likelihood
  array[N] real log_lik;

  for (i in 1:N) {
    // Generate predictions
    for (c in 1:2) {
      y_ilr_rep[i][c] = normal_rng(mu[i][c], se_ilr[i][c]);
    }

    // Back-transform to compositions
    composition_fitted[i] = ilr_inverse(mu[i]);
    composition_rep[i] = ilr_inverse(y_ilr_rep[i]);

    // Log-likelihood
    log_lik[i] = normal_lpdf(y_ilr[i][1] | mu[i][1], se_ilr[i][1]) +
                 normal_lpdf(y_ilr[i][2] | mu[i][2], se_ilr[i][2]);
  }

  // Province-level compositions (averaged over other dimensions)
  array[N_prov] vector[3] province_composition;
  for (p in 1:N_prov) {
    vector[2] prov_ilr;
    for (c in 1:2) {
      prov_ilr[c] = alpha[c] + S[c][p];
    }
    province_composition[p] = ilr_inverse(prov_ilr);
  }

  // Age-specific compositions
  array[N_age] vector[3] age_composition;
  for (a in 1:N_age) {
    vector[2] age_ilr;
    for (c in 1:2) {
      age_ilr[c] = alpha[c] + delta[c][a];
    }
    age_composition[a] = ilr_inverse(age_ilr);
  }
}
