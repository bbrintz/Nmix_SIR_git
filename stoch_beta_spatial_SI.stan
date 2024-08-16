data {
  int<lower=1> N_C;
  int<lower=0> TT;
  vector<lower=0>[N_C] pop_size;
  array[TT-1,N_C] real<lower=0> ii;
  matrix[N_C, N_C] D;
}
parameters {

  matrix[TT, N_C] u_t_logit_eta;
  matrix[TT, N_C] w_t_logit_eta;
  row_vector<lower=0, upper=1>[N_C] i0;
  real<lower=0, upper=1> p;
  vector<lower=0, upper=1>[N_C] gamma;
  real<lower=0> sigma; // Standard deviation of noise
  real<lower=0> rho; // Spatial range parameter
  real<lower=0> decay_rate_space; // Spatial decay rate
  vector[N_C] log_beta_diag; // Diagonal elements of log_beta
}

transformed parameters {
  matrix[TT, N_C] s_t;
  matrix[TT, N_C] si_t;
  matrix[TT, N_C] i_t;
  matrix[TT, N_C] ir_t;
  matrix[TT, N_C] u_t;
  matrix[TT, N_C] w_t;
  matrix[N_C, N_C] Sigma;
  matrix[N_C, N_C] log_beta;
  matrix[N_C, N_C] beta;

   // Create spatial covariance matrix
  for (i in 1:N_C) {
    for (j in 1:N_C) {
      Sigma[i, j] = square(sigma) * exp(-D[i, j] / rho);
    }
  }

  // Initialize log_beta with zeros
  log_beta = rep_matrix(0, N_C, N_C);

  // Set diagonal elements
  for (i in 1:N_C) {
    log_beta[i, i] = log_beta_diag[i];
  }

  // Set off-diagonal elements based on the spatial decay
  for (i in 1:(N_C - 1)) {
    for (j in (i + 1):N_C) {
      //real distance_from_center = abs(i - j);
      log_beta[i, j] = log_beta[i, i] - decay_rate_space * D[i, j];
      log_beta[j, i] = log_beta[i, j];
    }
  }

  // Exponentiate log_beta to get beta
  beta = exp(log_beta);


  si_t[1,] = rep_row_vector(0,N_C);
  i_t[1,] = i0;
  s_t[1] = 1 - i0;
  ir_t[1,] = rep_row_vector(0,N_C);
  for (n in 2:TT) {
    for (ct in 1:N_C) {

      real u_t_mean = exponential_cdf(dot_product(beta[ct,], i_t[n - 1,]) | 1);
      u_t[n,ct] = inv_logit(u_t_logit_eta[n,ct] / sqrt(pop_size[ct] * s_t[n-1,ct] * u_t_mean * (1 - u_t_mean)) + logit(u_t_mean));
      w_t[n,ct] = inv_logit(w_t_logit_eta[n,ct] / sqrt(pop_size[ct] * i_t[n-1,ct] * gamma[ct] * (1 - gamma[ct])) + logit(gamma[ct]));
      si_t[n,ct] = u_t[n,ct] * s_t[n-1,ct];
      s_t[n,ct] = s_t[n-1,ct] - si_t[n,ct];
      ir_t[n,ct] = w_t[n,ct] * i_t[n-1,ct];
      i_t[n,ct] = i_t[n-1,ct] + si_t[n,ct] - ir_t[n,ct];
    }
  }
}
model {
  // Priors
  sigma ~ gamma(2, 2);
  rho ~ gamma(2, 2);
  decay_rate_space ~ gamma(2, 2);

  // Prior for log_beta_diag
  log_beta_diag ~ multi_normal(rep_vector(0, N_C), Sigma);

 i0 ~ beta_proportion(0.01, 50);
 p ~ beta_proportion(0.75, 5);
 gamma ~ beta_proportion(0.8, 15);
 to_vector(u_t_logit_eta) ~ std_normal();
 to_vector(w_t_logit_eta) ~ std_normal();
  for (ct in 1:N_C) {
    for (i in 1:(TT-1)) {
        if (si_t[i,ct] > 0)
          ii[i,ct] ~ normal(p * pop_size[ct] * si_t[i+1,ct], sqrt(pop_size[ct] * p * si_t[i+1,ct] * (1 - p)));
    }
  }
}

