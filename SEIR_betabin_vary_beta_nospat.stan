data {
  int<lower=1> N_C;
  int<lower=0> TT;
  vector<lower=0>[N_C] pop_size;
  array[TT-1,N_C] real<lower=0> ii;
  matrix[N_C, N_C] D;
  int b_freq;
}

transformed data {
    int N_seg = TT %/% b_freq + (TT % b_freq > 0 ? 1 : 0);
}

parameters {

  matrix[TT, N_C] u_t_logit_eta;
  matrix[TT, N_C] w_t_logit_eta;
  matrix[TT, N_C] v_t_logit_eta;
  row_vector<lower=0, upper=1>[N_C] i0;
  row_vector<lower=0, upper=1>[N_C] e0;
  real<lower=0, upper=1> p;
  vector<lower=0, upper=1>[N_C] gamma;
  vector<lower=0, upper=1>[N_C] eta;

  real<lower=0> sigma; // Standard deviation of noise
  real<lower=0> rho; // Spatial range parameter
  real<lower=0,upper=1> rho_se; // Spatial range parameter
  real<lower=0,upper=1> rho_ir; // Spatial range parameter
  real<lower=0,upper=1> rho_ei; // Spatial range parameter
  real<lower=0> decay_rate_space; // Spatial decay rate
  vector[N_seg]  log_beta;
}

transformed parameters {
  matrix[TT, N_C] s_t;
  matrix[TT, N_C] se_t;
  matrix[TT, N_C] ei_t;
  matrix[TT, N_C] i_t;
  matrix[TT, N_C] e_t;
  matrix[TT, N_C] ir_t;
  matrix[TT, N_C] u_t;
  matrix[TT, N_C] v_t;
  matrix[TT, N_C] w_t;
  vector[N_seg] beta = exp(log_beta);



  i_t[1,] = i0;
  e_t[1,] = e0;
  s_t[1] = 1 - i0 - e0;
  se_t[1,] = rep_row_vector(0,N_C);
  ir_t[1,] = rep_row_vector(0,N_C);
  ei_t[1,] = rep_row_vector(0,N_C);
  for (n in 2:TT) {
    int n_idx = (n - 1) %/% b_freq + 1;

    for (ct in 1:N_C) {

      real u_t_mean = exponential_cdf(beta[n_idx]*i_t[n - 1,ct] | 1);
      u_t[n,ct] = inv_logit(u_t_logit_eta[n,ct]*(sqrt((1-rho_se) / (pop_size[ct] * s_t[n-1,ct] * u_t_mean * (1 - u_t_mean)) + rho_se/(u_t_mean*(1-u_t_mean)))) + logit(u_t_mean));
      v_t[n,ct] = inv_logit(v_t_logit_eta[n,ct]*(sqrt((1-rho_ei) / (pop_size[ct] * e_t[n-1,ct] * eta[ct] * (1 - eta[ct])) + rho_ei/(eta[ct]*(1-eta[ct])))) + logit(eta[ct]));
      w_t[n,ct] = inv_logit(w_t_logit_eta[n,ct]*(sqrt((1-rho_ir) / (pop_size[ct] * i_t[n-1,ct] * gamma[ct] * (1 - gamma[ct])) + rho_ir/(gamma[ct]*(1-gamma[ct])))) + logit(gamma[ct]));
      se_t[n,ct] = u_t[n,ct] * s_t[n-1,ct];
      ei_t[n,ct] = v_t[n,ct] * e_t[n-1,ct];
      ir_t[n,ct] = w_t[n,ct] * i_t[n-1,ct];
      s_t[n,ct] = s_t[n-1,ct] - se_t[n,ct];
      e_t[n,ct] = e_t[n-1,ct] + se_t[n,ct] - ei_t[n,ct];
      i_t[n,ct] = i_t[n-1,ct] + ei_t[n,ct] - ir_t[n,ct];
    }
  }
}
model {
  // Priors

  rho_ir ~ beta(1,3); //gamma(2, 2);
  rho_ei ~ beta(1,3); //gamma(2, 2);
  rho_se ~ beta(1,3);//gamma(2, 2);
  log_beta ~ normal(0, 0.5);
  // Prior for log_beta_diag

 i0 ~ beta_proportion(0.01, 50);
 e0 ~ beta_proportion(0.01, 50);

 p ~ beta_proportion(0.75, 5);
 gamma ~ beta_proportion(0.8, 15);
 eta ~ beta_proportion(0.8, 15);

 to_vector(u_t_logit_eta) ~ std_normal();
 to_vector(v_t_logit_eta) ~ std_normal();
 to_vector(w_t_logit_eta) ~ std_normal();

 for (i in 1:(TT-1)) {
  for (ct in 1:N_C) {
        if (se_t[i+1,ct] > 0)
          ii[i,ct] ~ normal(p * pop_size[ct] * ei_t[i+1,ct], sqrt(pop_size[ct] * p * ei_t[i+1,ct] * (1 - p)));
    }
  }
}

