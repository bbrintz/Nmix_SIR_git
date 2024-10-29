data {
  int<lower=1> N_C;
  int<lower=0> TT;
  vector<lower=0>[N_C] pop_size;
  array[TT-1,N_C] real<lower=0> ii;
}

parameters {

  matrix[TT, N_C] u_t_logit_eta;
  matrix[TT, N_C] w_t_logit_eta;
  row_vector<lower=0, upper=1>[N_C] i0;
  real<lower=0, upper=1> p;
  vector<lower=0, upper=1>[N_C] gamma;

  //real<lower=0,upper=1> rho_si; // Spatial range parameter
  real<lower=0,upper=1> rho_ir; // Spatial range parameter
  real<lower=-1,upper=1> phi;
  vector[TT] Z;
  real<lower=0> sigma; 
}

transformed parameters {
  matrix[TT, N_C] s_t;
  matrix[TT, N_C] si_t;
  matrix[TT, N_C] i_t;
  matrix[TT, N_C] ir_t;
  matrix[TT, N_C] u_t;
  matrix[TT, N_C] w_t;
  vector[TT] log_beta;
  vector<lower=0>[TT] beta;
  real<lower=0,upper=0> rho_si;
  rho_si = 0;
  i_t[1,] = i0;
  s_t[1] = 1 - i0;
  si_t[1,] = rep_row_vector(0,N_C);
  ir_t[1,] = rep_row_vector(0,N_C);
  log_beta[1] = Z[1];
  beta[1]=exp(log_beta[1]);
  for (n in 2:TT) {
    log_beta[n] = phi*log_beta[n-1] + Z[n];
    beta[n]=exp(log_beta[n]);
    for (ct in 1:N_C) {

      real u_t_mean = exponential_cdf(beta[n-1]*i_t[n - 1,ct] | 1);
      u_t[n,ct] = inv_logit(u_t_logit_eta[n,ct]*(sqrt((1-rho_si) / (pop_size[ct] * s_t[n-1,ct] * u_t_mean * (1 - u_t_mean)) + rho_si/(u_t_mean*(1-u_t_mean)))) + logit(u_t_mean));
      w_t[n,ct] = inv_logit(w_t_logit_eta[n,ct]*(sqrt((1-rho_ir) / (pop_size[ct] * i_t[n-1,ct] * gamma[ct] * (1 - gamma[ct])) + rho_ir/(gamma[ct]*(1-gamma[ct])))) + logit(gamma[ct]));
      si_t[n,ct] = u_t[n,ct] * s_t[n-1,ct];
      ir_t[n,ct] = w_t[n,ct] * i_t[n-1,ct];
      s_t[n,ct] = s_t[n-1,ct] - si_t[n,ct];
      i_t[n,ct] = i_t[n-1,ct] + si_t[n,ct] - ir_t[n,ct];
    }
  }
}
model {
  // Priors

  rho_ir ~ beta(1,3); //gamma(2, 2);
  //rho_si ~ beta(1,3);//gamma(2, 2);
  phi ~ uniform(-1,1);
  sigma ~ gamma(2,.5);
  Z[2:TT] ~ normal(0, sigma);
  Z[1] ~ normal(0, sigma / sqrt(1 - phi^2)); 
  i0 ~ beta_proportion(0.01, 50);

 p ~ beta_proportion(0.75, 5);
 gamma ~ beta_proportion(0.8, 15);

 to_vector(u_t_logit_eta) ~ std_normal();
 to_vector(w_t_logit_eta) ~ std_normal();

 for (i in 1:(TT-1)) {
  for (ct in 1:N_C) {
        if (si_t[i+1,ct] > 0)
          ii[i,ct] ~ normal(p * pop_size[ct] * si_t[i+1,ct], sqrt(pop_size[ct] * p * si_t[i+1,ct] * (1 - p)));
    }
  }
}

