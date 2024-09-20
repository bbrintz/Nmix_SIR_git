functions {
  real manual_binomial_lpdf(real x, real N, real p) {
    return lgamma(N + 1) - lgamma(N - x + 1) + x * log(p) + (N - x) * log1m(p);
  }
}
data {
  int<lower=0> TT;
  real<lower=0> pop_size;
  array[TT] real<lower=0> ii;
  int b_freq;
  // array[TT] real b;
}
transformed data {
  array[TT] real prop_ii;
  int N_seg = TT %/% b_freq + (TT % b_freq > 0 ? 1 : 0);
  // real mean_log_b = 0.1;
  // real<lower=0> sd_log_b = 0.2;
  // real<lower=0, upper=1> rho = 0.9;
  // real kappa = 100;
  // real gamma = 0.6;
  // real phi = 0.15;
  // real p = 0.9;
  for (t in 1:TT)
    prop_ii[t] = ii[t] / pop_size;
}
parameters {
  vector[TT] u_t_logit_eta;
  vector[TT] w_t_logit_eta;
  vector[TT] v_t_logit_eta;
  vector<lower=0>[N_seg] b;
  real<lower=0, upper=1> i0;
  real<lower=0, upper=1> p;
 real<lower=0, upper=1> gamma;
 real<lower=0, upper=1> phi;
 // real mean_log_b;
 // real<lower=0> sd_log_b;
 // real<lower=0, upper=1> rho;
 real<lower=0, upper=1> inv_kappa;
}
transformed parameters {
  vector[TT] s_t;
  vector[TT] i_t;
  vector[TT] r_t;
  vector[TT] si_t;
  vector[TT] ir_t;
  vector[TT] rs_t;
  vector[TT] u_t;
  vector[TT] w_t;
  vector[TT] v_t;
  // vector[TT] b;
  real kappa = 1 / inv_kappa;
 // b[1] = exp(mean_log_b / (1 - rho) + log_b[1] * sd_log_b / sqrt((1 - rho)*(1+rho)));
 // for (t in 2:N_seg)
 //   b[t] = exp(mean_log_b + log(b[t-1]) * rho + sd_log_b * log_b[t]);
  si_t[1] = 0;
  i_t[1] = i0;
  s_t[1] = 1 - i0;
  ir_t[1] = 0;
  rs_t[1] = 0;
  r_t[1] = 0;
  for (n in 2:TT) {
    int n_idx = (n - 1) %/% b_freq + 1;
    real u_t_mean = exponential_cdf(b[n_idx] * i_t[n - 1] | 1);
    u_t[n] = inv_logit(u_t_logit_eta[n] * sqrt(kappa + pop_size * s_t[n-1]) / sqrt(pop_size * s_t[n-1] * u_t_mean * (1 - u_t_mean) * (kappa + 1)) + logit(u_t_mean));
    w_t[n] = inv_logit(w_t_logit_eta[n] / sqrt(pop_size * i_t[n-1] * gamma * (1 - gamma)) + logit(gamma));
    if (r_t[n-1] > 0)
      v_t[n] = inv_logit(v_t_logit_eta[n] / sqrt(pop_size * r_t[n-1] * phi * (1 - phi)) + logit(phi));
    else 
      v_t[n] = 0;
    si_t[n] = u_t[n] * s_t[n-1];
    ir_t[n] = w_t[n] * i_t[n-1];
    rs_t[n] = v_t[n] * r_t[n-1];
    s_t[n] = s_t[n-1] - si_t[n] + rs_t[n];
    i_t[n] = i_t[n-1] + si_t[n] - ir_t[n];
    r_t[n] = r_t[n-1] + ir_t[n] - rs_t[n];
  }
}
model {
 // mean_log_b ~ normal(0.1, 0.1);
 // sd_log_b ~ normal(0, 5);
 // log_b ~ std_normal();
 b ~ normal(1, 0.5);
 i0 ~ beta_proportion(0.01, 50);
 p ~ beta_proportion(0.5, 3);
 gamma ~ beta_proportion(0.8, 5);
 phi ~ beta_proportion(0.05, 50);
 inv_kappa ~ beta_proportion(0.001, 100);
  u_t_logit_eta ~ std_normal();
  w_t_logit_eta ~ std_normal();
  v_t_logit_eta ~ std_normal();
  for (i in 2:TT) {
    ii[i] ~ normal(p * pop_size * si_t[i], sqrt(pop_size * si_t[i] * p * (1 - p)));// * (kappa + pop_size * si_t[i])/(kappa + 1)));
  }
}

