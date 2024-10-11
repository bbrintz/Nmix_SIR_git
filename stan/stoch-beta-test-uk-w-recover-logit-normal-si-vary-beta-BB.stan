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
}
transformed data {
  array[TT] real prop_ii;
  int N_seg = TT %/% b_freq + (TT % b_freq > 0 ? 1 : 0);
  for (t in 1:TT)
    prop_ii[t] = ii[t] / pop_size;
}
parameters {
  vector<lower=0>[N_seg] b;
  vector[TT] u_t_logit_eta;
  vector[TT] w_t_logit_eta;
  real<lower=0, upper=1> i0;
  real<lower=0, upper=1> p;
  real<lower=0, upper=1> gamma;
  real<lower=0> inv_kappa;
}
transformed parameters {
  vector[TT] s_t;
  vector[TT] si_t;
  vector[TT] i_t;
  vector[TT] ir_t;
  vector[TT] u_t;
  vector[TT] w_t;
  real kappa = 1 / inv_kappa;
  si_t[1] = 0;
  i_t[1] = i0;
  s_t[1] = 1 - i0;
  ir_t[1] = 0;
  for (n in 2:TT) {
    int n_idx = (n - 1) %/% b_freq + 1;
    real u_t_mean = exponential_cdf(b[n_idx] * i_t[n - 1] | 1);
    u_t[n] = inv_logit(u_t_logit_eta[n] * sqrt(kappa + 1) / sqrt(pop_size * s_t[n-1] * u_t_mean * (1 - u_t_mean) * (kappa + pop_size * s_t[n-1])) + logit(u_t_mean));
    si_t[n] = u_t[n] * s_t[n-1];
    s_t[n] = s_t[n-1] - si_t[n];
    w_t[n] = inv_logit(w_t_logit_eta[n] / sqrt(pop_size * i_t[n-1] * gamma * (1 - gamma)) + logit(gamma));
    ir_t[n] = w_t[n] * i_t[n-1];
    i_t[n] = i_t[n-1] + si_t[n] - ir_t[n];
  }
}
model {
  b ~ normal(1, 0.5);
  i0 ~ beta_proportion(0.01, 50);
  p ~ beta_proportion(0.8, 50);
  gamma ~ beta_proportion(0.8, 20);
  inv_kappa ~ normal(0, 10);
  u_t_logit_eta ~ std_normal();
  w_t_logit_eta ~ std_normal();
  for (i in 2:TT) {
    ii[i] ~ normal(p * pop_size * si_t[i], sqrt(pop_size * si_t[i] * p * (1 - p)));// * (kappa + pop_size * si_t[i])/(kappa + 1)));
  }
}

