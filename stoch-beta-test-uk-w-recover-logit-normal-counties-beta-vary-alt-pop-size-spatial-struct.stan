functions {
  real manual_binomial_lpdf(real x, real N, real p) {
    return lgamma(N + 1) - lgamma(N - x + 1) + x * log(p) + (N - x) * log1m(p);
  }
}
data {
  int<lower=1> N_C;
  int<lower=0> TT;
  vector<lower=0>[N_C] pop_size;
  array[TT,N_C] real<lower=0> ii;
}
parameters {
  vector<lower=0>[N_C] b_self;
  vector<lower=0>[N_C - 2] b_od;
  matrix[TT, N_C] u_t_logit_eta;
  matrix[TT, N_C] w_t_logit_eta;
  row_vector<lower=0, upper=1>[N_C] i0;
  real<lower=0, upper=1> p;
  vector<lower=0, upper=1>[N_C] gamma;
}
transformed parameters {
  matrix[TT, N_C] s_t;
  matrix[TT, N_C] si_t;
  matrix[TT, N_C] i_t;
  matrix[TT, N_C] ir_t;
  matrix[TT, N_C] u_t;
  matrix[TT, N_C] w_t;
  matrix[N_C, N_C] b;
  for (n in 1:N_C)
    b[n, n] = b_self[n];
  for (i in 1:(N_C-1))
    for (j in (i + 1):N_C) {
      b[i,j] = j - i < 3 ? b_od[j - i] : 0;
      b[j,i] = b[i,j];
    }
  si_t[1,] = rep_row_vector(0,N_C);
  i_t[1,] = i0;
  s_t[1] = 1 - i0;
  ir_t[1,] = rep_row_vector(0,N_C);
  for (n in 2:TT) {
    for (ct in 1:N_C) {
      // if (is_nan(b[ct] * i_t[n - 1,ct])) {
      //   print("n: ",n, ", ct: ",ct, ", b: ",b[ct], ", i: ",i_t[1:(n - 1),ct]);
      //   print("gamma: ",gamma[ct], ", i0: ",i0[ct],", s0: ",1 - i0[ct], ", u_t: ",u_t_logit_eta[1:n,ct], ", w_t: ",w_t_logit_eta[1:n,ct]);
      //   print("lp: ", exponential_cdf(b[ct] * i_t[n - 2,ct] | 1));
      //   print("logit: ",logit(exponential_cdf(b[ct] * i_t[n - 2,ct] | 1)));
      // }
      real u_t_mean = exponential_cdf(dot_product(b[ct,], i_t[n - 1,]) | 1);
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
 b_self ~ gamma(4, 4);
 b_od ~ gamma(4, 4);
 i0 ~ beta_proportion(0.01, 50);
 p ~ beta_proportion(0.75, 5);
// gamma ~ beta_proportion(0.9, 15);
 gamma ~ beta_proportion(0.8, 15);
 to_vector(u_t_logit_eta) ~ std_normal();
 to_vector(w_t_logit_eta) ~ std_normal();
  for (i in 1:TT) {
    for (ct in 1:N_C) {
      // if (sqrt(pop_size[ct] * p * i_t[i,ct] * (1 - p * i_t[i,ct])) == 0)
        // print("i: ",i,", ct: ",ct,", p: ", p, ", i_t: ",i_t[i,ct]);
        if (i_t[i,ct] > 0)
          ii[i,ct] ~ normal(p * pop_size[ct] * i_t[i,ct], sqrt(pop_size[ct] * p * i_t[i,ct] * (1 - p * i_t[i,ct])));
    }
  }
}

