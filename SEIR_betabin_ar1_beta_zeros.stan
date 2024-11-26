data {
  int<lower=1> N_C;
  int<lower=0> TT;
  vector<lower=0>[N_C] pop_size;
  array[TT-1,N_C] real<lower=0> ii;
  array[N_C] int<lower=0> first;

  }

parameters {

  matrix[TT, N_C] u_t_logit_eta;
  matrix[TT, N_C] w_t_logit_eta;
  matrix[TT, N_C] v_t_logit_eta;
  row_vector<lower=0, upper=1>[N_C] i0;
  real<lower=0, upper=1> p;
  vector<lower=0, upper=1>[N_C] gamma;
  vector<lower=0, upper=1>[N_C] eta;
  //real<lower=0,upper=1> rho_ir; // Spatial range parameter
  //real<lower=0,upper=1> rho_ei; // Spatial range parameter
  real<lower=-1,upper=1> phi;
  vector[TT] Z;
  real<lower=0> sigma; 
}

transformed parameters {
  matrix[TT, N_C] s_t;
  matrix[TT, N_C] se_t;
  matrix[TT, N_C] ei_t;
  matrix[TT, N_C] i_t;
  matrix[TT, N_C] e_t = rep_matrix(0,TT,N_C);
  matrix[TT, N_C] ir_t;
  matrix[TT, N_C] u_t;
  matrix[TT, N_C] v_t;
  matrix[TT, N_C] w_t;
  vector[TT] log_beta;
  vector<lower=0>[TT] beta;
  real<lower=0,upper=0> rho_se;
  real<lower=0,upper=0> rho_ei;
  real<lower=0,upper=0> rho_ir;
  real u_t_mean;

  rho_se = 0;
  rho_ei = 0;
  rho_ir = 0;

  se_t=rep_matrix(0,TT,N_C);
  ir_t=rep_matrix(0,TT,N_C);
  ei_t=rep_matrix(0,TT,N_C);


  for (ct in 1:N_C){
  i_t[1:(first[ct]-1),ct] = rep_vector(0,first[ct]-1);
  i_t[first[ct],ct]=i0[ct];

  s_t[1:(first[ct]-1),ct] = rep_vector(1,first[ct]-1);
  s_t[first[ct],ct] = 1-i0[ct];

  ei_t[first[ct],ct] = i0[ct];
  }


 

  log_beta[1] = Z[1];
  beta[1]=exp(log_beta[1]);

  for (n in 2:TT) {
        log_beta[n] = phi * log_beta[n - 1] + Z[n];
        beta[n] = exp(log_beta[n]);
  }



  for (ct in 1:N_C) {
    for (n in 2:TT) {
      if (n > first[ct]) {
            // Begin conditional block
            //print("n: ", n, "ct: ", ct, "i_t[n-1,ct]: ", i_t[n-1,ct], "e_t[n-1,ct]: ", e_t[n-1,ct], "s_t[n-1,ct]: ", s_t[n-1,ct]);
            u_t_mean = exponential_cdf(beta[n - 1] * i_t[n - 1, ct] | 1);
            u_t[n-1, ct] = inv_logit(u_t_logit_eta[n-1, ct] * (sqrt((1 - rho_se) / (pop_size[ct] * s_t[n - 1, ct] * u_t_mean * (1 - u_t_mean)) + rho_se / (u_t_mean * (1 - u_t_mean)))) + logit(u_t_mean));

            if (ei_t[n-1, ct] == 0) {
              v_t[n-1, ct] = inv_logit(v_t_logit_eta[n-1, ct] * (sqrt((1 - rho_ei) / (pop_size[ct] * e_t[n - 1, ct] * eta[ct] * (1 - eta[ct])) + rho_ei / (eta[ct] * (1 - eta[ct])))) + logit(eta[ct]));
              ei_t[n-1, ct] = v_t[n-1, ct] * e_t[n - 1, ct];
            }

            w_t[n-1, ct] = inv_logit(w_t_logit_eta[n-1, ct] * (sqrt((1 - rho_ir) / (pop_size[ct] * i_t[n - 1, ct] * gamma[ct] * (1 - gamma[ct])) + rho_ir / (gamma[ct] * (1 - gamma[ct])))) + logit(gamma[ct]));
            se_t[n-1, ct] = u_t[n-1, ct] * s_t[n - 1, ct];
            
            ir_t[n-1, ct] = w_t[n-1, ct] * i_t[n - 1, ct];
            s_t[n, ct] = s_t[n - 1, ct] - se_t[n-1, ct];
            if (n > first[ct]+1) {
              e_t[n, ct] = e_t[n - 1, ct] + se_t[n-1, ct] - ei_t[n-1, ct];
            } else {
              e_t[n, ct] = e_t[n - 1, ct] + se_t[n-1, ct];
            }
            i_t[n, ct] = i_t[n - 1, ct] + ei_t[n-1, ct] - ir_t[n-1, ct];
          } 
      }
    }
  }

model {
  // Priors
  //rho_se ~ beta(1,3);//gamma(2, 2);
  //rho_ir ~ beta(1,3); //gamma(2, 2);
  //rho_ei ~ beta(1,3); //gamma(2, 2);
  phi ~ uniform(-1,1);
  sigma ~ gamma(2, 2);
  Z[2:TT] ~ normal(0, sigma);
  Z[1] ~ normal(0, sigma / sqrt(1 - phi^2)); 


 i0 ~ beta_proportion(0.01, 5);

 p ~ beta_proportion(0.25, 5);
 gamma ~ beta_proportion(0.8, 15);
 eta ~ beta_proportion(0.9, 15);

 to_vector(u_t_logit_eta) ~ std_normal();
 to_vector(v_t_logit_eta) ~ std_normal();
 to_vector(w_t_logit_eta) ~ std_normal();

 for (ct in 1:N_C) {
  for (i in 1:(TT-1)) {
       //   print("ii[i,ct]: ", ii[i,ct],"ei_t[i,ct]: ", ei_t[i,ct], "i: ", i, "ct: ", ct);
    if (i >= first[ct])
      ii[i,ct] ~ normal(p * pop_size[ct] * ei_t[i,ct], sqrt(pop_size[ct] * p * ei_t[i,ct] * (1 - p)));
    }
  }
}
