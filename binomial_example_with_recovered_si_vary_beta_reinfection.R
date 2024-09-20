library(cmdstanr)
library(bayesplot)
seed <- sample(1e5,1)
seed <- 38376
#seed <- 77386
import_rate <- 
set.seed(seed)
TT <- 55
pop_size <- 2e5
I_0 <- 200
S_0 <- pop_size - I_0
ii <- RS <-  SI <- IR <- R <- S <- I <- rep(NA_real_,TT)
gamma <- 0.4
phi <- 0.15
p <- 0.9
S[1] <- S_0
RS[1] <- SI[1] <- 0
I[1] <- I_0
IR[1] <- R[1] <- 0 
ii[1] <- rbinom(1, I_0, 0.9)
beta <- 0.5 + (cos(0:(TT-1)/6) + 1) / 1.5
for (t in 2:TT) {
  mu <- 1 - exp(- beta[t] * I[t-1] / pop_size)
  kappa <- 100
  p_draw <- rbeta(1, mu * kappa, (1 - mu) * kappa)
  SI[t] <- rbinom(1, S[t-1], p_draw)
#  gamma_draw <- rbeta(1, kappa * gamma, (1 - gamma) * kappa)
  IR[t] <- rbinom(1, I[t-1], gamma)
  I[t] <- I[t-1] + SI[t] - IR[t]
#  phi_draw <- rbeta(1, kappa * phi, (1 - phi) * kappa)
  RS[t] <- rbinom(1, R[t-1], phi)
  S[t] <- S[t - 1] - SI[t] + RS[t]
  R[t] <- R[t-1] + IR[t] - RS[t]
  ii[t] <- rbinom(1, SI[t], p)
}
plot(ii,type="l")

plot(beta)

tt <- cmdstan_model("stoch-beta-test-uk-w-recover-logit-normal-si-vary-beta-BB-reinfection-all-cp.stan")
tt <- cmdstan_model("stoch-beta-test-uk-w-recover-logit-normal-si-vary-beta-BB-reinfection-all-cp-rho-par.stan")

dat <- 
  list(
    v_0 = S_0 / pop_size,
    ii = ii,
    TT = TT,
    pop_size = pop_size,
    b_freq = 5
  )

fit <- tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.80,
                 max_treedepth = 15,
                 init = \() {list(u_t_logit = qlogis(rbeta(dat$TT, 9, 9)),
                                  w_t_logit = qlogis(rbeta(dat$TT, 9, 9)),
                                  v_t_logit = qlogis(rbeta(dat$TT, 9, 9)),
                                  b = exp(rnorm(dat$TT %/% dat$b_freq  + ifelse(dat$TT %% dat$b_freq > 0,1,0)) * 0.5),
                                  p = rbeta(1, 4, 4),
                                  i0 = rbeta(1, 1, 10),
                                  gamma = rbeta(1, 4, 4),
                                  phi = rbeta(1, 2, 4),
                                  rho = rbeta(1, 1, 9),
                                  inv_kappa = rbeta(1, 0.01 * 1000, 0.99 * 1e3),
                                  sd_log_b = abs(rnorm(1)*0.3),
                                  mean_log_b = rnorm(1, 0.1,0.1))},
                 iter_warmup = 2000,
                 iter_sampling = 2000,
                 parallel_chains = 4)

opt <- tt$optimize(data = dat)

dag <- tt$diagnose(data = dat)

pars_sum_1 <- fit$sampler_diagnostics()
colSums(pars[,,"divergent__"])
colSums(pars[,,"treedepth__"] == 10)
colMeans(pars[,,"treedepth__"])
colMeans(pars_sum_1[,,"stepsize__"])

ess_sum_1 <- fit$summary()$ess_bulk |> summary()

np_fit <- nuts_params(fit)
mcmc_pairs(fit$draws(), np = np_fit, pars = c("b","p","gamma","u_t[7]","i_t[6]"),
           off_diag_args = list(size = 0.75))

hist(fit$draws("rho"))

b_draws <- fit$draws("b",format = "draws_matrix")
plot(colMeans(b_draws), ylim = range(c(
  apply(b_draws,2,quantile,0.05), apply(b_draws,2,quantile,0.95)
)))
lines(apply(b_draws,2,quantile,0.05),col="red",lty=2)
lines(apply(b_draws,2,quantile,0.95),col="red",lty=2)
lines(beta[seq(1,100,by=5)])
u_t_d <- fit$draws("si_t",format = "draws_matrix")
plot(SI/pop_size)
lines(colMeans(u_t_d),col="red")
lines(apply(u_t_d,2,quantile,0.05),col="red",lty=2)
lines(apply(u_t_d,2,quantile,0.95),col="red",lty=2)

u_t_d <- fit$draws("i_t",format = "draws_matrix")
plot(I)
lines(pop_size * colMeans(u_t_d),col="red")
lines(pop_size * apply(u_t_d,2,quantile,0.05),col="red",lty=2)
lines(pop_size * apply(u_t_d,2,quantile,0.95),col="red",lty=2)

u_t_d <- fit$draws("ir_t",format = "draws_matrix")
plot(IR)
lines(pop_size * colMeans(u_t_d),col="red")
lines(pop_size * apply(u_t_d,2,quantile,0.05),col="red",lty=2)
lines(pop_size * apply(u_t_d,2,quantile,0.95),col="red",lty=2)

u_t_d <- fit$draws("rs_t",format = "draws_matrix")
plot(RS)
lines(pop_size * colMeans(u_t_d),col="red")
lines(pop_size * apply(u_t_d,2,quantile,0.05),col="red",lty=2)
lines(pop_size * apply(u_t_d,2,quantile,0.95),col="red",lty=2)

u_t_d <- fit$draws("w_t",format = "draws_matrix")[,-1]
qs <- c(apply(u_t_d,2,quantile,0.05),apply(u_t_d,2,quantile,0.95))
plot(colMeans(u_t_d),col="red",ylim=range(qs))
lines(apply(u_t_d,2,quantile,0.05),col="red",lty=2)
lines(apply(u_t_d,2,quantile,0.95),col="red",lty=2)

fit$summary("gamma")
fit$summary("phi")

u_t_d <- fit$draws("w_t_logit_eta",format = "draws_matrix")
qs <- c(apply(u_t_d,2,quantile,0.05),apply(u_t_d,2,quantile,0.95))
plot(colMeans(u_t_d),col="red", ylim=range(qs))
lines(apply(u_t_d,2,quantile,0.05),col="red",lty=2)
lines(apply(u_t_d,2,quantile,0.95),col="red",lty=2)

u_t_d <- fit$draws("v_t_logit_eta",format = "draws_matrix")
qs <- c(apply(u_t_d,2,quantile,0.05),apply(u_t_d,2,quantile,0.95))
plot(colMeans(u_t_d),col="red", ylim=range(qs))
lines(apply(u_t_d,2,quantile,0.05),col="red",lty=2)
lines(apply(u_t_d,2,quantile,0.95),col="red",lty=2)


u_t_d <- fit$draws("u_t_logit",format = "draws_matrix")
qs <- c(apply(u_t_d,2,quantile,0.05),apply(u_t_d,2,quantile,0.95))
plot(colMeans(u_t_d),col="red",ylim=range(qs))
lines(apply(u_t_d,2,quantile,0.05),col="red",lty=2)
lines(apply(u_t_d,2,quantile,0.95),col="red",lty=2)

hist(u_t_d[,3])
vars <- apply(u_t_d,2,var)
jumps <- diff(I) / head(I,-1)
plot(jumps, tail(vars,-1))
length(vars)
length(jumps)
I[20:25]
jumps[20:25]

which.min(vars)

lims <- hist(fit$draws("p"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.5 * 3, 0.5 * 3),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("p"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = p,col="black")

lims <- hist(-log(1-fit$draws("gamma")),plot=FALSE)
ymax <- lims$density |> max()
hist(-log(1-rbeta(1e6,0.8 * 5, 0.2 * 5)),freq=FALSE,ylim=c(0,ymax+1), border = NA,breaks=100,
     main = "Prior (grey) vs. posterior (red) for rate of recovery",
     xlab = bquote(-log(1-gamma)))
hist(-log(1-fit$draws("gamma")),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = -log(1-gamma),col="black")

lims <- hist(-log(1-fit$draws("phi")),plot=FALSE)
ymax <- lims$density |> max()
hist(-log(1-rbeta(1e6,0.15 * 50, 0.85 * 50)),freq=FALSE,ylim=c(0,ymax+1), border = NA,breaks=100,
     main = "Prior (grey) vs. posterior (red) for rate of waning",
     xlab = bquote(-log(1-phi)))
hist(-log(1-fit$draws("phi")),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = -log(1-phi),col="black")


lims <- hist(fit$draws("i0"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.01 * 50, 0.99 * 50),freq=FALSE,ylim=c(0,ymax+1), border = NA,breaks=100,
     main = "Prior (grey) vs. posterior (red) for I_0",
     xlab = bquote(I[0]))
hist(fit$draws("i0"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = I_0/pop_size,col="black")

lims <- hist(fit$draws("kappa"),plot=FALSE)
ymax <- lims$density |> max()
hist(1/rbeta(1e6,0.01 * 100, 0.99 * 100),freq=FALSE,ylim=c(0,ymax+1), border = NA,breaks=100,
     main = "Prior (grey) vs. posterior (red) for kappa",
     xlim=c(0,1e3),
     xlab = bquote(kappa))
hist(fit$draws("kappa"),freq=FALSE,col=rgb(1,0,0,0.5),border = NA)
abline(v = kappa,col="black")

lims <- hist(fit$draws("p"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.8 * 50, 0.2 * 50),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("p"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = p,col="black")

hist(fit$draws("phi",format = "draws_array"),breaks=100)
abline(v = phi, col = "red")

hist(fit$draws("i0",format = "draws_array"),breaks=100)
abline(v = I_0/pop_size, col = "red")

hist(fit$draws("p",format = "draws_array"),breaks=100)
abline(v = 0.9, col = "red")

plot(fit$draws("p"), fit$draws("kappa"))
points(p, 100,col='red')

z_t_d <- fit$draws("i_t", format = "draws_matrix")
plot(I, xlab = "Time", ylab = "Prevalence", main = "Latent prevalence v. time")
lines(pop_size * colMeans(z_t_d), col = "red")
lines(pop_size * apply(z_t_d,2,quantile,0.05),col="red",lty=2)
lines(pop_size * apply(z_t_d,2,quantile,0.95),col="red",lty=2)

z_t_d <- fit$draws("si_t", format = "draws_matrix")
plot(SI, xlab = "Time", ylab = "Prevalence", main = "Latent prevalence v. time")
lines(pop_size * colMeans(z_t_d), col = "red")
lines(pop_size * apply(z_t_d,2,quantile,0.05),col="red",lty=2)
lines(pop_size * apply(z_t_d,2,quantile,0.95),col="red",lty=2)

z_t_d <- fit$draws("si_t", format = "draws_matrix") * 0.9
plot(ii, xlab = "Time", ylab = "Prevalence", main = "Latent prevalence v. time")
lines(pop_size * colMeans(z_t_d), col = "red")
lines(pop_size * apply(z_t_d,2,quantile,0.05),col="red",lty=2)
lines(pop_size * apply(z_t_d,2,quantile,0.95),col="red",lty=2)

z_t_d <- fit$draws("si_t", format = "draws_matrix") * 0.9
plot(, xlab = "Time", ylab = "Prevalence", main = "Latent prevalence v. time")
lines(pop_size * colMeans(z_t_d), col = "red")
lines(pop_size * apply(z_t_d,2,quantile,0.05),col="red",lty=2)
lines(pop_size * apply(z_t_d,2,quantile,0.95),col="red",lty=2)

lims <- hist(fit$draws("p"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.9 * 50, 0.1 * 50),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("p"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = 0.9,col="black")

lims <- hist(fit$draws("gamma"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.9 * 20, 0.1 * 20),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("gamma"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = 0.9,col="black")

lims <- hist(fit$draws("b"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,2, 0.75),freq=FALSE,ylim=c(0,ymax+1),breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for intra-county rate of infection in 1st cty",
     xlab = bquote(beta[11]))
hist(fit$draws("b"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = beta,col="black")

