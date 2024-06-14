library(cmdstanr)
library(bayesplot)
library(posterior)
#set.seed(123)
#beta <- 1.05
#TT <- 100
N_C <- 10 
set.seed(31)
beta <- diag(runif(N_C, min = 0.6, max = 0.9))
beta_od <- c(0.1, 0.05)
for (i in 1:(N_C-1))
  for (j in (i + 1):N_C) {
    beta[i,j] = ifelse(j - i < 3, beta_od[j - i], 0)
    beta[j,i] = beta[i,j]
  }
gamma <- runif(N_C, min = 0.7, max = 0.8)
TT <- 30
pop_size <- 1e3 * sample(1:10,N_C,TRUE)
I_0 <- sample(10:20,N_C,TRUE)
S_0 <- pop_size - I_0
ii <- R <- S <- I <- matrix(NA_real_,TT, N_C)
p_detect <- 0.6
S[1,] <- S_0
I[1,] <- I_0
R[1,] <- 0 
ii[1,] <- rbinom(N_C, I_0, p_detect)
for (t in 2:TT) {
  for (ct in  1:N_C) {
    if (t > 2)
      prev_t <- colSums(I[1:(t-1),]) - colSums(R[1:(t-1),])
    else
      prev_t <- I[1,] - R[1,]
    I[t,ct] <- rbinom(1, S[t-1,ct], 1 - exp(- sum(beta[ct,] * prev_t / pop_size)))
    S[t,ct] <- S[t - 1,ct] - I[t,ct]
    R[t,ct] <- rbinom(1, prev_t[ct], gamma[ct])
    ii[t,ct] <- rbinom(1, prev_t[ct] + I[t,ct] - R[t,ct], p_detect)
  }
}

matplot(ii, type="l")

tt <- cmdstan_model("stoch-beta-test-uk-w-recover-logit-normal-counties-beta-vary-alt-pop-size-spatial-struct.stan")
head((apply(I,2,cumsum) - apply(R,2,cumsum)))
tail((apply(I,2,cumsum) - apply(R,2,cumsum)))
(pop_size - S[TT,]) / pop_size

dat <- 
  list(
    ii = ii,
    TT = TT,
    N_C = N_C,
    pop_size = pop_size
  )

fit <- tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.95,
                 max_treedepth = 14,
                 init = \() {list(u_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  w_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  p = rbeta(1, 4, 4),
                                  b_od = rgamma(N_C-2, 4, 4),
                                  b_self = rgamma(N_C, 4, 4),
                                  i0 = rbeta(N_C, 0.01*50, 0.99*50),
                                  gamma = rbeta(N_C, 0.7 * 6, 0.3 * 6))},
                 iter_warmup = 1000,
                 iter_sampling = 1000, parallel_chains = 4,
                 step_size = 1.5e-3)


lims <- hist(fit$draws("p"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.75 * 5, 0.25 * 5),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("p"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = 0.6,col="black")

lims <- hist(fit$draws("b[1,1]"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,4, 4),freq=FALSE,ylim=c(0,ymax+1),breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for intra-county rate of infection in 1st cty",
     xlab = bquote(beta[11]))
hist(fit$draws("b[1,1]"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = beta[1,1],col="black")

lims <- hist(fit$draws("b[1,2]"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,4, 4),freq=FALSE,ylim=c(0,ymax+1),breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for inter-county rate of infection for NN",
     xlab = bquote(beta[ij]~", "~abs(i - j) == 1))
hist(fit$draws("b[1,2]"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = beta[1,2],col="black")

lims <- hist(fit$draws("b[1,3]"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,4, 4),freq=FALSE,ylim=c(0,ymax+1),breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for inter-county rate of infection for 2NN",
     xlab = bquote(beta[ij]~", "~abs(i - j) == 2))
hist(fit$draws("b[1,3]"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = beta[1,3],col="black")

lims <- hist(-log(1-fit$draws("gamma[1]")),plot=FALSE)
ymax <- lims$density |> max()
hist(-log(1-rbeta(1e6,0.8 * 15, 0.2 * 15)),freq=FALSE,ylim=c(0,ymax+1), border = NA,breaks=100,
     main = "Prior (grey) vs. posterior (red) for rate of recovery cty 1",
     xlab = bquote(-log(1-gamma[1])))
hist(-log(1-fit$draws("gamma[1]")),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = -log(1-gamma[1]),col="black")

lims <- hist(-log(1-fit$draws("gamma[2]")),plot=FALSE)
ymax <- lims$density |> max()
hist(-log(1-rbeta(1e6,0.8 * 15, 0.2 * 15)),freq=FALSE,ylim=c(0,ymax+1), border = NA,breaks=100,
     main = "Prior (grey) vs. posterior (red) for rate of recovery cty 2",
     xlab = bquote(-log(1-gamma[2])))
hist(-log(1-fit$draws("gamma[2]")),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = -log(1-gamma[2]),col="black")


plot(fit$summary("b")$mean, as.vector(beta))
abline(0,1)

plot(fit$summary("gamma")$mean, as.vector(gamma))
abline(0,1)

fit$summary("gamma")
sim_pars <- fit$sampler_diagnostics()
colSums(sim_pars[,,"treedepth__"] == 14)
colMeans(sim_pars[,,"treedepth__"])

colMeans(sim_pars[,,"stepsize__"])

meta_data <- fit$metadata()
names(meta_data)
meta_data$init

inits <- fit$init()
summary(inits[[1]]$b)
summary(inits[[2]]$b)
summary(inits[[3]]$b)
summary(inits[[4]]$b)

saveRDS(fit, file = "4-cty-spatial-struct.RDS")
fit$summary("i0")
fit$summary("p")
hist(fit$draws("p"))
abline(v = 0.6, col = "red")
fit$summary("b")
gamma
plot(fit$summary("gamma")$mean,gamma)
abline(0,1)
plot(fit$summary("b")$mean,beta)
abline(0,1)
plot(fit$draws("b[1]"),fit$draws("gamma[1]"))

hist(fit$draws("b[1]"))
abline(v = beta[1],col="red")

hist(fit$draws("b[2]"))
abline(v = beta[2],col="red")

hist(fit$draws("b[9]"))
abline(v = beta[9],col="red")

hist(fit$draws("gamma[9]"))
abline(v = gamma[9],col="red")

dim(fit$draws("gamma[9]"))
plot(fit$draws("gamma[9]")[,1,],fit$draws("b[9]")[,1,])
plot(fit$draws("gamma[9]")[,2,],fit$draws("b[9]")[,2,])
plot(fit$draws("gamma[9]")[,3,],fit$draws("b[9]")[,3,])
plot(fit$draws("gamma[9]")[,4,],fit$draws("b[9]")[,4,])


plot(fit$draws("i0[9]")[,1,],fit$draws("b[9]")[,1,])
plot(fit$draws("i0[9]")[,1,],fit$draws("b[9]")[,2,])
plot(fit$draws("i0[9]")[,1,],fit$draws("b[9]")[,3,])
plot(fit$draws("i0[9]")[,1,],fit$draws("b[9]")[,4,])

np_fit <- nuts_params(fit)
mcmc_pairs(fit$draws(c("p","i0","gamma","b_self")), np = np_fit, pars = c("p","gamma[3]","b_self[3]","gamma[4]","b_self[4]"),
           off_diag_args = list(size = 0.75))

z_t_d <- fit$draws("i_t", format = "draws_array") |> posterior::as_draws_rvars()
z_t_d <- z_t_d$i_t
qpt025 <- quantile(z_t_d,0.025)
qpt975 <- quantile(z_t_d,0.975)
idx <- 1:2
ylims <- range(c(
  range(sweep(qpt975[1,,],MARGIN = 2, STATS = pop_size, FUN = "*")[,idx]),
  range(sweep(qpt025[1,,],MARGIN = 2, STATS = pop_size, FUN = "*")[,idx])
))
matplot((apply(I,2,cumsum) - apply(R,2,cumsum))[,idx], xlab = "Time", 
     ylab = "Prevalence", main = "Latent prevalence v. time first 2 cties",
     pch = 19,ylim=ylims)
matlines(sweep(mean(z_t_d),MARGIN = 2, STATS = pop_size, FUN = "*")[,idx],lty=1)
matlines(sweep(qpt025[1,,],MARGIN = 2, STATS = pop_size, FUN = "*")[,idx],lty=2)
matlines(sweep(qpt975[1,,],MARGIN = 2, STATS = pop_size, FUN = "*")[,idx],lty=2)
