library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(MASS)
library(VGAM)

set.seed(32)
#set.seed(123)
#beta <- 1.05
#TT <- 100
N_C <- 10

#  p=plogis(peta0 + peta1 * tests/1000)
#  plogis(0)

#Simulation parameters:
# p = .25, .4, .7
# sigma = .15
# rho = 1 
# decay_rate_space = 1.1
# gamma[N_C] = .75 - .9
# beta

# Spatial parameters
# Parameters

#fit$summary("rho_ir")


fit=readRDS(paste0("./output_",fitsize,"/fit_10chn.rds"))
fit$summary("rho")

sigma <- .27#.05 # Standard deviation of noise
rho <- 1#1.2  # Spatial range parameter
decay_rate_space <- 1.1#2 # Spatial decay rate
rho_si = .0003
rho_ir = .35


DiffsMat=expand.grid(1:N_C, 1:N_C) #%>% mutate(dist=abs(Var1-Var2))
DiffsMat %>% mutate(dist=abs(Var1-Var2)) -> D

D <- dist#matrix(D$dist,nrow=N_C)
Sigma <- sigma^2 * exp(-(D/10) / rho)
log_beta_diag=mvrnorm(1, rep(0, N_C), Sigma);log_beta_diag
log_beta = matrix(0, N_C, N_C)
diag(log_beta)=log_beta_diag

for (i in 1:(N_C-1)) {
  for (j in (i+1):(N_C)) {
    distance_from_center <- D[i,j]/10#abs(i - j)
    log_beta[i, j] <- log_beta[i, i] - decay_rate_space * distance_from_center
    log_beta[j, i] <- log_beta[i, j]
  }
}
beta=exp(log_beta)

gamma <- runif(N_C, min = 0.75, max = 0.85)
TT <- 18
pop_size <- pop$Population_2020[1:10]#1e4 * sample(1:10,N_C,TRUE)
I_0 <- round(pop_size/1000)#ample(10:20,N_C,TRUE)
S_0 <- pop_size - I_0
R <- S <- I <- matrix(NA_real_,TT, N_C)
ii <- SI <- IR <- matrix(NA_real_,TT-1, N_C)
p_detect <- 0.4
S[1,] <- S_0
I[1,] <- I_0
R[1,] <- 0 

for (t in 1:(TT-1)){
    for (ct in  1:N_C) {
    SI[t,ct]=rbetabinom(1,S[t,ct],1-exp(- sum(beta[ct,] * I[t,]/ pop_size)),rho=rho_si)
    IR[t,ct]=rbetabinom(1,I[t,ct],gamma[ct],rho=rho_ir)
    I[t+1,ct]=I[t,ct]+SI[t,ct]-IR[t,ct]
    S[t+1,ct]=S[t,ct]-SI[t,ct]
    R[t+1,ct]=R[t,ct]+IR[t,ct]
  }
  }


ii=matrix(rbinom(N_C*(TT-1),SI,p_detect),nrow=TT-1);ii

matplot(ii, type="l")

tt <- cmdstan_model("stoch_beta_spatial_SI_betabin.stan")

dat <- 
  list(
    ii = ii,
    TT = TT,
    N_C = N_C,
    pop_size = pop_size,
    D=D/10
  )

saveRDS(dat,"tst_dat.rds")

seed(123)
fit <- tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.95,
                 max_treedepth = 14,
                 init = \() {list(u_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  w_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  p = rbeta(1, 4, 4),
                                  log_beta_diag = rnorm(N_C,0,.001),
                                  #b_od = rgamma(N_C-2, 4, 4),
                                  #b_self = rgamma(N_C, 4, 4),
                                  i0 = rbeta(N_C, 0.01*50, 0.99*50),
                                  rho_si = runif(1, 0.0001, 0.005),
                                  gamma = rbeta(N_C, 0.7 * 6, 0.3 * 6))},
                 iter_warmup = 1000,
                 iter_sampling = 1000, parallel_chains = 4,
                 step_size = 1.5e-3)


lims <- hist(fit$draws("p"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.5 * 4, 0.5 * 4),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("p"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = p_detect,col="black")

fit$draws("p") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()
fit$summary("p")

lims <- hist(fit$draws("beta[1,1]"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,4, 4),freq=FALSE,ylim=c(0,ymax+1),breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for intra-county rate of infection in 1st cty",
     xlab = bquote(beta[1,1]))
hist(fit$draws("beta[1,1]"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = beta[1,1],col="black")



lims <- hist(fit$draws("b[1,2]"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,4, 4),freq=FALSE,ylim=c(0,ymax+1),breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for inter-county rate of infection for NN",
     xlab = bquote(beta[ij]~", "~abs(i - j) == 1))
hist(fit$draws("beta[1,2]"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = beta[1,2],col="black")

lims <- hist(fit$draws("b[1,3]"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,20, 4),freq=FALSE,ylim=c(0,ymax+1),breaks=100, border = NA,
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

lims <- hist(fit$draws("decay_rate_space"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,10, 2),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("decay_rate_space"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = decay_rate_space,col="black")

lims <- hist(fit$draws("sigma"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,2, 2),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("sigma"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = sigma,col="black")

lims <- hist(fit$draws("rho"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,2, 2),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("rho"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = rho,col="black")

lims <- hist(fit$draws("rho_ir"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,2, 2),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("rho_ir"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = rho_si,col="black")

lims <- hist(fit$draws("rho_si"),plot=FALSE)
ymax <- lims$density |> max()
hist(rgamma(1e6,2, 2),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("rho_si"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = rho_ir,col="black")

fit$draws("rho_si") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()


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
mcmc_pairs(fit$draws(c("p","decay_rate_space","gamma","beta","rho_si")), np = np_fit, pars = c("p","beta[1,1]","decay_rate_space","beta[1,2]","rho_si"),
           off_diag_args = list(size = 0.75))



z_t_d <- fit$draws("si_t", format = "draws_array") |> posterior::as_draws_rvars()
z_t_d <- z_t_d$si_t
qpt025 <- quantile(z_t_d,0.025)
qpt975 <- quantile(z_t_d,0.975)

pop4=pop_size


obs=rbind(data.frame(value="mean",sweep(mean(z_t_d)[-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
data.frame(value="lwr",sweep(qpt025[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
data.frame(value="upr",sweep(qpt975[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")))

obs=obs %>% gather(County,Cases,-value) %>% pivot_wider(names_from=value,values_from=Cases) %>% unnest() %>% mutate(date=rep(1:(TT-1),N_C))



truth=SI %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
mutate(County= gsub("V", "X", County))

dat$ii %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
mutate(County= gsub("V", "X", County)) %>%
ggplot(aes(x=date,y=Cases,color=County)) + geom_point() + facet_wrap(~County) + 
geom_line(data=obs,aes(x=date,y=mean),color="black") + geom_line(data=obs,aes(x=date,y=lwr),linetype=2,color="black") + 
geom_line(data=obs,aes(x=date,y=upr),linetype=2,color="black") +
geom_point(data=truth,aes(x=date,y=Cases),color="#5757b8",alpha=.25) +
 theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
 
 ggsave(filename="ex_sim.png",width=8,height=6)



