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

counties=c(1,2,3,8)#c(1,2,3,4,5,6,7,8,10,11,12) # dist/10


N_C <- 4
TT <- 20
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

sigma <- 2#.05 # Standard deviation of noise
rho <- .9#1.2  # Spatial range parameter
decay_rate_space <- .15#2 # Spatial decay rate
rho_se = .00097
rho_ir = .657
rho_ei = .0015


DiffsMat=expand.grid(1:N_C, 1:N_C) #%>% mutate(dist=abs(Var1-Var2))
DiffsMat %>% mutate(dist=abs(Var1-Var2)) -> D

pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
pop=pop[counties,]
# get the distance between each county (Admin2) in dat_final
dist=matrix(0,nrow=N_C,ncol=N_C)
for(i in 1:N_C){
  for(j in 1:N_C){
    dist[i,j]=geosphere::distHaversine(as.matrix(pop %>% dplyr::select(lon=Longitude,lat=Latitude))[c(i,j),])/1000
  }
} 

betamn <- 0.5 + (cos(20:(TT+18)/6) + 1) / 2.5;betamn
logbetamn=log(betamn)

# create an empty list
beta_list = list()

for (each in 1:(TT-1)){
D <- dist#matrix(D$dist,nrow=N_C)
Sigma <- sigma^2 * exp(-(D/10) / rho)
log_beta_diag=mvrnorm(1, rep(logbetamn[each], N_C), Sigma);log_beta_diag
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
beta_list[[each]] <- beta
}

gamma <- runif(N_C, min = 0.75, max = 0.81)
eta <- runif(N_C, min = 0.72, max = 0.78)





pop_size <- pop$Population_2020[1:N_C]#1e4 * sample(1:10,N_C,TRUE)
E_0 <- round(c(11902,6600,1120,737))#round(.9*pop_size/1000)#ample(10:20,N_C,TRUE)
I_0 <- round(c(1530,3172,2761,670))#round(c(6632.3402, 9638.4035, 3168.7845,  619.8996, 1376.2502,  263.1645, 1004.1161,
#2101.1389, 1019.0036,  372.1891,  230.9937))#round(.1*pop_size/1000)#ample(10:20,N_C,TRUE)

S_0 <- pop_size - E_0 - I_0
R <- S <- I <- E <- matrix(NA_real_,TT, N_C)
ii <- SE <- IR <- EI <- matrix(NA_real_,TT-1, N_C)
p_detect <- 0.4
S[1,] <- S_0
E[1,] <- E_0
I[1,] <- I_0
R[1,] <- 0 

#imp_rate=50

for (t in 1:(TT-1)){
    for (ct in  1:N_C) {
    SE[t,ct]=rbetabinom(1,S[t,ct],1-exp(- sum(beta_list[[t]][ct,] * I[t,]/ pop_size)),rho=rho_se) 
    EI[t,ct]=rbetabinom(1,E[t,ct],eta[ct],rho=rho_ei) 
    IR[t,ct]=rbetabinom(1,I[t,ct],gamma[ct],rho=rho_ir)
    E[t+1,ct]=E[t,ct]+SE[t,ct]-EI[t,ct] #+ rpois(1,pop_size[ct]/1000)
    I[t+1,ct]=I[t,ct]+EI[t,ct]-IR[t,ct]
    S[t+1,ct]=S[t,ct]-SE[t,ct]
    R[t+1,ct]=R[t,ct]+IR[t,ct]
  }
  }


ii=matrix(rbinom(N_C*(TT-1),EI,p_detect),nrow=TT-1);ii

quartz()
ii %>% as_tibble %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
ggplot(aes(y=Cases,x=date,group=County,color=County)) + geom_line() + facet_wrap(~County,scales="free_y") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


tt <- cmdstan_model("SEIR_betabin_vary_beta.stan")

dat <- 
  list(
    ii = ii,
    TT = TT,
    N_C = N_C,
    pop_size = pop_size,
    D=D/10,
    b_freq=3
  )


seed(123)
fit <- tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.99,
                 max_treedepth = 15,
                 init = \() {list(u_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  v_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  w_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  p = rbeta(1, 4, 4),
                                  log_beta_diag = lapply(1:(((TT-1) %/% dat$b_freq)+1), function(t) rnorm(N_C, 0, 0.001)),
                                  i0 = rbeta(N_C, 0.01*50, 0.99*50),
                                  e0 = rbeta(N_C, 0.01*50, 0.99*50),
                                  rho_se = runif(1, 0.0001, 0.001),
                                  rho_ei = runif(1, 0.05, .25),
                                  rho_ir = runif(1, 0, 1),
                                  gamma = rbeta(N_C, 0.7 * 6, 0.3 * 6))},
                 iter_warmup = 100,
                 iter_sampling = 100, parallel_chains = 4,
                 step_size = 1.5e-3)


np_fit <- nuts_params(fit)

fit$summary("beta")

quartz()
mcmc_pairs(fit$draws(c("p","decay_rate_space","gamma","beta","rho_se")), np = np_fit, pars = c("p","beta[1,1,1]","decay_rate_space","rho_se"),
           off_diag_args = list(size = 0.75))

quartz()
mcmc_pairs(fit$draws(c("rho_ei","rho_ir","rho_se")), np = np_fit, pars = c("rho_ei","rho_ir","rho_se"),
           off_diag_args = list(size = 0.75))

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



z_t_d <- fit$draws("ei_t", format = "draws_array") |> posterior::as_draws_rvars()
z_t_d <- z_t_d$ei_t
qpt025 <- quantile(z_t_d,0.025)
qpt975 <- quantile(z_t_d,0.975)
pop4=pop_size
obs=rbind(data.frame(value="mean",sweep(mean(z_t_d)[-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
data.frame(value="lwr",sweep(qpt025[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
data.frame(value="upr",sweep(qpt975[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")))

obs=obs %>% gather(County,Cases,-value) %>% pivot_wider(names_from=value,values_from=Cases) %>% unnest() %>% mutate(date=rep(1:(TT-1),N_C))

truth=EI %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
mutate(County= gsub("V", "X", County))

quartz()
dat$ii %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
mutate(County= gsub("V", "X", County)) %>%
ggplot(aes(x=date,y=Cases,color=County)) + geom_point() + facet_wrap(~County) + 
geom_line(data=obs,aes(x=date,y=mean),color="black") + geom_line(data=obs,aes(x=date,y=lwr),linetype=2,color="black") + 
geom_line(data=obs,aes(x=date,y=upr),linetype=2,color="black") +
geom_point(data=truth,aes(x=date,y=Cases),color="#5757b8",alpha=.25) +
 theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
 
 ggsave(filename="ex_sim.png",width=8,height=6)



