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

#counties=c(1,2,3,8)#c(1,2,3,4,5,6,7,8,10,11,12) # dist/10


N_C <- 11
TT <- 50
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

#sigma <- 2#.05 # Standard deviation of noise
#rho <- .9#1.2  # Spatial range parameter
#decay_rate_space <- .15#2 # Spatial decay rate
rho_se = 0
rho_ir = 0.1
rho_ei = 0.1
phi1=.9

first=sample(3:8,N_C,replace=T);first

#pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
pop=1e4*sample(1:N_C)#pop[counties,]

counties=c(1,2,3,4,5,6,7,8,10,11,12) # dist/10
pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
pop=pop$Population_2020[counties]#1e4*sample(5:N_C)#pop[counties,]

sigma=.05
log_beta <- numeric(TT-1)
log_beta[1] <- rnorm(1,0,sigma/sqrt(1-phi1^2));log_beta[1]

for (t in 2:(TT-1)) {
  log_beta[t] <- phi1 * log_beta[t - 1] + rnorm(1,0,sigma)
}
exp(log_beta) %>% plot
beta=exp(log_beta);beta


# create an empty list

gamma <- runif(N_C, min = 0.75, max = 0.81)
eta <- runif(N_C, min = 0.72, max = 0.78)





pop_size <- pop#pop$Population_2020[1:N_C]
#E_0 <- round(.009*pop_size)#ample(10:20,N_C,TRUE)#round(c(11902,6600,1120,737))#
I_0 <- round(.005*pop_size)#

S_0 <- pop_size - I_0
R <- S <- I <- E <- matrix(0,TT, N_C)
ii <- SE <- IR <- EI <- matrix(0,TT-1, N_C)
p_detect <- 0.4
for (i in 1:N_C){
S[1:(first[i]-1),i]= 0     
S[first[i],i] <- S_0[i]

I[1:(first[i]-1),i] = 0
I[first[i],i] <- I_0[i]

EI[first[i],i]=I_0[i]
}


#imp_rate=50


for (ct in  1:N_C) {
    for (t in 2:(TT)){
    if (t > first[ct]) {
     SE[t-1,ct]=rbetabinom(1,S[t-1,ct],1-exp(- sum(beta[t-1] * I[t-1,ct]/ pop_size[ct])),rho=rho_se) 
     if (EI[t-1,ct] == 0) {
          EI[t-1,ct]=rbetabinom(1,E[t-1,ct],eta[ct],rho=rho_ei) 
     } 
     IR[t-1,ct]=rbetabinom(1,I[t-1,ct],gamma[ct],rho=rho_ir)
     E[t,ct]=E[t-1,ct]+SE[t-1,ct]-(EI[t-1,ct]*as.numeric(t> (first[ct]+1))) #+ rpois(1,pop_size[ct]/1000)
     I[t,ct]=I[t-1,ct]+EI[t-1,ct]-IR[t-1,ct]
     S[t,ct]=S[t-1,ct]-SE[t-1,ct]
     R[t,ct]=R[t-1,ct]+IR[t-1,ct]
  }
  }
}


ii=matrix(rbinom(N_C*(TT-1),EI,p_detect),nrow=TT-1);ii
#ii_0=rbinom(N_C,I_0,p_detect)

#for (i in 1:N_C){
#ii[first[i],i]=ii_0[i]
#}

#quartz()
ii %>% as_tibble %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
ggplot(aes(y=Cases,x=date,group=County,color=County)) + geom_line() + facet_wrap(~County,scales="free_y") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


tt <- cmdstan_model("SEIR_betabin_ar1_beta_zeros.stan")

dat <- 
  list(
    ii = ii,
    TT = TT,
    N_C = N_C,
    pop_size = pop_size,
    first=first
  )

#saveRDS(dat,file="dat2.rds")
dat=readRDS("dat2.rds")
#seed(123)
fit <- tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.99,
                 max_treedepth = 14,
                 init = \() {list(u_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  v_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  w_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  p = rbeta(1, 4, 4),
                                  phi=runif(1,0,1),
                                  Z=runif(TT,-.25,.25),
                                  sigma = runif(1, 0, 2),
                                  i0 = rbeta(N_C, 0.01*50, 0.99*50),
                                  #rho_se = runif(1, 0.0001, 0.01),
                                  #rho_ei = runif(1, 0.05, .25),
                                  #rho_ir = runif(1, 0, 1),
                                  gamma = rbeta(N_C, 0.7 * 6, 0.3 * 6))},
                 iter_warmup = 1000,
                 iter_sampling = 1000, parallel_chains = 4)#,
                 #step_size = 1.5e-3)


np_fit <- nuts_params(fit)
quartz()
mcmc_pairs(fit$draws(c("p","gamma","beta","rho_ei","phi","sigma","eta","i0")), np = np_fit, pars = c("p","beta[10]","eta[1]","gamma[1]","i0[2]"),
            off_diag_args = list(size = 0.75))


betas=fit$draws("beta", format = "draws_array") |> posterior::as_draws_rvars()

1:49 %>% purrr::map_df(function(x){
qs=quantile(betas$beta[x],c(0.025,0.975))
data.frame(est=mean(betas$beta[x]),low=qs[1],high=qs[2],week=x)
}
) %>%
ggplot(aes(x=week,y=est)) + geom_point() + geom_line() + geom_ribbon(aes(ymin=low,ymax=high),alpha=.25) + geom_line(aes(y=beta),color="red")



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
ggplot(aes(x=date,y=Cases,color=County)) + geom_point() + facet_wrap(~County,scales="free_y") + 
geom_line(data=obs,aes(x=date,y=mean),color="black") + geom_line(data=obs,aes(x=date,y=lwr),linetype=2,color="black") + 
geom_line(data=obs,aes(x=date,y=upr),linetype=2,color="black") +
geom_point(data=truth,aes(x=date,y=Cases),color="#5757b8",alpha=.25) +
 theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
 
 ggsave(filename="ex_sim.png",width=8,height=6)

p=.9
phi=10
rbeta(1e6,p * phi, (1-p) * phi) %>% hist

