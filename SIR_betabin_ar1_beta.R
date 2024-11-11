library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(MASS)
library(VGAM)

set.seed(32)
#set.seed(123)
#beta <- 1.05
N_C <- 11
TT <- 18

rho_si = 0
rho_ir = .852

phi1=.916

counties=c(1,2,3,4,5,6,7,8,10,11,12) # dist/10
pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))
pop=pop$Population_2020[counties]#1e4*sample(5:N_C)#pop[counties,]

sigma=.15
log_beta <- numeric(TT-1)
log_beta[1] <- rnorm(1,0,sigma/sqrt(1-phi1^2));log_beta[1]

for (t in 2:(TT-1)) {
  log_beta[t] <- phi1 * log_beta[t - 1] + rnorm(1,0,sigma)
}
exp(log_beta) %>% plot
beta=exp(log_beta);beta

gamma <- c(.65,.63,.7,.7,.73,.96,.8,.64,.76,.65,.89)#runif(N_C, min = 0.63, max = 0.96)
pop_size <- pop#1e4 * sample(1:10,N_C,TRUE)
I_0 <- c(133000,73650,8000,6000,12500,3500,3200,8700,7000,3200,1900)#round(pop_size/500)#ample(10:20,N_C,TRUE)
S_0 <- pop_size - I_0
R <- S <- I <- matrix(NA_real_,TT, N_C)
ii <- SI <- IR <- matrix(NA_real_,TT-1, N_C)
p_detect <- 0.77
S[1,] <- S_0
I[1,] <- I_0
R[1,] <- 0 

for (t in 1:(TT-1)){
    for (ct in  1:N_C) {
    SI[t,ct]=rbetabinom(1,S[t,ct],1-exp(- sum(beta[t] * I[t,ct]/ pop_size[ct])),rho=rho_si)
    IR[t,ct]=rbetabinom(1,I[t,ct],gamma[ct],rho=rho_ir)
    I[t+1,ct]=I[t,ct]+SI[t,ct]-IR[t,ct]
    S[t+1,ct]=S[t,ct]-SI[t,ct]
    R[t+1,ct]=R[t,ct]+IR[t,ct]
  }
  }


ii=matrix(rbinom(N_C*(TT-1),SI,p_detect),nrow=TT-1);ii

matplot(ii, type="l")

tt <- cmdstan_model("SIR_betabin_ar1_beta.stan")

dat <- 
  list(
    ii = ii,
    TT = TT,
    N_C = N_C,
    pop_size = pop_size
  )

saveRDS(dat, "dat.rds")

seed(123)
fit <- tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.95,
                 max_treedepth = 14,
                 init = \() {list(u_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  w_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  p = rbeta(1, 4, 4),
                                  phi=runif(1,.25,.5),
                                  Z=runif(TT,-.25,.25),
                                  sigma = runif(1, 0, 1),
                                  i0 = rbeta(N_C, 0.01*50, 0.99*50),
                                  #rho_si = runif(1, 0.0001, 0.005),
                                  gamma = rbeta(N_C, 0.7 * 6, 0.3 * 6))},
                 iter_warmup = 1000,
                 iter_sampling = 1000, parallel_chains = 4)#,
                 #step_size = .005)

fit$init()
fit$sampler_diagnostics()

lims <- hist(fit$draws("p"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.5 * 4, 0.5 * 4),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("p"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = p_detect,col="black")

fit$draws("p") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()
fit$summary("p")
fit$draws("beta[2]") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()

np_fit <- nuts_params(fit)
mcmc_pairs(fit$draws(c("p","gamma","beta","rho_ir")), np = np_fit, pars = c("p","beta[1]","beta[2]","beta[3]","rho_ir"),
           off_diag_args = list(size = 0.75))

betas=fit$draws("beta", format = "draws_array") |> posterior::as_draws_rvars()


1:30 %>% purrr::map_df(function(x){
qs=quantile(betas$beta[x],c(0.025,0.975))
data.frame(est=mean(betas$beta[x]),low=qs[1],high=qs[2],week=x)
}
) %>% mutate(truth=c(beta,1)) %>%
ggplot(aes(x=week,y=est)) + geom_point() + geom_line() + geom_ribbon(aes(ymin=low,ymax=high),alpha=.25) + geom_line(aes(y=truth),color="red")


fit$draws("si_t[2,1]") %>% as_tibble %>% gather() %>% ggplot(aes(y=value*pop_size[1],x=rep(1:1000,4),group=key,color=key)) + geom_line()



sim_pars <- fit$sampler_diagnostics()
colSums(sim_pars[,,"treedepth__"] == 14)
colMeans(sim_pars[,,"treedepth__"])

colMeans(sim_pars[,,"stepsize__"])

meta_data <- fit$metadata()
names(meta_data)
meta_data$data
meta_data$init


z_t_d <- fit$draws("si_t", format = "draws_array") |> posterior::as_draws_rvars()
z_t_d <- z_t_d$si_t
qpt025 <- quantile(z_t_d,0.025)
qpt975 <- quantile(z_t_d,0.975)

pop4=pop_size


obs=rbind(data.frame(value="mean",sweep(mean(z_t_d)[-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
data.frame(value="lwr",sweep(qpt025[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
data.frame(value="upr",sweep(qpt975[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")))

obs=obs %>% gather(County,Cases,-value) %>% pivot_wider(names_from=value,values_from=Cases) %>% unnest() %>% mutate(date=rep(1:(TT-1),N_C)) #%>%
#filter(County %in% unique(truth$County)[1:4])#



truth=SI %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
mutate(County= gsub("V", "X", County)) #%>% filter(County %in% unique(truth$County)[1:4])


obs %>% left_join(truth,by=c("County","date")) %>% 
ggplot(aes(x=date,y=Cases,color=County)) + geom_point() + 
geom_line(aes(y=mean),linetype="dashed") + geom_line(aes(y=lwr),linetype="dashed") + geom_line(aes(y=upr),linetype="dashed") +
facet_wrap(~County) 


dat$ii %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
mutate(County= gsub("V", "X", County)) %>% #filter(County %in% unique(truth$County)[1:4]) %>%
ggplot(aes(x=date,y=Cases,color=County)) + geom_point() + facet_wrap(~County) + 
geom_line(data=obs,aes(x=date,y=mean),color="black") + geom_line(data=obs,aes(x=date,y=lwr),linetype=2,color="black") + 
geom_line(data=obs,aes(x=date,y=upr),linetype=2,color="black") +
geom_point(data=truth,aes(x=date,y=Cases),color="#5757b8",alpha=.25) +
 theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
 
 ggsave(filename="ex_sim.png",width=8,height=6)



# Given mean (mu) and precision (phi)
mu <- 0.25    # Example mean
phi <- 5      # Example precision

# Convert to alpha and beta
alpha <- mu * phi
beta <- (1 - mu) * phi
simulated_data <- rbeta(n = 1000, shape1 = alpha, shape2 = beta)

# Visualize the simulated data
hist(simulated_data, breaks = 30, main = "Simulated Beta Distribution", xlab = "Value")
