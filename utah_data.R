library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)
case = read_csv("./data/time_series_covid19_confirmed_US.csv") %>% filter(`Province_State`=="Utah")
test = read_csv("./data/time_series_covid19_US_testing_by_state.csv") %>% filter(state=="UT")
vax = read_csv("./data/COVID-19_Vaccinations_in_the_United_States_County.csv") %>% filter(Recip_State=="UT")
test=test %>% mutate(date=mdy(date)) %>% arrange(date) %>% dplyr::select(date,tests_combined_total)
vax=vax %>% mutate(date=mdy(Date)) %>% arrange(date) %>% dplyr::select(date,vax=Series_Complete_Yes) %>% 
  group_by(date) %>% summarize(vax=sum(vax))

case=case %>% dplyr::select(-c(1:5,7:11)) %>% gather("date","cases",-Admin2) %>% mutate(date=mdy(date)) %>% group_by(Admin2,date) %>% summarize(cases=sum(as.numeric(cases)))

#case=test %>% right_join(case,by="date") %>% left_join(vax,by="date") %>% mutate(vax=replace_na(vax,0))

ggplot(case,aes(x=date,y=cases)) + geom_line() + facet_wrap(~Admin2) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

case=case %>% ungroup() %>% mutate(date=(date = 7 * (as.numeric(date - min(date)) %/% 7) + min(date))) %>% 
  group_by(date,Admin2) %>% summarize(cases=sum(cases)) %>%#,tests=sum(tests_combined_total-lag(tests_combined_total),na.rm=T))#,
                               #vax=sum(vax-lag(vax),na.rm=T))
                               filter(date>=ymd("2020-04-15")) %>% filter(date<ymd("2021-06-02"))
                               #filter(date>=ymd("2020-10-06")) %>% filter(date<ymd("2021-02-02"))


data.frame(case %>% group_by(date) %>% summarize(sum(cases>0)))


case %>% filter(Admin2 %in% readRDS("final_counties.rds")) %>%# filter(any(cases>1000)) %>%
ggplot(aes(x=date,y=cases)) + facet_wrap(~Admin2) + geom_line() + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

data.frame(case %>% group_by(Admin2) %>% filter(any(cases>100)) %>% 
ungroup() %>% group_by(date) %>% summarize(sum(cases>0)))
###
case=case %>% group_by(Admin2) %>% mutate(new_cases=cases-lag(cases)) %>% arrange(Admin2,date) %>% filter(!is.na(new_cases))# %>% filter(Admin2 %in% readRDS("final_counties.rds")) %>%# filter(any(cases>1000)) %>%
#saveRDS(unique(dat_final$Admin2),"final_counties.rds")

case %>% filter(Admin2 %in% readRDS("final_counties.rds")) %>%# filter(any(cases>1000)) %>%
ggplot(aes(x=date,y=new_cases)) + facet_wrap(~Admin2) + geom_line() + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


dat_final=case %>% filter(Admin2 %in% readRDS("final_counties.rds"))#group_by(Admin2) %>% filter(any(cases>100))

pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))

# get the distance between each county (Admin2) in dat_final
dist=matrix(0,nrow=length(unique(dat_final$Admin2)),ncol=length(unique(dat_final$Admin2)))
for(i in 1:12){
  for(j in 1:12){
    dist[i,j]=geosphere::distHaversine(as.matrix(pop %>% dplyr::select(lon=Longitude,lat=Latitude))[c(i,j),])/1000
  }
} 

dat_final  %>% 
filter(date<ymd("2020/8/20")) %>%#, 
#Admin2 %in% c("Salt Lake","Utah","Davis","Weber-Morgan")) %>% 
ggplot(aes(x=date,y=new_cases)) + geom_line() + facet_wrap(~Admin2) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

d1=dat_final %>% filter(date<ymd("2020/8/19")) %>% 
ungroup() %>% rename(County="Admin2") %>% left_join(pop,by="County") %>% arrange(desc(Population_2020),date) %>%
dplyr::select(-Population_2020,-Latitude,-Longitude,-cases) %>% pivot_wider(names_from=County,values_from=new_cases) %>%
dplyr::select(-date)
#d1[3,12]=0

tt <- cmdstan_model("SEIR_betabin.stan")#"stoch_beta_spatial_SI_utah_betabin.stan")

d1[6,9]=1
#d1=d1 %>% dplyr::select(`Salt Lake`,Utah,Davis,`Weber-Morgan`)

    TT = nrow(d1)+1

#1,2,3,8 for 4 counties # dist/10
counties=c(1,2,3,4,5,6,7,8,10,11,12) # dist/10
counties=c(1,2,3,8)
N_C = length(counties)#ncol(d1)
dat <- 
  list(
    ii = as.matrix(d1)[,counties],
    TT = TT,
    N_C = N_C,
    pop_size = pop$Population_2020[counties],
    D=dist[counties,counties]/10
  )

? cmdstanr::cmdstan_model

set.seed(123)
fit <- tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.99,
                 max_treedepth = 15,
                 init = \() {list(u_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  v_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  w_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  p = rbeta(1, 4, 4),
                                  log_beta_diag = rnorm(N_C,0,.001),
                                  i0 = rbeta(N_C, 0.01*50, 0.99*50),
                                  e0 = rbeta(N_C, 0.01*50, 0.99*50),
                                  rho_se = runif(1, 0.0001, 0.005),
                                  rho_ei = runif(1, 0.5, 1),
                                  rho_ir = runif(1, 0.5, 1),
                                  gamma = rbeta(N_C, 0.7 * 6, 0.3 * 6))},
                 iter_warmup = 1500,
                 iter_sampling = 1500, parallel_chains = 4,
                 step_size = 1.5e-3)
                 #output_dir = "./output_4")
saveRDS(fit,"./output_4/fit_10chn.rds")

fit=readRDS("./output_11/fit_10chn.rds")
#fit=readRDS("./output_4/fit.rds")

fit$summary("sigma") 
fit$summary("p")
fit$diagnostic_summary()


lims <- hist(fit$draws("p"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.5 * 4, 0.5 * 4),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("p"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)

png("four_county_p.png")
ggplot() + 
geom_histogram(data=data.frame(p=rbeta(1e6,0.5 * 4, 0.5 * 4)),
aes(x=p,..density..),alpha=.25,fill="grey",linetype=2) +
geom_histogram(data=fit$draws("p") %>% as_tibble %>% gather() %>% mutate(key=factor(key)), 
aes(x=value,y=..density..,fill=1),alpha=0.75) + 
scale_fill_viridis_c() +
theme_minimal() +
theme(legend.position="none")
dev.off()


fit$draws("p") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1500,10),group=key,color=key)) + geom_line()
fit$draws("beta[3,3]") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()

fit$draws("beta[1,1]") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()
fit$draws("beta[1,2]") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()


fit$draws("rho_ir") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1500,4),group=key,color=key)) + geom_line()
fit$draws("rho") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()
fit$draws("decay_rate_space") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()


fit$draws("si_t") %>% as_tibble()

z_t_d <- fit$draws("si_t", format = "draws_array") |> posterior::as_draws_rvars()
z_t_d <- z_t_d$si_t
qpt025 <- quantile(z_t_d,0.025)
qpt975 <- quantile(z_t_d,0.975)
idx <- 1:4
ylims <- range(c(
  range(sweep(qpt975[1,,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")[,idx]),
  range(sweep(qpt025[1,,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")[,idx])
))
matplot(dat$ii[,idx], xlab = "Time", 
     ylab = "Prevalence", main = "Latent prevalence v. time first 2 cties",
     pch = 19,ylim=ylims) 
matlines(sweep(mean(z_t_d)[-1,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")[,idx],lty=1)
matlines(sweep(qpt025[1,-1,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")[,idx],lty=2)
matlines(sweep(qpt975[1,-1,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")[,idx],lty=2)

obs=rbind(data.frame(value="mean",sweep(mean(z_t_d)[-1,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")[,idx]),
data.frame(value="lwr",sweep(qpt025[1,-1,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")[,idx]),
data.frame(value="upr",sweep(qpt975[1,-1,],MARGIN = 2, STATS = dat$pop_size, FUN = "*")[,idx])) %>%
rename(`Salt Lake`=X1,Utah=X2,Davis=X3,Summit=X4) %>% gather(County,Cases,-value) %>% 
pivot_wider(names_from=value,values_from=Cases) %>% unnest() %>% mutate(date=rep(dts,4))

dts=unique(dat_final %>% filter(date<ymd("2020/8/19")) %>% pull(date))
png("inf_est_4.png")
dat$ii %>% as_tibble() %>% mutate(date=dts) %>%  gather(County,Cases,-date) %>% 
ggplot(aes(x=date,y=Cases,color=County)) + geom_point() + facet_wrap(~County) + 
geom_line(data=obs,aes(x=date,y=mean),color="black") + geom_line(data=obs,aes(x=date,y=lwr),linetype=2,color="black") + 
geom_line(data=obs,aes(x=date,y=upr),linetype=2,color="black") +
 theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


np_fit <- nuts_params(fit)
mcmc_pairs(fit$draws(c("p","decay_rate_space","gamma","beta","rho_se")), np = np_fit, pars = c("p","beta[1,1]","decay_rate_space","beta[1,2]","rho_se"),
           off_diag_args = list(size = 0.75))


hist(rgamma(1e6,5, 1),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,.1))


mu=.5
kappa=4

alp=mu*kappa
bet=(1-mu)*kappa

hist(rbeta(1e6,alp,bet),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))

