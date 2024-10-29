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
case=case %>% ungroup() %>% mutate(date=(date = 7 * (as.numeric(date - min(date)) %/% 7) + min(date))) %>% 
  group_by(date,Admin2) %>% summarize(cases=sum(cases)) %>%#,tests=sum(tests_combined_total-lag(tests_combined_total),na.rm=T))#,
                               #vax=sum(vax-lag(vax),na.rm=T))
                               filter(date>=ymd("2020-04-15")) %>% filter(date<ymd("2021-06-02"))
                               #filter(date>=ymd("2020-10-06")) %>% filter(date<ymd("2021-02-02"))


case=case %>% group_by(Admin2) %>% mutate(new_cases=cases-lag(cases)) %>% arrange(Admin2,date) %>% filter(!is.na(new_cases))# %>% filter(Admin2 %in% readRDS("final_counties.rds")) %>%# filter(any(cases>1000)) %>%
#saveRDS(unique(dat_final$Admin2),"final_counties.rds")


dat_final=case %>% filter(Admin2 %in% readRDS("final_counties.rds"))#group_by(Admin2) %>% filter(any(cases>100))
pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))

# get the distance between each county (Admin2) in dat_final
dist=matrix(0,nrow=length(unique(dat_final$Admin2)),ncol=length(unique(dat_final$Admin2)))
for(i in 1:12){
  for(j in 1:12){
    dist[i,j]=geosphere::distHaversine(as.matrix(pop %>% dplyr::select(lon=Longitude,lat=Latitude))[c(i,j),])/1000
  }
} 

d1=dat_final %>% filter(date<ymd("2020/8/19")) %>% 
ungroup() %>% rename(County="Admin2") %>% left_join(pop,by="County") %>% arrange(desc(Population_2020),date) %>%
dplyr::select(-Population_2020,-Latitude,-Longitude,-cases) %>% pivot_wider(names_from=County,values_from=new_cases) %>%
dplyr::select(-date)
#d1[3,12]=0

tt <- cmdstan_model("SEIR_betabin_ar1_beta.stan")#cmdstan_model("SEIR_betabin_vary_beta_nospat.stan")#"stoch_beta_spatial_SI_utah_betabin.stan")

d1[6,9]=1
#d1=d1 %>% dplyr::select(`Salt Lake`,Utah,Davis,`Weber-Morgan`)

    TT = nrow(d1)+1

#1,2,3,8 for 4 counties # dist/10
counties=c(1,2,3,4,5,6,7,8,10,11,12) # dist/10
#counties=c(1,2,3,8)
N_C = length(counties)#ncol(d1)
dat <- 
  list(
    ii = as.matrix(d1)[,counties],
    TT = TT,
    N_C = N_C,
    pop_size = pop$Population_2020[counties],
    D=dist[counties,counties]/10,
    b_freq=2
  )

? cmdstanr::cmdstan_model

set.seed(123)
fit <- tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.95,
                 max_treedepth = 15,
                 init = \() {list(u_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  v_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  w_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  p = rbeta(1, 4, 4),
                                  phi=runif(1,0,1),
                                  Z=runif(TT,-.25,.25),
                                  sigma = runif(1, 0, 1),
                                  i0 = rbeta(N_C, 0.01*50, 0.99*50),
                                  #rho_si = runif(1, 0.0001, 0.005),
                                  gamma = rbeta(N_C, 0.7 * 6, 0.3 * 6))},
                 iter_warmup = 1000,
                 iter_sampling = 1000, parallel_chains = 4)#,
                 #step_size = .005)

saveRDS(fit,"./output_4/fit_4chn_SEIR.rds")

fit=readRDS("./output_11/fit_10chn.rds")
#fit=readRDS("./output_4/fit.rds")

fit$summary("sigma") 
fit$summary("p")
fit$diagnostic_summary()
fit$sampler_diagnostics()

#diags <- fit$sampler_diagnostics()
#diags[1,,"stepsize__"]

tt <- fit$inv_metric(matrix = FALSE)
summary(tt$`3`)

lims <- hist(fit$draws("p"),plot=FALSE)
ymax <- lims$density |> max()
hist(rbeta(1e6,0.5 * 4, 0.5 * 4),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("p"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)
abline(v = p_detect,col="black")

lims <- hist(fit$draws("phi"),plot=FALSE)
ymax <- lims$density |> max()
hist(runif(1e6,-1, 1),freq=FALSE,breaks=100, border = NA,
     main = "Prior (grey) vs. posterior (red) for prob of detection",
     xlab = "p", ylim = c(0,ymax+1))
hist(fit$draws("phi"),freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,border = NA)

fit$draws("p") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()
fit$summary("p")
fit$draws("beta[2]") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()

np_fit <- nuts_params(fit)
quartz()
mcmc_pairs(fit$draws(c("p","gamma","beta","rho_ei","phi","sigma","eta")), np = np_fit, pars = c("p","beta[2]","beta[3]","eta[1]","phi","rho_ei"),
           off_diag_args = list(size = 0.75))

betas=fit$draws("beta", format = "draws_array") |> posterior::as_draws_rvars()


1:17 %>% purrr::map_df(function(x){
qs=quantile(betas$beta[x],c(0.025,0.975))
data.frame(est=mean(betas$beta[x]),low=qs[1],high=qs[2],week=x)
}
) %>%
ggplot(aes(x=week,y=est)) + geom_point() + geom_line() + geom_ribbon(aes(ymin=low,ymax=high),alpha=.25) #+ geom_line(aes(y=truth),color="red")





sim_pars <- fit$sampler_diagnostics()
colSums(sim_pars[,,"treedepth__"] == 14)
colMeans(sim_pars[,,"treedepth__"])

sim_pars[,,"treedepth__"] %>% as_tibble

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



#truth=SI %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
#mutate(County= gsub("V", "X", County)) #%>% filter(County %in% unique(truth$County)[1:4])

fit$summary("p")

truth=dat$ii %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) 
truth$County

unique(truth$County)
quartz()
obs %>% mutate(County=rep(unique(truth$County),each=17)) %>% left_join(truth,by=c("County","date")) %>% 
ggplot(aes(x=date,y=mean,color=County)) + geom_point() + 
geom_line(aes(y=lwr),linetype="dashed") + geom_line(aes(y=upr),linetype="dashed") + 
#geom_point(aes(y=Cases),color="red") +
facet_wrap(~County) 


dat$ii %>% as_tibble() %>% mutate(date=1:(TT-1)) %>% gather(County,Cases,-date) %>% 
mutate(County= gsub("V", "X", County)) %>% #filter(County %in% unique(truth$County)[1:4]) %>%
ggplot(aes(x=date,y=Cases,color=County)) + geom_point() + facet_wrap(~County) + 
geom_line(data=obs,aes(x=date,y=mean),color="black") + geom_line(data=obs,aes(x=date,y=lwr),linetype=2,color="black") + 
geom_line(data=obs,aes(x=date,y=upr),linetype=2,color="black") +
geom_point(data=truth,aes(x=date,y=Cases),color="#5757b8",alpha=.25) +
 theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

fit$summary() %>% as_tibble() %>% filter(!grepl("u_t",variable),!grepl("w_t",variable),!grepl("i0",variable),!grepl("gamma",variable),
!grepl("Z",variable),,!grepl("s_t",variable),!grepl("si_t",variable),!grepl("ir_t",variable),!grepl("i_t",variable))

fit$summary() %>% as_tibble() %>% filter(grepl("Z",variable))

