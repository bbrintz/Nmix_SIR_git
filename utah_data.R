library(tidyverse)
library(cmdstanr)
case = read_csv("./data/time_series_covid19_confirmed_US.csv") %>% filter(`Province_State`=="Utah")
test = read_csv("./data/time_series_covid19_US_testing_by_state.csv") %>% filter(state=="UT")
vax = read_csv("./data/COVID-19_Vaccinations_in_the_United_States_County.csv") %>% filter(Recip_State=="UT")
test=test %>% mutate(date=mdy(date)) %>% arrange(date) %>% select(date,tests_combined_total)
vax=vax %>% mutate(date=mdy(Date)) %>% arrange(date) %>% select(date,vax=Series_Complete_Yes) %>% 
  group_by(date) %>% summarize(vax=sum(vax))

case=case %>% select(-c(1:5,7:11)) %>% gather("date","cases",-Admin2) %>% mutate(date=mdy(date)) %>% group_by(Admin2,date) %>% summarize(cases=sum(as.numeric(cases)))

#case=test %>% right_join(case,by="date") %>% left_join(vax,by="date") %>% mutate(vax=replace_na(vax,0))

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

#saveRDS(unique(dat_final$Admin2),"final_counties.rds")
dat_final=case %>% filter(Admin2 %in% readRDS("final_counties.rds"))#group_by(Admin2) %>% filter(any(cases>100))

pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))

# get the distance between each county (Admin2) in dat_final
dist=matrix(0,nrow=length(unique(dat_final$Admin2)),ncol=length(unique(dat_final$Admin2)))
for(i in 1:12){
  for(j in 1:12){
    dist[i,j]=geosphere::distHaversine(as.matrix(pop %>% select(lon=Longitude,lat=Latitude))[c(i,j),])/1000
  }
} 


d1=dat_final %>% ungroup() %>% rename(County="Admin2") %>% left_join(pop,by="County") %>% arrange(desc(Population_2020),date) %>%
select(-Population_2020,-Latitude,-Longitude) %>% pivot_wider(names_from=County,values_from=cases) %>%
select(-date)
#d1[3,12]=0

tt <- cmdstan_model("stoch_beta_spatial_cum.stan")



    TT = nrow(d1)
    N_C = ncol(d1)

dat <- 
  list(
    ii = as.matrix(d1),
    TT = nrow(d1),
    N_C = ncol(d1),
    pop_size = pop$Population_2020,
    D=dist
  )

fit <- tt$sample(data = dat, chains = 4,
                 adapt_delta = 0.95,
                 max_treedepth = 12,
                 init = \() {list(u_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  w_t_logit_eta = matrix(qlogis(rbeta(TT*N_C, 1, 9)), TT, N_C),
                                  p = rbeta(1, 4, 4),
                                  log_beta_diag = rnorm(N_C,0,.001),
                                  #b_od = rgamma(N_C-2, 4, 4),
                                  #b_self = rgamma(N_C, 4, 4),
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


fit$draws("p") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()
fit$draws("beta[3,3]") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()

fit$draws("gamma[12]") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()
fit$draws("rho") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()
fit$draws("decay_rate_space") %>% as_tibble %>% gather() %>% ggplot(aes(y=value,x=rep(1:1000,4),group=key,color=key)) + geom_line()
