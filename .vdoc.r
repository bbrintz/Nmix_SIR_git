#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
library(tidyverse)
library(patchwork)
#
#
#

case = read_csv("./data/time_series_covid19_confirmed_US.csv") %>% filter(`Province_State`=="Utah")
case=case %>% dplyr::select(-c(1:5,7:11)) %>% gather("date","cases",-Admin2) %>% mutate(date=mdy(date)) %>% group_by(Admin2,date) %>% summarize(cases=sum(as.numeric(cases)))
case=case %>% ungroup() %>% mutate(date=(date = 7 * (as.numeric(date - min(date)) %/% 7) + min(date))) %>% 
  group_by(date,Admin2) %>% summarize(cases=sum(cases)) %>%
                               filter(date>=ymd("2020-04-15")) %>% filter(date<ymd("2021-06-02"))
case=case %>% group_by(Admin2) %>% mutate(new_cases=cases-lag(cases)) %>% arrange(Admin2,date) %>% filter(!is.na(new_cases))# %>% filter(Admin2 %in% readRDS("final_counties.rds")) %>%# filter(any(cases>1000)) %>%

dat_final=case %>% filter(Admin2 %in% readRDS("final_counties.rds"))#group_by(Admin2) %>% filter(any(cases>100))

pop=read_csv("./data/utah_counties_pop_coord.csv") %>% arrange(desc(Population_2020))

d1=dat_final %>% filter(date<ymd("2020/8/19")) %>% 
ungroup() %>% rename(County="Admin2") %>% left_join(pop,by="County") %>% arrange(desc(Population_2020),date) %>%
dplyr::select(-Population_2020,-Latitude,-Longitude,-cases) %>% pivot_wider(names_from=County,values_from=new_cases)

#fit=readRDS("./output_4/fit_10chn.rds")
#fit2=readRDS("./output_11/fit_10chn.rds")

counties=c(1,2,3,8)

get_plot <- function(fitsize,counties){
fit=readRDS(paste0("./output_",fitsize,"/fit_10chn.rds"))
z_t_d <- fit$draws("si_t", format = "draws_array") |> posterior::as_draws_rvars()
z_t_d <- z_t_d$si_t
qpt025 <- quantile(z_t_d,0.025)
qpt975 <- quantile(z_t_d,0.975)

pop4=pop$Population_2020[counties]

colnames(d1)[counties+1]

obs=rbind(data.frame(value="mean",sweep(mean(z_t_d)[-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
data.frame(value="lwr",sweep(qpt025[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")),
data.frame(value="upr",sweep(qpt975[1,-1,],MARGIN = 2, STATS = pop4, FUN = "*")))

colnames(obs)[-1] = colnames(d1)[counties+1]


obs=obs %>% gather(County,Cases,-value) %>% pivot_wider(names_from=value,values_from=Cases) %>% unnest() %>% mutate(date=rep(d1$date,fitsize))


d1[,c(1,counties+1)] %>%  gather(County,Cases,-date) %>% 
ggplot(aes(x=date,y=Cases,color=County)) + geom_point() + facet_wrap(~County) + 
geom_line(data=obs,aes(x=date,y=mean),color="black") + geom_line(data=obs,aes(x=date,y=lwr),linetype=2,color="black") + 
geom_line(data=obs,aes(x=date,y=upr),linetype=2,color="black") +
 theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

a=get_plot(4,c(1,2,3,8))
b=get_plot(11,c(1,2,3,4,5,6,7,8,10,11,12))

a/b

#
#
#
#
#
#

fit=readRDS(paste0("./output_",4,"/fit_10chn.rds"))
fit2=readRDS(paste0("./output_",11,"/fit_10chn.rds"))


a=ggplot() + 
geom_histogram(data=data.frame(p=rbeta(1e6,0.5 * 4, 0.5 * 4)),
aes(x=p,..density..),alpha=.25,fill="grey",linetype=2) +
geom_histogram(data=fit$draws("p") %>% as_tibble %>% gather() %>% mutate(key=factor(key)), 
aes(x=value,y=..density..,fill=1),alpha=0.75) + 
scale_fill_viridis_c() +
theme_minimal() +
theme(legend.position="none")

b=ggplot() + 
geom_histogram(data=data.frame(p=rbeta(1e6,0.5 * 4, 0.5 * 4)),
aes(x=p,..density..),alpha=.25,fill="grey",linetype=2) +
geom_histogram(data=fit2$draws("p") %>% as_tibble %>% gather() %>% mutate(key=factor(key)), 
aes(x=value,y=..density..,fill=1),alpha=0.75) + 
scale_fill_viridis_c() +
theme_minimal() +
theme(legend.position="none")

a+b
#
#
#
#
#
#
#
#
#
#
