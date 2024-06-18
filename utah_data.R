library(tidyverse)

case = read_csv("./data/time_series_covid19_confirmed_US.csv") %>% filter(Province_State=="Utah")
test = read_csv("./data/time_series_covid19_US_testing_by_state.csv") %>% filter(state=="UT")
vax = read_csv("./data/COVID-19_Vaccinations_in_the_United_States_County.csv") %>% filter(Recip_State=="UT")
test=test %>% mutate(date=mdy(date)) %>% arrange(date) %>% select(date,tests_combined_total)
vax=vax %>% mutate(date=mdy(Date)) %>% arrange(date) %>% select(date,vax=Series_Complete_Yes) %>% 
  group_by(date) %>% summarize(vax=sum(vax))

case=case %>% select(-c(1:5,7:11)) %>% gather("date","cases",-Admin2) %>% mutate(date=mdy(date)) %>% group_by(Admin2,date) %>% summarize(cases=sum(as.numeric(cases)))

#case=test %>% right_join(case,by="date") %>% left_join(vax,by="date") %>% mutate(vax=replace_na(vax,0))

case=case %>% ungroup() %>% mutate(date=(date = 14 * (as.numeric(date - min(date)) %/% 14) + min(date))) %>% 
  group_by(date,Admin2) %>% summarize(cases=sum(cases-lag(cases),na.rm=T)) %>%#,tests=sum(tests_combined_total-lag(tests_combined_total),na.rm=T))#,
                               #vax=sum(vax-lag(vax),na.rm=T))
                               filter(date>=ymd("2020-03-06")) %>% filter(date<ymd("2021-06-02"))
                               #filter(date>=ymd("2020-10-06")) %>% filter(date<ymd("2021-02-02"))




case %>% group_by(Admin2) %>% filter(any(cases>100000)) %>%
ggplot(aes(x=date,y=cases)) + facet_wrap(~Admin2) + geom_line() + theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
