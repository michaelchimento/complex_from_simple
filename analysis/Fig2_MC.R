library(tidyverse)
library(ggpubr)
library(ggnewscale)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#all together with heights proportional to solutions
load(file="../data/df_solves_cleaned.Rda")

df_solves = df_solves %>% group_by(ID) %>% mutate(DEMO=ifelse(sum(DEMO)>0,TRUE,FALSE),solution_category=relevel(solution_category,"complex")) %>% ungroup() %>% filter(DEMO==FALSE,PATCH %in% c("CP","BB","MA","GW","BO","MP"), total_solutions_byexp>=3) %>% droplevels()
levels(df_solves$experiment) = c("Dial diffusion", "Complex 1st gen.", "Complex 2nd gen.")
levels(df_solves$PATCH) = c("BB","Control","CP","GW","MA","MP")
df_solves = df_solves %>% mutate(PATCH=ifelse(experiment=="Complex 2nd gen." & PATCH=="MA","Control",as.character(PATCH)))
df_solves$PATCH = factor(df_solves$PATCH, levels=c("GW","BB","CP","MP","MA","Control"))

df_solves = df_solves %>% filter((experiment=="Dial diffusion" & solution_category=="dial") | (experiment=="Complex 1st gen." & solution_category=="complex") | (experiment=="Complex 2nd gen." & solution_category=="complex")) %>% ungroup() %>% droplevels()

df_sum = df_solves %>% group_by(experiment,PATCH,ABSDATE,solution) %>% summarize(soln_count=n())
df_sum = df_sum %>% group_by(experiment,PATCH,ABSDATE) %>% mutate(day_count=sum(soln_count)) %>% ungroup()
df_sum = df_sum %>% mutate(prop=soln_count/day_count)
df_sum = df_sum %>% group_by(experiment) %>% mutate(max_day_count=max(day_count))

df_sum_DD = df_solves %>% filter(experiment=="Dial diffusion") %>% mutate(simple_solution=solution)
df_sum = df_solves %>% filter(experiment!="Dial diffusion") %>% mutate(complex_solution=solution) %>% droplevels()
df_sum$complex_solution = factor(df_sum$complex_solution, levels=c("dial_l_slide_l","slide_l_dial_l", "dial_l_slide_r","slide_r_dial_l", "dial_r_slide_l","slide_l_dial_r","dial_r_slide_r","slide_r_dial_r"))

ggplot() +
  facet_grid(experiment~PATCH, scales = "free_y")+
  geom_bar(data=df_sum_DD, position="stack", aes(x=as.factor(ABSDATE),fill=simple_solution))+
  scale_fill_viridis_d(option="A",begin=.4,end=.7, direction=-1)+
  new_scale_fill() +
  geom_bar(data=df_sum,position="stack", aes(x=as.factor(ABSDATE),fill=complex_solution))+
  scale_fill_viridis_d(option="D", direction = -1,end=.9)+
  scale_alpha(range=c(0.3,1), breaks=seq(0,1,.1))+
  guides(alpha=F)+
  scale_x_discrete(breaks=c(1,5,10,15,20))+
  labs(x="Experimental day", y="Solutions (count)")+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA))

ggsave("../output/Fig2.pdf",width=17.8,height=9,units="cm",scale=1.5)
