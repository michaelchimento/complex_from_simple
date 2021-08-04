library(tidyverse)
library(ggpubr)
library(ggnewscale)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#all together with heights proportional to solutions
load(file="../data/df_solves_cleaned.Rda")
df_solves = df_solves %>%
  group_by(ID) %>%
  mutate(DEMO=ifelse(sum(DEMO)>0,TRUE,FALSE),
         solution_category=relevel(solution_category,"complex")) %>%
  ungroup() %>%
  filter(DEMO==FALSE,
         PATCH %in% c("CP","BB","MA","GW","BO","MP")) %>%
  droplevels()

levels(df_solves$experiment) = c("Dial diffusion", "Complex 1st gen.", "Complex 2nd gen.")
levels(df_solves$PATCH) = c("BB","Control","CP","GW","MA","MP")
df_solves = df_solves %>% mutate(PATCH=ifelse(experiment=="Complex 2nd gen." & PATCH=="MA","Control",as.character(PATCH)))
df_solves$PATCH = factor(df_solves$PATCH, levels=c("GW","BB","CP","MP","MA","Control"))

#only include dial/complex solutions
df_solves = df_solves %>% filter((experiment=="Dial diffusion" & solution_category=="dial") | (experiment=="Complex 1st gen." & solution_category=="complex") | (experiment=="Complex 2nd gen." & solution_category=="complex")) %>% ungroup() %>% droplevels()
df_sum$solution = factor(df_sum$solution, levels=c("dial_l_slide_l","slide_l_dial_l", "dial_l_slide_r","slide_r_dial_l", "dial_r_slide_l","slide_l_dial_r","dial_r_slide_r","slide_r_dial_r"))

#exclude those who produced less than 3 solutions per experiment
ID_include = df_sum %>% group_by(experiment,ID) %>% summarize(count = n()) %>% filter(count>=3) %>% ungroup() %>% distinct(ID)
df_sum = df_sum %>% filter(ID %in% ID_include$ID)
ID_include = df_sum_DD %>% group_by(experiment,ID) %>% summarize(count = n()) %>% filter(count>=3) %>% ungroup() %>% distinct(ID)
df_sum_DD = df_sum_DD %>% filter(ID %in% ID_include$ID)

ggplot() +
  facet_grid(experiment~PATCH, scales = "free_y")+
  geom_bar(data=df_sum_DD, position="stack", aes(x=as.factor(ABSDATE),fill=simple_solution))+
  scale_fill_viridis_d(option="A",begin=.4,end=.7, direction=-1,  guide = guide_legend(order = 1))+
  labs(fill="One-step variant")+
  new_scale_fill() +
  geom_bar(data=df_sum,position="stack", aes(x=as.factor(ABSDATE),fill=complex_solution))+
  scale_fill_viridis_d(option="D", direction = -1,end=.9, guide = guide_legend(order = 2))+
  labs(fill="Two-step variant")+
  scale_x_discrete(breaks=c(1,5,10,15,20))+
  labs(x="Experimental day", y="Solutions (count)")+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA))

ggsave("../output/Fig2.pdf",width=17.8,height=9,units="cm",scale=1.5)
