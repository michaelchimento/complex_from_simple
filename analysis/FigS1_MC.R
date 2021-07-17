library(tidyverse)
library(scales)
library(lme4)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


load(file="../data/df_solves_cleaned.Rda")

df_solves = df_solves %>% filter(DEMO==FALSE,SPECIES=="GRETI") %>% group_by(solution) %>% mutate(scaled_ind_solve_index_bytype = scale(ind_solve_index_bytype))

df_model = df_solves %>% filter(DEMO==FALSE,SPECIES=="GRETI")

m1 = lmer(log(TTS+1) ~ solution + (1|SITE / ID), data=df_model)
summary(m1)
exp(fixef(m1))

df_summary = df_solves %>% group_by(solution_category,solution) %>% filter(!is.na(solution),DEMO==FALSE,SPECIES=="GRETI") %>% summarize(freq=n(),avgTTS=mean(TTS,na.rm = T),unique_solvers=n_distinct(ID))

df_summary

ggplot(data=df_summary, aes(x=freq,y=avgTTS))+
  facet_wrap(~solution_category,scales="free_x")+
  geom_point(aes(color=solution),size=2)+
  scale_y_continuous()+
  labs(x="Frequency",y="Average TTS (seconds)")+
  scale_color_viridis_d(direction=-1, end = .9)+
  theme(axis.text.x = element_text(angle = 90),legend.key.width = unit(2, 'cm'),legend.key.height = unit(1, 'cm'),)+
  theme(legend.position="bottom")+
  theme_classic()

ggsave("../output/FigS1.png",width=11.8,height=5,units="cm",scale=2)

df_solves %>% group_by(solution_category) %>% filter(!is.na(solution),DEMO==FALSE,SPECIES=="GRETI") %>% summarize(freq=n(),avgTTS=mean(TTS,na.rm = T),unique_solvers=n_distinct(ID))
