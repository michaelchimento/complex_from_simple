library(tidyverse)
library(entropy)
library(ggridges)

#### compare real vs randomized data for BB and CP
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### INDIVIDUAL LEVEL SIMS ####
load(file="../data/df_solves_cleaned.Rda")

df_solves = df_solves %>% filter(SPECIES=="GRETI",
                                 DEMO==FALSE,
                                 experiment=="nextgen",
                                 !is.na(solution),
                                 PATCH %in% c("CP","BB"),
                                 solution_category=="complex") %>% group_by(ID) %>% filter(max(ind_solve_index_bycat_byexp)>= 3) %>% droplevels()

#create real entropy scores for individuals
df_observed_ind = df_solves %>%
  select(PATCH,ID,solution) %>%
  group_by(PATCH,ID) %>%
  summarize(entropy=entropy(table(solution)),obs=n(),number_solutions_used=n_distinct(solution),num_solvers = n_distinct(ID)) %>%
  mutate(exp_entropy=exp(entropy), pool="Individual", SITE="Individual (observed)",simulated="observed") %>% dplyr::select(PATCH,ID,exp_entropy,pool,SITE,simulated)

#simulate entropy scores
df_sim_ind <- data.frame(PATCH=character(0),ID=character(0),exp_entropy=integer(0),pool=character(0))
reps = 10000
for (i in 1:reps){
  simulated_data = df_solves %>%
     select(PATCH,ID,solution) %>%
     group_by(PATCH,ID) %>%
     mutate(solution = sample(unique(solution), size=n(), replace = T)) %>%
     summarize(entropy=entropy(table(solution)),obs=n(),number_solutions_used=n_distinct(solution), num_solvers = n_distinct(ID)) %>%
     mutate(exp_entropy=exp(entropy), pool="Individual", SITE="Individual (simulated)",simulated="simulated") %>% dplyr::select(PATCH,ID,exp_entropy,pool,SITE,simulated)
   df_sim_ind = bind_rows(df_sim_ind, simulated_data)
 }
save(df_sim_ind,file="../data/df_sim_ind.Rda")

#### SITE LEVEL SIMS####
load(file="../data/df_solves_cleaned.Rda")

df_solves = df_solves %>% filter(SPECIES=="GRETI",
                                 DEMO==FALSE,
                                 experiment=="nextgen",
                                 !is.na(solution),
                                 PATCH %in% c("CP","BB"),
                                 solution_category=="complex") %>% group_by(ID) %>% filter(max(ind_solve_index_bycat_byexp)>= 3) %>% droplevels()
possible_solutions= levels(df_solves$solution)

df_observed_site = df_solves %>%
  ungroup() %>%
  arrange(time_stamp) %>%
  group_by(PATCH,SITE) %>%
  droplevels() %>%
  summarize(entropy=entropy(table(solution))) %>%
  mutate(exp_entropy=exp(entropy), pool="Site", simulated="observed")

df_sim_site <- data.frame(PATCH=character(0),SITE=character(0),exp_entropy=integer(0),pool=character(0),simulated=character(0))
for (i in 1:reps){
  simulated_data = df_solves %>%
    group_by(ID,solution) %>%
    mutate(solution=sample(possible_solutions[!possible_solutions %in% solution],size=1)) %>%
    ungroup() %>%
    group_by(PATCH,SITE) %>%
    droplevels() %>%
    summarize(entropy=entropy(table(solution))) %>%
    mutate(exp_entropy=exp(entropy), pool="Site",simulated="simulated")

  df_sim_site = bind_rows(df_sim_site,simulated_data)
}
save(df_sim_site,file="../data/df_sim_site.Rda")

#### SUBPOP LEVEL SIMS####
load(file="../data/df_solves_cleaned.Rda")
df_solves = df_solves %>% filter(SPECIES=="GRETI",
                                 DEMO==FALSE,
                                 experiment=="nextgen",
                                 !is.na(solution),
                                 PATCH %in% c("CP","BB"),
                                 solution_category=="complex") %>% group_by(ID) %>% filter(max(ind_solve_index_bycat_byexp)>= 3) %>% droplevels()
possible_solutions= levels(df_solves$solution)
df_observed_subpop = df_solves %>%
  ungroup() %>%
  arrange(time_stamp) %>%
  group_by(PATCH) %>%
  droplevels() %>%
  summarize(entropy=entropy(table(solution))) %>%
  mutate(exp_entropy=exp(entropy), pool="Subpopulation", simulated="observed")

df_sim_subpop <- data.frame(PATCH=character(0),exp_entropy=integer(0),pool=character(0),simulated=character(0))
for (i in 1:reps){
  simulated_data = df_solves %>%
    group_by(ID,solution) %>%
    mutate(solution=sample(possible_solutions[!possible_solutions %in% solution],size=1)) %>%
    ungroup() %>%
    group_by(PATCH) %>%
    droplevels() %>%
    summarize(entropy=entropy(table(solution))) %>%
    mutate(exp_entropy=exp(entropy), pool="Subpopulation",simulated="simulated")

  df_sim_subpop = bind_rows(df_sim_subpop,simulated_data)
}
save(df_sim_subpop,file="../data/df_sim_subpop.Rda")


#### load data to create figure ####
load("../data/df_sim_ind.Rda")
load("../data/df_sim_site.Rda")
load("../data/df_sim_subpop.Rda")


df_sim = bind_rows(df_sim_subpop, df_sim_site,df_sim_ind,df_observed_ind)
df_sim$SITE = as.factor(ifelse(is.na(df_sim$SITE), "Subpopulation", df_sim$SITE))
levels(df_sim$SITE) <- c("Site1","Site2","Site3","Site1","Site2","Site3","Individual (observed)","Individual (simulated)","Subpopulation")


df_obs = bind_rows(df_observed_ind,df_observed_site,df_observed_subpop)
df_obs$SITE = as.factor(ifelse(is.na(df_obs$SITE), "Subpopulation", df_obs$SITE))
levels(df_obs$SITE) <- c("Site1","Site2","Site3","Site1","Site2","Site3","Individual (observed)","Subpopulation")

df_sim$pool = factor(df_sim$pool, levels=c("Individual","Site","Subpopulation"))
df_obs$pool = factor(df_obs$pool, levels=c("Individual","Site","Subpopulation"))

df_text = bind_rows(df_obs,df_sim) %>% filter(pool!="Individual") %>% ungroup() %>% droplevels() %>% mutate(PATCH=as.factor(PATCH))
df_text = df_text %>% group_by(pool,PATCH,SITE) %>% filter(is.na(ID)) %>% mutate(is_lessthan = exp_entropy <= exp_entropy[simulated=="observed"]) %>% summarize(exp_entropy=exp_entropy[simulated=="observed"], tail_prob = sum(is_lessthan)/10000)


p6 = ggplot()+
  facet_grid(pool~PATCH,scales="free_y")+
  stat_density_ridges(data=df_sim,
                      aes(x=exp_entropy,
                          y=SITE,
                          fill = 0.5 - abs(0.5-stat(ecdf))),
                      calc_ecdf = TRUE,
                      geom="density_ridges_gradient",
                      quantiles = c(0.05, 0.95),
                      scale=.9) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1,option="A") +
  geom_point(data=df_obs, aes(x=exp_entropy, y=as.character(SITE)))+
  geom_text(data=df_text, aes(x=exp_entropy, y=as.character(SITE), label= paste0(round(exp_entropy,2)," (",round(tail_prob,2),")")),vjust=2,hjust=.6)+
  labs(x="Behavioural diversity (exp(H))", y = "Level of analysis")+
  scale_x_continuous(limits=c(.99,8.01),breaks=c(1:8))+
  theme_classic()

#complex gen 2
load(file="../data/df_solves_cleaned.Rda")
df_solves = df_solves %>% filter(SPECIES=="GRETI",
                                 DEMO==FALSE,
                                 experiment=="nextgen",
                                 !is.na(solution),
                                 PATCH %in% c("CP","BB"),
                                 solution_category=="complex") %>% group_by(ID) %>% filter(max(ind_solve_index_bycat_byexp)>= 3) %>% droplevels()

levels(df_solves$experiment) = c("Complex Gen. 2")
#levels(df_solves$SITE) = c("Site1","Site2","Site3","Site1","Site2","Site3")
df_solves$solution = factor(df_solves$solution, levels=c("dial_l_slide_l","slide_l_dial_l", "dial_l_slide_r","slide_r_dial_l", "dial_r_slide_l","slide_l_dial_r","dial_r_slide_r","slide_r_dial_r"))

df_sum = df_solves %>% group_by(PATCH,solution) %>% summarize(soln_count=n()) %>% ungroup()
df_sum = df_sum %>% group_by(PATCH) %>% mutate(day_count=sum(soln_count)) %>% ungroup()
df_sum = df_sum %>% mutate(prop=soln_count/day_count)

df_sum = df_sum %>% ungroup() %>%
  group_by(PATCH) %>%
  arrange(PATCH, 1-prop) %>%
  mutate(order=row_number())

p7=ggplot(df_sum,aes(x=as.factor(order),y=prop,fill=solution)) +
  facet_grid(~PATCH,scales="free_x")+
  geom_bar(stat="identity")+
  scale_fill_viridis_d(option="D", direction = -1,end=.9)+
  scale_alpha(range=c(0.5,1))+
  labs(x="Solution variants", y="Prop.",fill="Variant")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
library(cowplot)
g1=plot_grid(p6,p7,ncol=1,labels=c("A","B"),align = "v",axis = "l",rel_heights = c(1,.4))

ggsave(g1,file="../output/Fig5.pdf",width=11.8,height=12,units="cm",scale=2)
