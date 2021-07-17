#Additional Statistical Analyses

#load packages
library(tidyverse)
library(entropy)
library(ggpubr)
library(lme4)
library(lmerTest)
library(arm)
library(stargazer)

#set options
options(scipen=5)
options(digits=5)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Section C: Prior knowledge state for acquiring the complex behaviour ####
### Experiment 2 ###
load(file="../data/df_solves_cleaned.Rda")
df_solves = df_solves %>% filter(SPECIES=="GRETI", DEMO==FALSE, !is.na(solution),total_solutions>=3, !(PATCH %in% c("BO")), experiment=="complex_prog") %>% droplevels()
load(file="../data/df_visitors_nonsolvers.Rda")
df_visitors = df_visitors %>% filter(PATCH!="BO",experiment=="complex_prog") %>% ungroup()

df_exp2 = df_solves %>%
  filter(total_solutions_byexp>=3) %>%
  mutate(soln_cat = solution_category) %>%
  pivot_wider(names_from = solution_category, values_from=ind_solve_index_bycat_byexp) %>%
  group_by(ID) %>%
  fill(dial,slide,complex) %>%
  mutate(dial=ifelse(is.na(dial),0,dial),
         slide=ifelse(is.na(slide),0,slide),
         complex=ifelse(is.na(complex),0,complex))

#create dataframe of birds who did not learn complex by the end of experiment 2
df_exp2_complex_dummies = df_exp2 %>% group_by(ID) %>% slice_tail(n=1) %>% filter(complex<3) %>% dplyr::select(ID,dial,slide,complex,experiment) %>% mutate(learned_complex=FALSE)
#create dataframe of birds who visited but did not produce any solutions
df_exp2_ultra_dummies = df_visitors %>% ungroup() %>% dplyr::select(ID,experiment) %>% distinct() %>% mutate(dial=0,slide=0,complex=0,learned_complex=FALSE)
n_distinct(df_visitors$ID)
#combine all dummies together
df_exp2_complex_dummies = bind_rows(df_exp2_complex_dummies,df_exp2_ultra_dummies)

#label birds who produced more than 3 complex solutions as knowledgeable
df_exp2_complex_smart = df_exp2 %>% group_by(ID) %>% filter(max(complex)>=3, soln_cat=="complex") %>% slice_head(n=1) %>% dplyr::select(ID,dial,slide,complex,experiment) %>% mutate(learned_complex=TRUE)

#combine lists of non-learners and learners
df_dumb_smart = bind_rows(df_exp2_complex_dummies,df_exp2_complex_smart)

#label knowledge states for simple solution components, with cutoff of 3 solutions for each category
df_dumb_smart = df_dumb_smart %>%
  mutate(knowledgable_neither=ifelse(dial<3 & slide<3,TRUE,FALSE),
         knowledgable_dial=ifelse(dial>=3,TRUE,FALSE),
         knowledgable_slide=ifelse(slide>=3,TRUE,FALSE),
         knowledgable_both = ifelse(knowledgable_dial & knowledgable_slide, TRUE, FALSE)) %>%
  mutate(knowledgable_only_dial=ifelse(knowledgable_dial & !knowledgable_both,TRUE, FALSE),
         knowledgable_only_slide=ifelse(knowledgable_slide & !knowledgable_both,TRUE, FALSE))

#save for later analysis
df_exp2 = df_dumb_smart

num_solvers = n_distinct(df_dumb_smart$ID)

table_exp2 = df_dumb_smart %>% group_by(learned_complex=as.character(learned_complex)) %>% summarize(neither = sum(knowledgable_neither), prop.0 = neither/num_solvers, only_dial = sum(knowledgable_only_dial), prop.1=only_dial/num_solvers, only_slide = sum(knowledgable_only_slide), prop.2=only_slide/num_solvers, both = sum(knowledgable_both), prop.3=both/num_solvers, Total=neither+only_dial+only_slide+both) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

table_exp2
write.table(table_exp2, file = "../output/TableS2Exp2.txt", sep = ",", quote = FALSE, row.names = F)



### Experiment 3 ###
#save a list of IDs which learned complex during exp2
df_exp2_complex_solvers = df_dumb_smart %>% filter(learned_complex==TRUE) %>% distinct(ID)

#reload solvers DF
load(file="../data/df_solves_cleaned.Rda")
df_solves = df_solves %>% filter(SPECIES=="GRETI", DEMO==FALSE, !is.na(solution),total_solutions>=3, !(PATCH %in% c("BO","MA")), experiment=="nextgen") %>% droplevels()
load(file="../data/df_visitors_nonsolvers.Rda")
df_visitors = df_visitors %>% filter(!(PATCH %in% c("BO","MA")),experiment=="nextgen") %>% ungroup()

df_exp3 = df_solves %>%
  filter(total_solutions_byexp>=3) %>%
  mutate(soln_cat = solution_category) %>%
  pivot_wider(names_from = solution_category, values_from=ind_solve_index_bycat_byexp) %>%
  group_by(ID) %>%
  fill(dial,slide,complex) %>%
  mutate(dial=ifelse(is.na(dial),0,dial),
         slide=ifelse(is.na(slide),0,slide),
         complex=ifelse(is.na(complex),0,complex))

#create dataframe of birds who did not learn complex by the end of experiment 2
df_exp3_complex_dummies = df_exp3 %>% group_by(ID) %>% slice_tail(n=1) %>% filter(complex<3) %>% dplyr::select(ID,dial,slide,complex,experiment) %>% mutate(learned_complex=FALSE)
#create dataframe of birds who visited but did not produce any solutions
df_exp3_ultra_dummies = df_visitors %>% ungroup() %>% dplyr::select(ID,experiment) %>% distinct() %>% mutate(dial=0,slide=0,complex=0,learned_complex=FALSE)
n_distinct(df_visitors$ID)
#combine all dummies together
df_exp3_complex_dummies = bind_rows(df_exp3_complex_dummies,df_exp3_ultra_dummies)

#label birds who produced more than 3 complex solutions as knowledgeable, remove those who had already learned during exp2
df_exp3_complex_smart = df_exp3 %>% group_by(ID) %>% filter(max(complex)>=3, soln_cat=="complex") %>% slice_head(n=1) %>% filter(!(ID %in% df_exp2_complex_solvers$ID)) %>% dplyr::select(ID,dial,slide,complex,experiment) %>% mutate(learned_complex=TRUE)

df_dumb_smart = bind_rows(df_exp3_complex_dummies,df_exp3_complex_smart)

df_dumb_smart = df_dumb_smart %>%
  mutate(knowledgable_neither=ifelse(dial<3 & slide<3,TRUE,FALSE),
         knowledgable_dial=ifelse(dial>=3,TRUE,FALSE),
         knowledgable_slide=ifelse(slide>=3,TRUE,FALSE),
         knowledgable_both = ifelse(knowledgable_dial & knowledgable_slide, TRUE, FALSE)) %>%
  mutate(knowledgable_only_dial=ifelse(knowledgable_dial & !knowledgable_both,TRUE, FALSE),
         knowledgable_only_slide=ifelse(knowledgable_slide & !knowledgable_both,TRUE, FALSE))

df_exp3 = df_dumb_smart

num_solvers = n_distinct(df_dumb_smart$ID)

table_exp3 = df_dumb_smart %>% group_by(learned_complex=as.character(learned_complex)) %>% summarize(neither = sum(knowledgable_neither), prop.0 = neither/num_solvers, only_dial = sum(knowledgable_only_dial), prop.1=only_dial/num_solvers, only_slide = sum(knowledgable_only_slide), prop.2=only_slide/num_solvers, both = sum(knowledgable_both), prop.3=both/num_solvers, Total=neither+only_dial+only_slide+both) %>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total")))

table_exp3
write.table(table_exp3, file = "../output/TableS2_exp_3.txt", sep = ",", quote = FALSE, row.names = F)


df_combined = bind_rows(df_exp2,df_exp3) %>% droplevels()
m1 = glmer(learned_complex ~ knowledgable_dial*knowledgable_slide + (1 | experiment), data=df_combined, family="binomial")
m2 = glm(learned_complex ~ knowledgable_dial*knowledgable_slide, data=df_combined, family="binomial")
anova(m1,m2)

log_to_prob <- function(log_odds){
  odds <- exp(log_odds)
  prob <- odds / (1 + odds)
  return(prob)
}

summary(m2)
stargazer(m2, dep.var.labels = c("Learned complex"), covariate.labels = c("Knowledgable (dial only)", "Knowledgable (slide only)", "Knowledgable (both D,S)", "Naive"), title="Prior knowledge state to learning complex", font.size = "small",report="vcstp*",single.row=T,out="../output/TableS3.html")

#intercept
log_to_prob(coef(m2)[1])
log_to_prob(coef(m2)[1]-se.coef(m2)[1])
log_to_prob(coef(m2)[1]+se.coef(m2)[1])

#dial only
log_to_prob(coef(m2)[1]+coef(m2)[2])
log_to_prob(coef(m2)[1]+se.coef(m2)[1]+coef(m2)[2]+se.coef(m2)[2])
log_to_prob(coef(m2)[1]-se.coef(m2)[1]+coef(m2)[2]-se.coef(m2)[2])
#slide only
log_to_prob(coef(m2)[1]+coef(m2)[3])
log_to_prob(coef(m2)[1]+se.coef(m2)[1]+coef(m2)[3]+se.coef(m2)[3])
log_to_prob(coef(m2)[1]-se.coef(m2)[1]+coef(m2)[3]-se.coef(m2)[3])

#knowledge of both
log_to_prob(coef(m2)[1]+coef(m2)[2]+coef(m2)[3]+coef(m2)[4])
log_to_prob(coef(m2)[1]+se.coef(m2)[1]+coef(m2)[2]+se.coef(m2)[2]+coef(m2)[3]+se.coef(m2)[3]+coef(m2)[4]+se.coef(m2)[4])
log_to_prob(coef(m2)[1]-se.coef(m2)[1]+coef(m2)[2]-se.coef(m2)[2]+coef(m2)[3]-se.coef(m2)[3]+coef(m2)[4]-se.coef(m2)[4])

#### Section D: Solution choice within individuals, sites and subpopulations ####
load(file="../data/df_solves_cleaned.Rda")
summary(df_solves)
df_solves = df_solves %>% filter(SPECIES=="GRETI", DEMO==FALSE, !is.na(solution), (PATCH %in% c("CP","BB","GW","MP","MA")), experiment !="dial_diffusion")
df_solves = df_solves %>%
  mutate(soln_cat = solution_category) %>%
  pivot_wider(names_from = soln_cat, values_from=ind_solve_index_bycat_byexp) %>%
  group_by(ID) %>%
  fill(dial,slide,complex) %>%
  mutate(dial=ifelse(is.na(dial),0,dial),
         slide=ifelse(is.na(slide),0,slide),
         complex=ifelse(is.na(complex),0,complex)) %>% mutate(learned_complex=ifelse(complex>0,TRUE,FALSE))

df_entropy = df_solves %>% ungroup() %>% arrange(time_stamp) %>% group_by(experiment,PATCH,ID,learned_complex,solution_category) %>% droplevels() %>% summarize(entropy=entropy(table(solution)),obs=n(),number_solutions_used=n_distinct(solution)) %>% filter(obs>=3) %>% mutate(exp_entropy=exp(entropy))
levels(df_entropy$experiment) = c("Complex gen. 1","Complex gen. 2")

#no significant effect of subpop on solution category
model_learned_complex = lmer(exp_entropy ~ learned_complex*solution_category + (1| experiment) + (1 | ID), data=df_entropy)
summary(model_learned_complex)

#no significant effect of subpop on solution category
model_patch = lmer(exp_entropy ~ PATCH*solution_category + (1 | ID), data=df_entropy)
summary(model_patch)

#no significant effect of experiment on solution category
model_experiment = lmer(exp_entropy ~ experiment*solution_category + (1 | ID), data=df_entropy)
summary(model_experiment)

#ok, given that there are no sig differences, let's remake the dataframe
#behavioural diversity within solution category
load(file="../data/df_solves_cleaned.Rda")
summary(df_solves)
df_solves = df_solves %>% filter(SPECIES=="GRETI", DEMO==FALSE, !is.na(solution), (PATCH %in% c("CP","BB","GW","MP","MA")), experiment !="dial_diffusion")
df_entropy = df_solves %>% ungroup() %>% arrange(time_stamp) %>% group_by(experiment,PATCH,ID,solution_category) %>% droplevels() %>% summarize(entropy=entropy(table(solution)),obs=n(),number_solutions_used=n_distinct(solution)) %>% filter(obs>=3) %>% mutate(exp_entropy=exp(entropy))
model_solution_category = lmer(exp_entropy ~ solution_category + (1| experiment) + (1 | ID), data=df_entropy)
summary(model_solution_category)

#point estimate and SE for complex solvers
fixef(model_solution_category)[1] + fixef(model_solution_category)[3]
fixef(model_solution_category)[1]+se.coef(model_solution_category)$fixef[1] + fixef(model_solution_category)[3]+se.coef(model_solution_category)$fixef[3]
fixef(model_solution_category)[1]-se.coef(model_solution_category)$fixef[1] + fixef(model_solution_category)[3]-se.coef(model_solution_category)$fixef[3]

#effect of number known variants on behavioural diversity
model_df = df_entropy %>% filter(solution_category=="complex", number_solutions_used>1) %>% mutate(number_solutions_used = number_solutions_used-2)
model_behav_variants = lmer(exp_entropy ~ number_solutions_used + (1| experiment), data=model_df)

#export to stargazer
class(model_solution_category) <- "lmerMod"
class(model_behav_variants) <- "lmerMod"

stargazer(model_solution_category,model_behav_variants, dep.var.labels = c("Behavioural Diversity","Behavioural Diversity"), covariate.labels = c("Solution category (slide)", "Solution category (complex)", "Known complex variants","Intercept"), title="Behavioural Diversity Models", font.size = "small",report="vcstp*",single.row=T,out="../output/TableS4.html")

#what proportion of complex solutions contain an individuals most frequently produced simple solution?
load(file="../data/df_solves_cleaned.Rda")

df_complex_solvers = df_solves %>% filter(SPECIES=="GRETI",
                                 DEMO==FALSE,
                                 experiment!="dial_diffusion",
                                 !is.na(solution),
                                 PATCH %in% c("CP","BB","GW","MA","MP"),
                                 solution_category=="complex") %>%
  group_by(experiment,ID) %>%
  filter(max(ind_solve_index_bycat_byexp)>= 3) %>%
  ungroup() %>%
  droplevels() %>% dplyr::select(ID) %>% distinct()

df_solves = df_solves %>% filter(SPECIES=="GRETI",
                     DEMO==FALSE,
                     experiment!="dial_diffusion",
                     !is.na(solution),
                     PATCH %in% c("CP","BB","GW","MA","MP"))

df = df_solves %>% filter(ID %in% df_complex_solvers$ID)
df = df %>%
  mutate(soln_cat = solution_category) %>%
  pivot_wider(names_from = soln_cat, values_from=ind_solve_index_bycat_byexp) %>%
  group_by(ID) %>%
  fill(dial,slide,complex) %>%
  mutate(dial=ifelse(is.na(dial),0,dial),
         slide=ifelse(is.na(slide),0,slide),
         complex=ifelse(is.na(complex),0,complex)) %>%
  ungroup() %>%
  dplyr::select(experiment,ID,solution_category,solution,dial,slide,complex)

df_simple_pref = df %>% filter(complex==0) %>% group_by(experiment,ID,solution_category,solution) %>% summarize(count=n()) %>% ungroup()

## First look at simple preference, ignoring whether it's dial or slide
df_simple_pref = df_simple_pref %>% filter(solution_category!= "complex") %>% group_by(experiment,ID) %>% summarize(simple_pref = solution[count==max(count)])
df = df %>% filter(solution_category=="complex")
df = left_join(df,df_simple_pref)

df_prop = df %>% ungroup() %>% rowwise() %>% mutate(is_component= grepl(simple_pref, solution,fixed = TRUE))
df_prop = df_prop %>% group_by(experiment,ID) %>% summarize(numerator=sum(is_component),denom=n(), prop_true=numerator/denom) %>% filter(!is.na(numerator))

sum(df_prop$numerator,na.rm = T)/sum(df_prop$denom, na.rm = T)

ggplot(df_prop,aes(x=experiment,y=prop_true))+
  geom_jitter()+
  geom_boxplot(alpha=0.5, outlier.alpha = 0)+
  labs(y="Prop. complex solutions containing preffered simple solution")

# Now look at simple preference within individuals that have produced both dial AND slide, and see how many of those preferences are found in their favorite solutions
df_all_knowledge_states = bind_rows(df_exp2,df_exp3)
df_knewboth = df_all_knowledge_states %>% filter(knowledgable_both, learned_complex) %>% select(ID)


df_simple_pref = df_simple_pref %>% droplevels() %>% filter(ID %in% df_knewboth$ID, solution_category!= "complex") %>% group_by(experiment,ID,solution_category) %>% summarize(simple_pref = solution[count==max(count)][1])

df_simple_pref %>% print(n=Inf)

df_simple_pref = df_simple_pref %>% mutate(simple_pref=as.character(simple_pref)) %>% pivot_wider(values_from = simple_pref, names_from = solution_category) %>% drop_na()

df = df %>% filter(solution_category=="complex") %>% select(!c(dial,slide,complex))
df = left_join(df,df_simple_pref)

df = df %>% drop_na()

df_prop = df %>% ungroup() %>% rowwise() %>% mutate(is_component_dial= grepl(dial, solution,fixed = TRUE), is_component_slide= grepl(slide, solution,fixed = TRUE), is_component=is_component_dial & is_component_slide)
df_prop = df_prop %>% group_by(experiment,ID) %>% summarize(numerator=sum(is_component),denom=n(), prop_true=numerator/denom) %>% filter(!is.na(numerator))

sum(df_prop$numerator,na.rm = T)/sum(df_prop$denom, na.rm = T)


df_knewboth
