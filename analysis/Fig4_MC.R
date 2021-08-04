library(tidyverse)
library(entropy)
library(ggpubr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

gini <- function(x) {
  x <- unlist(x)
  x <- sort(x, decreasing = TRUE)
  p <- x/sum(x)
  n <- length(p)
  i <- 1:n
  1 - 2 * (sum(p * i) - 1)/(n - 1)
}

####Diversity within solution category ####
load(file="../data/df_solves_cleaned.Rda")
df_solves = df_solves %>% filter(SPECIES=="GRETI", DEMO==FALSE, !is.na(solution),(PATCH %in% c("CP","BB","GW","MP","MA")), experiment !="dial_diffusion")

df_entropy = df_solves %>% ungroup() %>% arrange(time_stamp) %>% group_by(experiment,PATCH,ID,solution_category) %>% droplevels() %>% summarize(entropy=entropy(table(solution)), evenness = 1-gini(table(solution)), obs=n(),number_solutions_used=n_distinct(solution)) %>% filter(obs>=3) %>% mutate(exp_entropy=exp(entropy))

levels(df_entropy$experiment) = c("Complex gen. 1","Complex gen. 2")
levels(df_entropy$solution_category) = c("dial","slide","two-step")

df_entropy %>% ungroup() %>% group_by(solution_category) %>% summarize(median(exp_entropy))

p1 = ggplot(df_entropy,aes(x=solution_category,y=evenness))+
  geom_jitter(alpha=0.3,aes(size=obs),height=0,width=0.1)+
  geom_boxplot(aes(fill=solution_category), outlier.alpha = 0,alpha=0.5)+
  scale_fill_manual(values=c("#78DDFF","#1BDD7C","#660099"))+
  guides(color=FALSE)+
  labs(x="Solution category",y="Individual behavioural diversity (exp(H))",size="Total\nproductions",fill="Solution category")+
  theme_classic()

p2 = ggplot(df_entropy %>% filter(solution_category=="two-step"),aes(y=exp_entropy,x=as.factor(number_solutions_used)))+
  geom_jitter(aes(size=obs),alpha=0.3,height = 0, width=.1)+
  geom_abline(slope=1,linetype=3)+
  geom_boxplot(show.legend = F,alpha=0.5, fill="#660099")+
  scale_fill_manual(values=c("#1BDD7C"))+
  scale_y_continuous(limits=c(0.9,8),breaks = c(1:8))+
  labs(x="Known two-step variants",y="",size="Total\nproductions")+
  theme_classic()
g1 = ggarrange(p1,p2, labels = "AUTO", nrow = 1, common.legend = T,legend="right")
g1
ggsave(g1, file="../output/Fig4.png",width=11.8,height=6,units="cm",scale=2)
