library(here)
library(tidyverse)

df.sim1 <- readRDS(file=here('output', 'sim', 'df.sim1.rds'))
df.sim2 <- readRDS(file=here('output', 'sim', 'df.sim2.rds'))

# figure
df.sim1$Simulation <- 1
df.sim2$Simulation <- 2

rbind(df.sim1, df.sim2) %>%
  mutate(method = recode(method, 
                         ebcd_pl = 'EBCD-pl',
                         ebcd_l = 'EBCD-l',
                         ebmf_pl = 'EBMF-n/pl',
                         ebmf_l = 'EBMF-n/l',
                         l1ppca = 'L1-penalized PCA',
                         spc = 'SPC',
                         ebpca = 'EB-PCA',
                         gpower = 'GPower',
                         pca = 'PCA', 
  ),
  dist.type2 = recode(dist.type,
                      d1 = 'd[1]', 
                      d2 = 'd[2]', 
                      d3 = 'd[3]',
                      dOR= 'd[or]', 
                      dcov = 'd[cov]')) %>%
  mutate(method = factor(method, levels=c('EBCD-pl', 'EBCD-l', 'EBMF-n/pl', 'EBMF-n/l', 'L1-penalized PCA', 'SPC', 
                                          'EB-PCA', 'GPower', 'PCA'))
  ) -> df.sim

# dummy for a dummy slot in the figure
df.sim <- rbind(df.sim, data.frame(method = 'PCA', seed=1, dist.type='d3',
                                   dist=NA, Simulation = 1, dist.type2='d[3]'))


df.sim %>%
  ggplot()+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'bottom')+
  geom_boxplot(aes(y=dist, col=method), outlier.size=0.5)+
  ylab('')+
  facet_wrap(vars(Simulation, dist.type2), scales = 'free_y', nrow=2,
             labeller=labeller(Simulation = label_both,
                               dist.type2 = label_parsed)) -> fig.sim


print(fig.sim)
ggsave(filename=here('output', 'sim', 'fig_sim.pdf'), device='pdf', width=6, height=4)


df.sim %>% 
  filter(Simulation==1) %>%
  group_by(dist.type, method) %>%
  summarize(mean.dist = mean(dist)) %>%
  pivot_wider(names_from = dist.type, values_from = mean.dist) %>%
  arrange(dcov) 

df.sim %>% 
  filter(Simulation==2) %>%
  group_by(dist.type, method) %>%
  summarize(mean.dist = mean(dist)) %>%
  pivot_wider(names_from = dist.type, values_from = mean.dist) %>%
  arrange(dcov) 
