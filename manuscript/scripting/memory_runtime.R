library(tidyverse)
library(ggplot2)

setwd('/Users/bradentemp/Library/CloudStorage/GoogleDrive-braden@twofrontiers.org/Shared drives/The Two Frontiers Project/STUDIES/PUBLICATION_2024-11-12_xtree_pvc')

kraken = read.delim('data_packet/benchmarktiming/ram_usage_kraken2',header=F,sep=',') %>% mutate(type = 'Kraken2')
kraken$V2 = kraken$V2 - kraken$V2[[1]]
kraken$V2[kraken$V2<0]=0
kraken$V1 = kraken$V1 - kraken$V1[[1]]
kraken = kraken %>% mutate(type = 'kraken (3 processes)',type = if_else(V1 > 9868 & V1<17993 ,'kraken (2 processes)',if_else(V1>=17993,'kraken (1 process)',type)))
kraken = kraken %>% group_by(type) %>% mutate(minval = min(V1)) %>% mutate(V1 = V1 - minval)

xtree = read.delim('data_packet/benchmarktiming/ram_usage_xtree_mb',header=F,sep=',')
xtree$V2 = xtree$V2 - xtree$V2[[1]] + 260000
#xtree$V1 = xtree$V1 - xtree$V1[[1]]
#xtree = xtree %>% filter(V1>=28)
xtree = xtree %>% group_by(V3) %>% mutate(minval = min(V1)) %>% mutate(V1 = V1 - minval) %>% ungroup %>% filter((V1/60)<40) 
colnames(xtree)[3] = 'type'
xtree = xtree %>% mutate(V2 = if_else(V1 == 0,0,V2)) 
xtree = xtree %>% group_by(type) %>% mutate(maxval = max(V1))
xtree = xtree %>% mutate(V2 = if_else(V1 == maxval,0,V2)) %>% select(-maxval)

dat = bind_rows(kraken,xtree)
dat$V2 = dat$V2/1000
dat$V1 = dat$V1/60
colnames(dat) = c('Minutes','Memory (GB)','type','minval')

dat$type = gsub('kraken','Kraken2',dat$type)
dat$type = gsub('Xtree','XTree',dat$type)

ggplot(dat %>% mutate(type2 = if_else(grepl('XTree',type),'XTree','Kraken2')), aes(x = Minutes, y = `Memory (GB)`, color = type)) + geom_line(linewidth = 1.5, alpha = 0.9) + theme_cowplot() + scale_color_manual(values = c("Kraken2 (1 process)" = "#8B0000", "Kraken2 (2 processes)" = "#B22222", "Kraken2 (3 processes)" = "#FF6347", "XTree (1 process)" = "#00008B", "XTree (2 processes)" = "#4169E1", "XTree (3 processes)" = "#87CEFA")) + ggtitle('Performance with parallelization') + theme(legend.position = 'bottom') + facet_grid(type2 ~ .)
ggsave('plots/performance.pdf',width=6,height=4)



#ggplot(xtree,aes(x = V1/100,y=V2/1000,fill = V3, color= V3)) + geom_line(linewidth=1) + theme_cowplot() + scale_fill_brewer(palette = 'Set3') + ggtitle('Performance with parallelization') + theme(legend.position = 'bottom')









