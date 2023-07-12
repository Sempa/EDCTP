dt01 <- as.data.frame(cbind(`Unstimulated and 24hr MTB stimulated RANTES (ng/ml)` = c('RANTES 0-5ng/ml',
                                                                        'RANTES 0-5ng/ml',
                                                                        'RANTES >5ng/ml',
                                                                        'RANTES >5ng/ml',
                                                                        'RANTES >10ng/ml',
                                                                        'RANTES >10ng/ml'
),
Category = c('TB/HIV', 'HIV/C', 'TB/HIV', 'HIV/C', 'TB/HIV', 'HIV/C'
),
`Frequency (% patients)` = c(23, 68, 77, 31, 52, 26),
id = rep(c(1,2), 3)
))
library(ggplot2)
library(tidyverse)
# library(m)
plot1 <- ggplot(dt01 %>% 
                  mutate(`Unstimulated and 24hr MTB stimulated RANTES (ng/ml)` = factor(
                    `Unstimulated and 24hr MTB stimulated RANTES (ng/ml)`, levels = c('RANTES 0-5ng/ml',
                                                                                   # 'RANTES 0-5ng/ml',
                                                                                   'RANTES >5ng/ml',
                                                                                   # 'RANTES>5ng/ml',
                                                                                   # 'RANTES >10ng/ml',
                                                                                   'RANTES >10ng/ml'
                  )), 
                         Category = factor(Category, levels = c('TB/HIV', 'HIV/C')),
                  `Frequency (% patients)` = as.numeric(`Frequency (% patients)`)),                                      # Grouped barplot using ggplot2
       aes(x = `Unstimulated and 24hr MTB stimulated RANTES (ng/ml)`,
           y = `Frequency (% patients)`,
           group = id,
           fill = Category)) +
  geom_bar(stat = "identity",
           position = "dodge") +
  scale_fill_manual(values = c('#998ec3', '#f1a340')) +
  scale_y_continuous(labels = , breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90), expand = c(0,0)) +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    legend.title=element_blank()
  )

plot1
tiff("boxplot_Mammy.tiff", units="in", width = 12, height=7, res=300)
plot1
dev.off()