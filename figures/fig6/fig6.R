library(tidyverse)
library(ggpubr)
library(latex2exp)
# library(ungeviz)

load("data/Result.RData")

pred.dfr.2 %>% 
  ggplot(aes(x = as.character(date), color = method)) +
  geom_errorbar(aes(ymin = DeltaC.ll, ymax = DeltaC.ul), position = "dodge") +
  geom_point(aes(y = N_fit), position= position_dodge(width = 0.9)) +
  # ungeviz::geom_hpline(aes(y = DeltaC), size = 1, width = 0.9, linetype = "dashed", color = "black") +
  facet_wrap(~state, scales = "free_y") + theme_light() +
  labs(x = "Date", y = "Daily Cases", color = "Model") +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
