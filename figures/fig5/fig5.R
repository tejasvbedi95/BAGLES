library(tidyverse)
library(ggpubr)
library(latex2exp)

load("data/Result.RData")

intervals.dfr %>%
  select(-C.ll, -C.ul) %>%
  mutate(state = c(rep("California", 340), rep("New York", 340)),
         N_fit = c(res_rj_cali$N_fit[1:340], res_rj_ny$N_fit[1:340]),
         date = rep(cali.data$date[1:340], 2),
         DeltaC = c(cali.data$DeltaC[1:340], ny.data$DeltaC[1:340])) %>%
  ggplot(aes(x = date, y = DeltaC)) +
  geom_point(size = 0.5) + 
  geom_line(aes(y = N_fit), linetype = "dashed", size = 1, color = "red") +
  #geom_ribbon(aes(ymin = DeltaC.ll, ymax = DeltaC.ul), fill = "red", alpha = 0.3) +
  geom_point(data = pred.dfr, aes(x = date, y = DeltaC), size = 0.5) +
  geom_line(data = pred.dfr, aes(y = N_fit, color = method), linetype = "solid", size = 1) +
  geom_ribbon(data = pred.dfr, aes(ymin = DeltaC.ll, ymax = DeltaC.ul, fill = method), alpha = 0.2) +
  facet_wrap(~state, scales = "free") + labs(y = "Daily Cases") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  theme_light() + labs(x = "", color = "Method", fill = "Method", y = "Daily confirmed case number") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "bottom",
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6),
        text = element_text(size = 20)
  ) 
