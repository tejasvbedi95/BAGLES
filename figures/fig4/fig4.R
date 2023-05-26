library(tidyverse)
library(ggpubr)
library(latex2exp)

load("data/fig4.RData")


g1 <- ggplot(data = C_dfr, aes(x = Days, y = P)) +
  geom_line(linetype = "dotted", color = "grey50", size = 1) +
  
  geom_vline(data = comb_dfr, aes(xintercept = cali.data$date[CI_cali[1,1]]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_dfr, aes(xintercept = cali.data$date[CI_cali[2,1]]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_dfr, aes(xintercept = cali.data$date[CI_cali[4,1]]), color = "orange", linetype = "dashed", size = 1) +

  geom_ribbon(data = comb_dfr, aes(xmin = cali.data$date[CI_cali[1,2]], xmax = cali.data$date[CI_cali[1,3]]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_dfr, aes(xmin = cali.data$date[CI_cali[2,2]], xmax = cali.data$date[CI_cali[2,3]]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_dfr, aes(xmin = cali.data$date[CI_cali[4,2]], xmax = cali.data$date[CI_cali[4,3]]), fill = "orange", alpha = 0.3, size = 1) +

  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  geom_segment(data = comb_dfr, aes(x=Days, xend=Days, y=0, yend=P),color = "black", size = 1) +
  # geom_point(data = comb_dfr, size = 1, aes(color = P_col)) +
  scale_color_manual(values = "cadetblue1") +
  theme_light() + labs(x = "", y = "Posterior probability of inclusion (PPI)", title = "California") + theme(legend.position = "none",
                                                                                                             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 20))
g3 <- ggplot(data = C_ny_dfr, aes(x = Days, y = P)) +
  
  geom_vline(data = comb_ny_dfr, aes(xintercept = ny.data$date[CI_ny[1,1]]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_ny_dfr, aes(xintercept = ny.data$date[CI_ny[2,1]]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_ny_dfr, aes(xintercept = ny.data$date[CI_ny[3,1]]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_ny_dfr, aes(xintercept = ny.data$date[CI_ny[4,1]]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_ny_dfr, aes(xintercept = ny.data$date[CI_ny[5,1]]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_ny_dfr, aes(xintercept = ny.data$date[CI_ny[6,1]]), color = "orange", linetype = "dashed", size = 1) +

  geom_ribbon(data = comb_ny_dfr, aes(xmin = ny.data$date[CI_ny[1,2]], xmax = ny.data$date[CI_ny[1,3]]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_ny_dfr, aes(xmin = ny.data$date[CI_ny[2,2]], xmax = ny.data$date[CI_ny[2,3]]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_ny_dfr, aes(xmin = ny.data$date[CI_ny[3,2]], xmax = ny.data$date[CI_ny[3,3]]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_ny_dfr, aes(xmin = ny.data$date[CI_ny[4,2]], xmax = ny.data$date[CI_ny[4,3]]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_ny_dfr, aes(xmin = ny.data$date[CI_ny[5,2]], xmax = ny.data$date[CI_ny[5,3]]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_ny_dfr, aes(xmin = ny.data$date[CI_ny[6,2]], xmax = ny.data$date[CI_ny[6,3]]), fill = "orange", alpha = 0.3, size = 1) +

  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  geom_line(aes(x = Days, y = P), linetype = "dotted", color = "grey50", size = 1) +
  geom_segment(data = comb_ny_dfr, aes(x=Days, xend=Days, y=0, yend=P),color = "black") +
  # geom_point(data = comb_ny_dfr, size = 1, aes(color = P_col)) +
  scale_color_manual(values = "cadetblue1") +
  theme_light() + labs(x = "", y = "Posterior probability of inclusion (PPI)", title = "New York") + theme(legend.position = "none",
                                                                                                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 20))

ggpubr::ggarrange(g1, g3, ncol = 2, labels = c("(a)", "(b)"))


