library(tidyverse)
library(ggpubr)
library(latex2exp)

load("data/fig2.RData")

Data.labs <- c(Growth = "With GLC dynamics", SIR = "With SIR dynamics")

ggplot(data = C_dfr, aes(x = Days, y = P)) +
  geom_vline(data = comb_dfr %>% filter(Data == "Growth" & Method == "BayesSMEG (auto)"), aes(xintercept = 1), color = "green", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "Growth"), aes(xintercept = 52), color = "green", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "Growth"), aes(xintercept = 103), color = "green", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "SIR"), aes(xintercept = 1), color = "green", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "SIR"), aes(xintercept = 31), color = "green", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "SIR"), aes(xintercept = 61), color = "green", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "SIR"), aes(xintercept = 91), color = "green", size = 1) +
  
  geom_vline(data = comb_dfr %>% filter(Data == "Growth" & Method == "BayesSMEG (auto)"), aes(xintercept = CI[1,1]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "Growth" & Method == "BayesSMEG (auto)"), aes(xintercept = CI[2,1]), color = "orange", linetype = "dashed", size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "Growth" & Method == "BayesSMEG (auto)"), aes(xmin = CI[1,2], xmax = CI[1,3]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "Growth" & Method == "BayesSMEG (auto)"), aes(xmin = CI[2,2], xmax = CI[2,3]), fill = "orange", alpha = 0.3, size = 1) +
  
  geom_vline(data = comb_dfr %>% filter(Data == "Growth" & Method == "BayesSMEG (manual)"), aes(xintercept = CI_3[1,1]), color = "red", linetype = "dashed", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "Growth" & Method == "BayesSMEG (manual)"), aes(xintercept = CI_3[2,1]), color = "red", linetype = "dashed", size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "Growth" & Method == "BayesSMEG (manual)"), aes(xmin = CI_3[1,2], xmax = CI_3[1,3]), fill = "red", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "Growth" & Method == "BayesSMEG (manual)"), aes(xmin = CI_3[2,2], xmax = CI_3[2,3]), fill = "red", alpha = 0.3, size = 1) +
  
  geom_vline(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (auto)"), aes(xintercept = CI_2[1,1]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (auto)"), aes(xintercept = CI_2[2,1]), color = "orange", linetype = "dashed", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (auto)"), aes(xintercept = CI_2[3,1]), color = "orange", linetype = "dashed", size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (auto)"), aes(xmin = CI_2[1,2], xmax = CI_2[1,3]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (auto)"), aes(xmin = CI_2[2,2], xmax = CI_2[2,3]), fill = "orange", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (auto)"), aes(xmin = CI_2[3,2], xmax = CI_2[3,3]), fill = "orange", alpha = 0.3, size = 1) +
  
  geom_vline(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (manual)"), aes(xintercept = CI_4[1,1]), color = "red", linetype = "dashed", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (manual)"), aes(xintercept = CI_4[2,1]), color = "red", linetype = "dashed", size = 1) +
  geom_vline(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (manual)"), aes(xintercept = CI_4[3,1]), color = "red", linetype = "dashed", size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (manual)"), aes(xmin = CI_4[1,2], xmax = CI_4[1,3]), fill = "red", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (manual)"), aes(xmin = CI_4[2,2], xmax = CI_4[2,3]), fill = "red", alpha = 0.3, size = 1) +
  geom_ribbon(data = comb_dfr %>% filter(Data == "SIR" & Method == "BayesSMEG (manual)"), aes(xmin = CI_4[3,2], xmax = CI_4[3,3]), fill = "red", alpha = 0.3, size = 1) +
  geom_line(aes(x = Days, y = P), linetype = "dotted", color = "grey50", size = 1) +
  geom_segment(data = comb_dfr, aes(x=Days, xend=Days, y=0, yend=P),color = "black", size = 1) +
  facet_grid(Method ~ Data, scales = "free",
             labeller = labeller(Data = Data.labs)) +
  # geom_point(data = comb_dfr, size = 1, aes(color = P_col)) +
  scale_color_manual(values = "cadetblue1") +
  theme_light() + labs(x = TeX(r"(Time point $\textit{t}$)"), y = "Posterior probability of inclusion (PPI)") + theme(legend.position = "none", text = element_text(size = 20))

