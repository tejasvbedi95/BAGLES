library(tidyverse)
library(ggpubr)
library(latex2exp)

load("data/fig3.RData")

cols <- c("95%" = "yellow", "80%" = "orange", "50%" = "red")

Model.labs <- c(`SMEG (M = 3)` = "BayesSMEG (manual, M = 3)",
                `SMEG (M = 2)` = "BayesSMEG (manual, M = 2)",
                `GLC` = "GLC (M = 1)")

K.dfr %>%
  ggplot(aes(x = X, y = K, group = Model)) +
  geom_line(alpha = 0.5) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K1_CI_1[1], X < K1_CI_1[2], Model == "SMEG (M = 3)"), 
              aes(ymax = K, fill = "95%"), ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K1_CI_2[1], X < K1_CI_2[2], Model == "SMEG (M = 3)"), 
              aes(ymax = K, fill = "80%"), ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K1_CI_3[1], X < K1_CI_3[2], Model == "SMEG (M = 3)"), 
              aes(ymax = K, fill = "50%"), ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K2_CI_1[1], X < K2_CI_1[2], Model == "SMEG (M = 2)"), 
              aes(ymax = K, fill = "95%"), ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K2_CI_2[1], X < K2_CI_2[2], Model == "SMEG (M = 2)"), 
              aes(ymax = K, fill = "80%"), ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K2_CI_3[1], X < K2_CI_3[2], Model == "SMEG (M = 2)"), 
              aes(ymax = K, fill = "50%"), ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K3_CI_1[1], X < K3_CI_1[2], Model == "GLC"), 
              aes(ymax = K, fill = "95%"), ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K3_CI_2[1], X < K3_CI_2[2], Model == "GLC"), 
              aes(ymax = K, fill = "80%"), ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K3_CI_3[1], X < K3_CI_3[2], Model == "GLC"), 
              aes(ymax = K, fill = "50%"), ymin = 0) +
  geom_vline(xintercept = 15000, color = "green") +
  facet_wrap(~Model, nrow = 3, labeller = labeller(Model = Model.labs)) + theme_light() +
  labs(x = TeX(r"(Final Epidemic Size $\textit{K}$)"), y = "Density", fill = "Credible intervals") +
  theme(legend.position = "bottom", text = element_text(size = 20)) +
  scale_fill_manual(values = cols)
