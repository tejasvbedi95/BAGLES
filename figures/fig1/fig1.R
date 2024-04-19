library(tidyverse)
library(ggpubr)
##
load("data/fig1.RData")

gg1 <- ggplot(data = out.growth.temp3, aes(x = Method, y = Value)) +
  geom_jitter(aes(color = Method), width = 0.1, height = 0.1, alpha = 0.7, size = 0.2) +
  geom_boxplot(aes(fill = Method), alpha = 0.3, color = "grey50") +
  facet_grid(phi~ Measure, scales = "free", labeller = label_parsed) +
  theme_light() + labs(title = "With GLC dynamics", x = '', y = '', color = NULL, fill = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(size = 20), legend.position="none") + 
  scale_fill_manual(values = c("red", "orange", "green", "blue", "purple")) +
  scale_color_manual(values = c("red", "orange", "green", "blue", "purple"))

gg2 <- ggplot(data = out.sir.temp3, aes(x = Method, y = Value)) +
  geom_jitter(aes(color = Method), width = 0.1, height = 0.1, alpha = 0.7, size = 0.2) +
  geom_boxplot(aes(fill = Method), alpha = 0.3, color = "grey50") +
  facet_grid(phi~ Measure, scales = "free", labeller = label_parsed) +
  theme_light() + labs(title = "With SIR dynamics", x = '', y = '', color = NULL, fill = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(size = 20), legend.position="bottom") + 
  scale_fill_manual(values = c("red", "orange", "green", "blue", "purple")) + 
  scale_color_manual(values = c("red", "orange", "green", "blue", "purple"))

ggpubr::ggarrange( gg1, gg2, nrow = 2, labels=c("(a)","(b)"), common.legend = T, legend = "bottom")

