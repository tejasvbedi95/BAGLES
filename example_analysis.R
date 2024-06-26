# Loading libraries

library(tidyverse)
library(Rcpp)

# Loading functions

sourceCpp('scripts/core_fixed.cpp')
sourceCpp('scripts/core_auto.cpp')
sourceCpp('scripts/Delta_arrange.cpp')
source('scripts/functions.R')

### Loading data

pop.data <- read.csv("data/pop.data.csv",header=T,row.names=1)
us.data <- read.csv("data/us-states.csv", header = T)


# Simulating growth epidemiological data

T = 150
M = 3
lambda <- c(0.1, 0.06, 0.08)
p <- c(0.9, 0.85, 0.9)
alpha <- c(1, 1, 1)
K <- c(10000, 9000, 15000)
phi <- c(10, 10, 10)
delta <- delta_arrange(T, M)
N <- 200000
rho <- 0.3

sim <- cp_simulator( T, M, phi, lambda, p, alpha, K, delta, C0 = 100, seed = 7)

## BayesSMEG with fixed M = 3
res_fixed <- growth_cp(sim$C, M = 3, is_p = -1.0, is_alpha = 1.0, POP = ceiling(N*rho), T_fin = 0, w = c(0.2, 0.4, 0.4), store = T)
res_mat_fixed <- apply(res_fixed$CI_mat[-1,], 2, mean)
CI_fixed <- CI_cp(res_fixed, which(res_fixed$map$delta_map == 1))

## BayesSMEG with unknown M (alpha = 0.001)
res_auto <- growth_cp_rj(sim$DeltaC, M_max = 20, POP = ceiling(N*rho), T_fin = 0, alpha = 0.001, store = T)
res_mat_auto <- apply(res_auto$CI_mat[-1,], 2, mean)
CI_auto <- CI_cp(res_auto, which(res_auto$map$delta_map == 1))

### Change point detection plots with PPI and cumulative cases mapped together

comb_dfr <- data.frame(Days = rep(1:150, 2),
                       P = c(res_mat_fixed, res_mat_auto),
                       Method = c(rep("fixed", 150), rep("auto", 150)))

comb_dfr <- comb_dfr %>%
  mutate(Measure = as.factor(Method)) %>%
  mutate(P_col = ifelse(P > 0, "red", NA))

C_dfr <- data.frame(Days = rep(1:150, 2),
                    P = rep(c(sim$C/max(sim$C)), 2))

ggplot(data = C_dfr, aes(x = Days, y = P)) +
  
  geom_vline(data = comb_dfr, aes(xintercept = 52), color = "green") +
  geom_vline(data = comb_dfr, aes(xintercept = 103), color = "green") +
  
  geom_vline(data = comb_dfr %>% filter(Method == "auto"), aes(xintercept = CI_auto[1,1]), color = "orange", linetype = "dashed") +
  geom_vline(data = comb_dfr %>% filter(Method == "auto"), aes(xintercept = CI_auto[2,1]), color = "orange", linetype = "dashed") +
  geom_ribbon(data = comb_dfr %>% filter(Method == "auto"), aes(xmin = CI_auto[1,2], xmax = CI_auto[1,3]), fill = "orange", alpha = 0.3) +
  geom_ribbon(data = comb_dfr %>% filter(Method == "auto"), aes(xmin = CI_auto[2,2], xmax = CI_auto[2,3]), fill = "orange", alpha = 0.3) +
  
  geom_vline(data = comb_dfr %>% filter(Method == "fixed"), aes(xintercept = CI_fixed[1,1]), color = "red", linetype = "dashed") +
  geom_vline(data = comb_dfr %>% filter(Method == "fixed"), aes(xintercept = CI_fixed[2,1]), color = "red", linetype = "dashed") +
  geom_ribbon(data = comb_dfr %>% filter(Method == "fixed"), aes(xmin = CI_fixed[1,2], xmax = CI_fixed[1,3]), fill = "red", alpha = 0.3) +
  geom_ribbon(data = comb_dfr %>% filter(Method == "fixed"), aes(xmin = CI_fixed[2,2], xmax = CI_fixed[2,3]), fill = "red", alpha = 0.3) +
  
  geom_line(aes(x = Days, y = P), linetype = "dotted", color = "grey50") +
  geom_segment(data = comb_dfr, aes(x=Days, xend=Days, y=0, yend=P),color = "black") +
  facet_wrap(~ Method, scales = "free") +
  geom_point(data = comb_dfr, size = 0.5, aes(color = P_col)) +
  theme_light() + labs(y = "PPI") + theme(legend.position = "none")


### Estimation of final epidemic size with density and intervals

res_fixed <- growth_cp(sim$C, M = 3, is_p = -1.0, is_alpha = 1.0, POP = ceiling(N*rho), T_fin = 0, w = c(0.2, 0.4, 0.4), store = T)
#res_fixed <- growth_cp(sim$C, M = 1, is_p = -1.0, is_alpha = 1.0, POP = ceiling(N*rho), T_fin = 0, w = c(0.2, 0.4, 0.4), store = T)


K_dist <- res_fixed$store_list$K_store[3,]


K.dfr <- data.frame(X = density(K_dist)$x,
                    K = density(K_dist)$y)

K_CI_1 <- quantile(K_dist, probs = c(0.025, 0.975))
K_CI_2 <- quantile(K_dist, probs = c(0.1, 0.9))
K_CI_3 <- quantile(K_dist, probs = c(0.25, 0.75))


cols <- c("0.95" = "yellow", "0.80" = "orange", "0.50" = "red")


K.dfr %>%
  ggplot(aes(x = X, y = K)) +
  geom_line(alpha = 0.5) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K_CI_1[1], X < K_CI_1[2]), 
              aes(ymax = K), fill = "yellow", ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K_CI_2[1], X < K_CI_2[2]), 
              aes(ymax = K), fill = "orange", ymin = 0) +
  geom_ribbon(data = K.dfr %>% 
                filter(X > K_CI_3[1], X < K_CI_3[2]), 
              aes(ymax = K), fill = "red", ymin = 0) +
  geom_vline(xintercept = 15000, color = "green") +
  theme_light() + labs(x = "K", y = "Density") +
  scale_fill_manual(values = cols, aesthetics = c("fill"), name = "HPD Intervals")


## Long term prediction of New York daily cases via BayesSMEG (auto)

ny.data <- us.data %>%
  select(-fips, -deaths) %>%
  filter(state == "New York", cases > 100) %>%
  rename(C = cases) %>%
  mutate(DeltaC = c(NA, diff(C)), date = as.Date(date))

N_ny <- pop.data[which(rownames(pop.data) == "New York"),]

res_rj_ny <- growth_cp_rj(ny.data$DeltaC[-1][1:340], M_max = 50,
                            POP = N_ny*0.3, T_fin = 150, alpha = 0.000001, store = T)

pred.dfr <- data.frame(DeltaC = ny.data$DeltaC[341:490],
                       N_fit = res_rj_ny$N_pred_mean,
                       DeltaC.ll = res_rj_ny$N_pred_lwr,
                       DeltaC.ul = res_rj_ny$N_pred_upp,
                       date = ny.data$date[341:490])

## For better scaling of plots, the intervals for GLC and SIR are bounded
pred.dfr[pred.dfr$state == "New York", "DeltaC.ul"] <- ifelse(pred.dfr[pred.dfr$state == "New York", "DeltaC.ul"] > 25000, 25000,
                                                                pred.dfr[pred.dfr$state == "New York", "DeltaC.ul"])


data.frame(N_fit = res_rj_ny$N_fit[1:340],
         date = ny.data$date[1:340],
         DeltaC = ny.data$DeltaC[1:340]) %>%
  ggplot(aes(x = date, y = DeltaC)) +
  geom_point(size = 0.5) + 
  geom_line(aes(y = N_fit), linetype = "dashed", size = 1, color = "red") +
  #geom_ribbon(aes(ymin = DeltaC.ll, ymax = DeltaC.ul), fill = "red", alpha = 0.3) +
  geom_point(data = pred.dfr, aes(x = date, y = DeltaC), size = 0.5) +
  geom_line(data = pred.dfr, aes(y = N_fit), linetype = "solid", size = 1, color = "red") +
  geom_ribbon(data = pred.dfr, aes(ymin = DeltaC.ll, ymax = DeltaC.ul), alpha = 0.2, fill = "red") +
  labs(y = "New York Daily Cases") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  theme_light() + labs(x = "Months", color = "Model", fill = "Model") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(.55, .95),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(6, 6, 6, 6)
  ) 








