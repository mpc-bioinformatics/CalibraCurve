

### Info from Gu at al 2014

# the weighting factor of 1 should be used if a linear relationship between σ and x^0 is observed (σ and σ2 are constants)

# The weighting factor of 1/x should be used if a linear relationship between σ and x0.5 (or between σ2 and x1)
# and a rightward quadratic relationship between σ and x1 are observed.

# Similarly, the weighting factor of 1/x2 should be used if a linear relationship between σ and x1 and
# a concave upward quadratic relationship between σ2 and x are observed.


####################################################################################################

library(tidyverse)

# D: data set
chooseWeightingMethod <- function(D) {

  D_var_sd <- D %>%
    group_by(Concentration) %>%
    summarize(sd = sd(Measurement, na.rm = TRUE), var = var(Measurement, na.rm = TRUE))

  D_var_sd_long <- pivot_longer(D_var_sd, cols = c(sd, var))

  pl <- ggplot(data = D_var_sd_long, mapping = aes(x = Concentration, y = value, colour = name, group = name)) +
    geom_point() +
    geom_line() +
    facet_wrap(vars(name), scales = "free") +
    theme_bw()

  return(pl)
}






####################################################################################################

data("D_MFAP4")

D <- D_MFAP4

library(tidyverse)

D_var_sd <- D_MFAP4 %>%
  group_by(Concentration) %>%
  summarize(sd = sd(Measurement, na.rm = TRUE), var = var(Measurement, na.rm = TRUE))


plot(D_var_sd$Concentration, D_var_sd$var, type = "b")
lines(D_var_sd$Concentration, D_var_sd$sd, col = "red", type = "b")


plot(D_var_sd$Concentration, D_var_sd$var, type = "b", xlim = c(0,5), ylim = c(0, 1e10))

plot(D_var_sd$Concentration, D_var_sd$sd, type = "b")

D_var_sd <- D_var_sd %>%
  mutate(RSD0 = sd, RSD05 = sd/sqrt(Concentration), RSD1 = sd/Concentration)

sd(D_var_sd$RSD0)/mean(D_var_sd$RSD0)*100
sd(D_var_sd$RSD05)/mean(D_var_sd$RSD05)*100
sd(D_var_sd$RSD1)/mean(D_var_sd$RSD1)*100




data("D_ALB")

library(tidyverse)

D_var_sd <- D_ALB %>%
  group_by(Concentration) %>%
  summarize(sd = sd(Measurement, na.rm = TRUE), var = var(Measurement, na.rm = TRUE))


plot(D_var_sd$Concentration, D_var_sd$var, type = "b", ylim = c(0, max(D_var_sd$sd)))
lines(D_var_sd$Concentration, D_var_sd$sd, col = "red", type = "b")

plot(D_var_sd$Concentration, D_var_sd$var, type = "b")

D_var_sd <- D_var_sd %>%
  mutate(RSD0 = sd, RSD05 = sd/sqrt(Concentration), RSD1 = sd/Concentration)

sd(D_var_sd$RSD0)/mean(D_var_sd$RSD0)*100
sd(D_var_sd$RSD05)/mean(D_var_sd$RSD05)*100
sd(D_var_sd$RSD1)/mean(D_var_sd$RSD1)*100




data("D_Apolipoprotein")

library(tidyverse)

D_var_sd <- D_Apolipoprotein %>%
  group_by(Concentration) %>%
  summarize(sd = sd(Measurement, na.rm = TRUE), var = var(Measurement, na.rm = TRUE))


plot(D_var_sd$Concentration, D_var_sd$var, type = "b", ylim = c(0, max(D_var_sd$var)))
lines(D_var_sd$Concentration, D_var_sd$sd, col = "red", type = "b")

plot(D_var_sd$Concentration, D_var_sd$var, type = "b")

D_var_sd <- D_var_sd %>%
  mutate(RSD0 = sd, RSD05 = sd/sqrt(Concentration), RSD1 = sd/Concentration)

sd(D_var_sd$RSD0)/mean(D_var_sd$RSD0)*100
sd(D_var_sd$RSD05)/mean(D_var_sd$RSD05)*100
sd(D_var_sd$RSD1)/mean(D_var_sd$RSD1)*100



x <- c(0.0050, 0.0079, 0.0112, 0.0158, 0.0354, 0.0500, 0.0791, 0.0935, 0.1118)

sd(x)/mean(x)


#88.8
