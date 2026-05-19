############ ClineDescription - Shell Color Across Ecology ###############

# BASIC INFO ----

# Researcher: [AuthorName] (main analyses) & [AuthorName2] (SLiM simulations)
# Section aim: Are shell colors associated with specific ecosystem types? If so, can non-selective mechanisms explain the association?
# Last update: 2026-05-19

# SETUP SECTION ----

## 1.- Packages ----
library(tidyverse) # Basic plots and statistical analyses
library(vegan) #  Bray-Curtis
library(reshape2) # For formats in phenotypic distances
library(geosphere) # For calculating geographical distances
library(glmmTMB) # Beta-binomial regression
library(emmeans) # Extract confidence intervals for predicted line in beta-binomial models
library(DHARMa) # Model diagnostics
library(ecodist) # Mantel, MRM
library(ape) # Mantel
library(cowplot) # Composite figures
library(grid) # Composite figures 

## 2.- Dataset loading ----
setwd("") # Set the directory of the dataset file
read.csv("CD_DA_RD_PopulationLevel.csv")  
read.csv("CD_DA_RD_TemporalSeries.csv") 
readr::read_tsv("CD_DA_RD_Simulations_ReplicatesDefault.tsv")
readr::read_tsv("CD_DA_RD_Simulations_ConsecutiveGenerations.tsv", show_col_types = FALSE)

## 3.- Data handling and overview ----

# PopulationLevel dataset:
PopLev <- read.csv("CD_DA_RD_PopulationLevel.csv")  %>%
  mutate(
    Ria = factor(Ria, levels = c("MurosNoia", "Arousa","Pontevedra", "Vigo")), # Set an order from north to south
    Coast = factor(Coast, levels = c("South", "North")),
    Cluster = factor(Cluster, levels = c("WaveSheltered", "WaveIntermediate", "WaveExposed"))) # Convert variables into adequate formats
str(PopLev) # Check the structure of the data

# TemporalSeries dataset:
TempSer <- read.csv("CD_DA_RD_TemporalSeries.csv") %>% 
  mutate(Year = factor(Year, levels = c("1979", "2022"))) # Convert the variable into a factor
str(TempSer) # Check the structure of the data

# ReplicatesDefault dataset (for simulations):
RepDef <- readr::read_tsv("CD_DA_RD_Simulations_ReplicatesDefault.tsv")
str(RepDef) # Check the structure of the data

# ConsecutiveGenerations dataset (for simulations):
ConsGen <- readr::read_tsv("CD_DA_RD_Simulations_ConsecutiveGenerations.tsv", show_col_types = FALSE)
str(ConsGen) # Check the structure of the data

# Expected frequencies dataset (for simulations):
Expected <- readr::read_tsv("CD_DA_RD_Simulations_ExpectedDefault.tsv", show_col_types = FALSE)
str(Expected) # Check the structure of the data

# Alternative scenarios datasets (drift and migration) (for simulations):
Drift <- readr::read_tsv("CD_DA_RD_Simulations_ReplicatesDrift.tsv", show_col_types = FALSE) %>% 
  dplyr::mutate(scenario = "drift")

Migration <- readr::read_tsv("CD_DA_RD_Simulations_ReplicatesMigration.tsv", show_col_types = FALSE) %>% 
  dplyr::mutate(scenario = "migration")

AltScen <- dplyr::bind_rows(Drift, Migration)
str(AltScen) # Check the structure of the data

## 4.- Themes for figures ----

# Theme for figures, v.1:
theme_cd <- function(BaseSize = 24, BaseFamily = "sans") {
  theme_classic(base_size = BaseSize, base_family = BaseFamily) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = BaseSize + 2, colour = "black"),
      axis.text = element_text(size = BaseSize, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks.length = unit(2.2, "mm"),
      legend.title = element_blank(),
      legend.text = element_text(size = BaseSize),
      legend.key.size = unit(4, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = BaseSize + 1, colour = "black"),
      plot.margin = margin(3, 3, 3, 3, unit = "mm")
    )
}

# Theme for tags:
tags_cd <- function(p, tag, size = 22, x = 0, y = 1) {
  p +
    labs(tag = tag) +
    theme(
      plot.tag = element_text(size = size, face = "bold"),
      plot.tag.position = c(x, y)
    )
}

# ANALYSIS SECTION ----

## Main analyses ----

### 1.- Morph frequencies across ecological profiles ----

# Subsection aim: Does color vary with EcolPC1? How can we model that relationship?

# ///


#### 1.1.- Pie charts in map figure ----

# NOTE: Once the pie charts and color gradients were obtained, the corresponding map figure was created in InkScape

###### 1.3.1.- Pie charts of Lineata, Fulva, and the rest ----

# Select the variables that we will use:
PopLev_PieCharts_limited <- PopLev %>% 
  dplyr::select(Ria, LocNumCumulative, LocName, Fulva_Frequency, 
                Lineata_Frequency, 
                Rest_Frequency
  )

# Set the variables into a long format:
PopLev_PieCharts_limited_long <- PopLev_PieCharts_limited %>% 
  pivot_longer(
    cols = ends_with("_Frequency"),
    names_to = "Morph",
    values_to = "Frequency"
  )

# Plot with aesthetics:
PieLFR <- ggplot(PopLev_PieCharts_limited_long, aes(x = "", y = Frequency, fill = Morph)) +
  geom_bar(stat = "identity", width = 1, color = "black",linewidth = 0.2) +
  coord_polar(theta = "y") +
  facet_wrap(~ LocNumCumulative) +
  scale_fill_manual(
    values = c(
      "Fulva_Frequency" = "#D6D2A9",
      "Lineata_Frequency" = "#A09D97",
      "Rest_Frequency" = "#95B877"
    )
  ) +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold")
  )

PieLFR # View the pie charts

###### 1.3.2.- Color palette of EcolPC1 ----

# Create a new variable to identify localities precisely (so that errors are prevented):
PopLev <- PopLev %>% 
  mutate(Ria_Coast_Num = paste(Ria, Coast, LocNumRiaCoast, sep = "_"))

# Select the variables:
EcolPC1_Colors <- PopLev %>% 
  dplyr::select(LocName, LocNumCumulative, Ria_Coast_Num, EcolPC1)

# Check the range of values:
range(PopLev$EcolPC1, na.rm = TRUE)

# Plot with aesthetics:
EcolPC1_Palette <- ggplot(EcolPC1_Colors, aes(x = LocNumCumulative, y = reorder(LocNumCumulative, EcolPC1), color = EcolPC1)) +
  geom_point(size = 6) +
  geom_text(aes(label = Ria_Coast_Num), color = "black",
            size = 3.5) +
  #scale_color_viridis(option = "plasma", direction = 1) +
  scale_color_gradient(
    low = "purple",
    high = "yellow",
  ) +
  theme_void() +
  theme(
    legend.key.height = unit(0.8, "cm"),
    legend.title = element_blank(),
    legend.ticks = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

EcolPC1_Palette # View the palette

#### 1.2.- Fulva snails ~ EcolPC1 ----

# Run a beta-binomial model:
Fulva_Beta <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                           EcolPC1,
                         family = betabinomial(link = "logit"),
                         data = PopLev)

# Run diagnostics: 
SimulatedResiduals_Fulva_Beta <- simulateResiduals(Fulva_Beta)
plot(SimulatedResiduals_Fulva_Beta)
# OUTPUT: The fit looks reasonable

# See the results table: 
summary(Fulva_Beta)

# Model to test whether the pattern differs between rias:
Fulva_Beta_Rias <- glmmTMB(cbind(Fulva_Counts, n-Fulva_Counts) ~ EcolPC1 * Ria, 
                          family = betabinomial(link = "logit"),
                          data = PopLev) 

# Run diagnostics: 
SimulatedResiduals_Fulva_Beta_Rias <- simulateResiduals(Fulva_Beta_Rias)
plot(SimulatedResiduals_Fulva_Beta_Rias)
# OUTPUT: The fit looks good but DHARMa cannot guarantee it is checking the model fit correctly. 

# See the summary table:
summary(Fulva_Beta_Rias) # Inspect predictors to see which interactions are different

# Likelihood ratio test:
Fulva_Beta_Reduced <- glmmTMB(cbind(Fulva_Counts, n-Fulva_Counts) ~ EcolPC1 + Ria, 
                         family = betabinomial(link = "logit"),
                         data = PopLev) # Model without interaction

anova(Fulva_Beta_Rias, Fulva_Beta_Reduced) # Formally assess the differences between models
# OUTPUT: The pattern (slope) differs between rias.

###### 1.2.1.- Fulva in Muros-Noia ----

# Fix the "MurosNoia" name by adding a "-" between them:
PopLev <- PopLev %>% 
  mutate(Ria = fct_recode(Ria, 
                          "Muros-Noia" = "MurosNoia"))

# Filter the Muros-Noia ria:
PopLev_MN <- PopLev %>% 
  filter(Ria %in% c("Muros-Noia"))

# Run the model for Muros-Noia:
Fulva_Beta_MN <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                        EcolPC1,
                      family = betabinomial(link = "logit"),
                      data = PopLev_MN)

# Run diagnostics: 
SimulatedResiduals_Fulva_Beta_MN <- simulateResiduals(Fulva_Beta_MN)
plot(SimulatedResiduals_Fulva_Beta_MN)
# OUTPUT: The fit looks reasonable, but not without problems

# See the results table: 
summary(Fulva_Beta_MN)

# Prediction grid: 
NewData_Fulva_MN <- data.frame(
  EcolPC1 = seq(
    min(PopLev_MN$EcolPC1, na.rm = TRUE),
    max(PopLev_MN$EcolPC1, na.rm = TRUE),
    length.out = 200))

# Predicted probabilities for regression line:
Prediction_Fulva_MN <- predict(Fulva_Beta_MN,
                         newdata = NewData_Fulva_MN,
                         type = "response",
                         se.fit = TRUE) 

NewData_Fulva_MN$Fit <- Prediction_Fulva_MN$fit
NewData_Fulva_MN$SE <- Prediction_Fulva_MN$se.fit

# Confidence intervals:
NewData_Fulva_MN <- NewData_Fulva_MN |>
  mutate(
    Lower = pmax(0, Fit - 1.96 * SE),
    Upper = pmin(1, Fit + 1.96 * SE))

# Plot with aesthetics:
Fig_4A_1 <- ggplot(PopLev_MN, aes(x = EcolPC1, y = Fulva_Frequency)) +
  geom_point(size = 4, shape = 21, fill = "#D6D2A9", color = "black") +
  geom_ribbon(data = NewData_Fulva_MN, aes(x = EcolPC1, ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = NewData_Fulva_MN, aes (y = Fit), linewidth = 1.5, color = "black") +
  labs(title = "Muros-Noia",
      x = "Ecological PC1",
      y = expression(paste("Freq. of ", italic("Fulva")))) +
  theme_cd(18) + 
  theme(plot.title = element_text(color = "black", size = 18)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  #scale_x_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Same for the x axis
  scale_x_reverse(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_4A_1 # See the plot

###### 1.2.2.- Fulva in Arousa ----

# Filter the Arousa ria:
PopLev_A <- PopLev %>% 
  filter(Ria %in% c("Arousa"))

# Run the model for Arousa:
Fulva_Beta_A <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                           EcolPC1,
                         family = betabinomial(link = "logit"),
                         data = PopLev_A)

# Run diagnostics: 
SimulatedResiduals_Fulva_Beta_A <- simulateResiduals(Fulva_Beta_A)
plot(SimulatedResiduals_Fulva_Beta_A)
# OUTPUT: The fit looks good

# See the results table: 
summary(Fulva_Beta_A)

# Prediction grid: 
NewData_Fulva_A <- data.frame(
  EcolPC1 = seq(
    min(PopLev_A$EcolPC1, na.rm = TRUE),
    max(PopLev_A$EcolPC1, na.rm = TRUE),
    length.out = 200))

# Predicted probabilities for regression line:
Prediction_Fulva_A <- predict(Fulva_Beta_A,
                         newdata = NewData_Fulva_A,
                         type = "response",
                         se.fit = TRUE) 

NewData_Fulva_A$Fit <- Prediction_Fulva_A$fit
NewData_Fulva_A$SE <- Prediction_Fulva_A$se.fit

# Confidence intervals:
NewData_Fulva_A <- NewData_Fulva_A |>
  mutate(
    Lower = pmax(0, Fit - 1.96 * SE),
    Upper = pmin(1, Fit + 1.96 * SE))

# Plot with aesthetics:
Fig_4A_2 <- ggplot(PopLev_A, aes(x = EcolPC1, y = Fulva_Frequency)) +
  geom_point(size = 4, shape = 21, fill = "#D6D2A9", color = "black") +
  geom_ribbon(data = NewData_Fulva_A, aes(x = EcolPC1, ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = NewData_Fulva_A, aes (y = Fit), linewidth = 1.5, color = "black") +
  labs(title = "Arousa",
       x = "Ecological PC1",
       y = expression(paste("Freq. of ", italic("Fulva")))) +
  theme_cd(18) + 
  theme(plot.title = element_text(color = "black", size = 18)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  #scale_x_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Same for the x axis
  scale_x_reverse(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_4A_2 # See the plot

###### 1.2.3.- Fulva in Pontevedra ----

# Filter the Pontevedra ria:
PopLev_P <- PopLev %>% 
  filter(Ria %in% c("Pontevedra"))

# Run the model for Pontevedra:
Fulva_Beta_P <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                          EcolPC1,
                        family = betabinomial(link = "logit"),
                        data = PopLev_P)

# Run diagnostics: 
SimulatedResiduals_Fulva_Beta_P <- simulateResiduals(Fulva_Beta_P)
plot(SimulatedResiduals_Fulva_Beta_P)
# OUTPUT: The fit doesn't look particularly good for the residuals

# See the results table: 
summary(Fulva_Beta_P)
Fulva_Beta_P$fit$convergence # The model converged properly

# Prediction grid: 
NewData_Fulva_P <- data.frame(
  EcolPC1 = seq(
    min(PopLev_P$EcolPC1, na.rm = TRUE),
    max(PopLev_P$EcolPC1, na.rm = TRUE),
    length.out = 200))

# Predicted probabilities for regression line:
Prediction_Fulva_P <- predict(Fulva_Beta_P,
                        newdata = NewData_Fulva_P,
                        type = "response",
                        se.fit = TRUE) 

NewData_Fulva_P$Fit <- Prediction_Fulva_P$fit
NewData_Fulva_P$SE <- Prediction_Fulva_P$se.fit

# Confidence intervals:
NewData_Fulva_P <- NewData_Fulva_P |>
  mutate(
    Lower = pmax(0, Fit - 1.96 * SE),
    Upper = pmin(1, Fit + 1.96 * SE))

# Plot with aesthetics:
Fig_4A_3 <- ggplot(PopLev_P, aes(x = EcolPC1, y = Fulva_Frequency)) +
  geom_point(size = 4, shape = 21, fill = "#D6D2A9", color = "black") +
  geom_ribbon(data = NewData_Fulva_P, aes(x = EcolPC1, ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = NewData_Fulva_P, aes (y = Fit), linewidth = 1.5, color = "black") +
  labs(title = "Pontevedra",
       x = "Ecological PC1",
       y = expression(paste("Freq. of ", italic("Fulva")))) +
  theme_cd(18) + 
  theme(plot.title = element_text(color = "black", size = 18)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  #scale_x_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Same for the x axis
  scale_x_reverse(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_4A_3 # See the plot

###### 1.2.4.- Fulva in Vigo ----

# Filter the Vigo ria:
PopLev_V <- PopLev %>% 
  filter(Ria %in% c("Vigo"))

# Run the model for Vigo:
Fulva_Beta_V <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                          EcolPC1,
                        family = betabinomial(link = "logit"),
                        data = PopLev_V)

# Run diagnostics: 
SimulatedResiduals_Fulva_Beta_V <- simulateResiduals(Fulva_Beta_V)
plot(SimulatedResiduals_Fulva_Beta_V)
# OUTPUT: The fit looks good

# See the results table: 
summary(Fulva_Beta_V)

# Prediction grid: 
NewData_Fulva_V <- data.frame(
  EcolPC1 = seq(
    min(PopLev_V$EcolPC1, na.rm = TRUE),
    max(PopLev_V$EcolPC1, na.rm = TRUE),
    length.out = 200))

# Predicted probabilities for regression line:
Prediction_Fulva_V <- predict(Fulva_Beta_V,
                        newdata = NewData_Fulva_V,
                        type = "response",
                        se.fit = TRUE) 

NewData_Fulva_V$Fit <- Prediction_Fulva_V$fit
NewData_Fulva_V$SE <- Prediction_Fulva_V$se.fit

# Confidence intervals:
NewData_Fulva_V <- NewData_Fulva_V |>
  mutate(
    Lower = pmax(0, Fit - 1.96 * SE),
    Upper = pmin(1, Fit + 1.96 * SE))

# Plot with aesthetics:
Fig_4A_4 <- ggplot(PopLev_V, aes(x = EcolPC1, y = Fulva_Frequency)) +
  geom_point(size = 4, shape = 21, fill = "#D6D2A9", color = "black") +
  geom_ribbon(data = NewData_Fulva_V, aes(x = EcolPC1, ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = NewData_Fulva_V, aes (y = Fit), linewidth = 1.5, color = "black") +
  labs(title = "Vigo",
       x = "Ecological PC1",
       y = expression(paste("Freq. of ", italic("Fulva")))) +
  theme_cd(18) + 
  theme(plot.title = element_text(color = "black", size = 18)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  #scale_x_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Same for the x axis
  scale_x_reverse(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_4A_4 # See the plot

#### 1.3.- Lineata snails ~ EcolPC1 ----

# Run a beta-binomial model:
Lineata_Beta <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                             EcolPC1,
                           family = betabinomial(link = "logit"),
                           data = PopLev)

# Run diagnostics: 
SimulatedResiduals_Lineata_Beta <- simulateResiduals(Lineata_Beta)
plot(SimulatedResiduals_Lineata_Beta)
# OUTPUT: Dispersion test is not good but the fit overall looks good

# See the results table: 
summary(Lineata_Beta)

# Model to test whether the pattern differs between rias:
Lineata_Beta_Rias <- glmmTMB(cbind(Lineata_Counts, n-Lineata_Counts) ~ EcolPC1 * Ria, 
                            family = betabinomial(link = "logit"),
                            data = PopLev) 

# Run diagnostics: 
SimulatedResiduals_Lineata_Beta_Rias <- simulateResiduals(Lineata_Beta_Rias)
plot(SimulatedResiduals_Lineata_Beta_Rias)
# OUTPUT: The fit is good 

# See the summary table:
summary(Lineata_Beta_Rias) # Inspect predictors to see which interactions are different

# Likelihood ratio test:
Lineata_Beta_Reduced <- glmmTMB(cbind(Lineata_Counts, n-Lineata_Counts) ~ EcolPC1 + Ria, 
                           family = betabinomial(link = "logit"),
                           data = PopLev) # Model without interaction

anova(Lineata_Beta_Rias, Lineata_Beta_Reduced) # Formally assess the differences between models
# OUTPUT: The pattern (slope) does not differ between rias.

###### 1.3.1.- Lineata in Muros-Noia ----

# Fix the "MurosNoia" name by adding a "-" between them:
PopLev <- PopLev %>% 
  mutate(Ria = fct_recode(Ria, 
                          "Muros-Noia" = "MurosNoia"))

# Filter the Muros-Noia ria:
PopLev_MN <- PopLev %>% 
  filter(Ria %in% c("Muros-Noia"))

# Run the model for Muros-Noia:
Lineata_Beta_MN <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                           EcolPC1,
                         family = betabinomial(link = "logit"),
                         data = PopLev_MN)

# Run diagnostics: 
SimulatedResiduals_Lineata_Beta_MN <- simulateResiduals(Lineata_Beta_MN)
plot(SimulatedResiduals_Lineata_Beta_MN)
# OUTPUT: The fit looks reasonable, but there is signifivant dispersion (although not worrisome upon inspection)

# See the results table: 
summary(Lineata_Beta_MN)

# Prediction grid: 
NewData_Lineata_MN <- data.frame(
  EcolPC1 = seq(
    min(PopLev_MN$EcolPC1, na.rm = TRUE),
    max(PopLev_MN$EcolPC1, na.rm = TRUE),
    length.out = 200))

# Predicted probabilities for regression line:
Prediction_Lineata_MN <- predict(Lineata_Beta_MN,
                         newdata = NewData_Lineata_MN,
                         type = "response",
                         se.fit = TRUE) 

NewData_Lineata_MN$Fit <- Prediction_Lineata_MN$fit
NewData_Lineata_MN$SE <- Prediction_Lineata_MN$se.fit

# Confidence intervals:
NewData_Lineata_MN <- NewData_Lineata_MN |>
  mutate(
    Lower = pmax(0, Fit - 1.96 * SE),
    Upper = pmin(1, Fit + 1.96 * SE))

# Plot with aesthetics:
Fig_4A_5 <- ggplot(PopLev_MN, aes(x = EcolPC1, y = Lineata_Frequency)) +
  geom_point(size = 4, shape = 21, fill = "#A09D97", color = "black") +
  geom_ribbon(data = NewData_Lineata_MN, aes(x = EcolPC1, ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = NewData_Lineata_MN, aes (y = Fit), linewidth = 1.5, color = "black") +
  labs(title = "Muros-Noia",
       x = "Ecological PC1",
       y = expression(paste("Freq. of ", italic("Lineata")))) +
  theme_cd(18) + 
  theme(plot.title = element_text(color = "black", size = 18)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  #scale_x_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Same for the x axis
  scale_x_reverse(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_4A_5 # See the plot

###### 1.3.2.- Lineata in Arousa ----

# Filter the Arousa ria:
PopLev_A <- PopLev %>% 
  filter(Ria %in% c("Arousa"))

# Run the model for Arousa:
Lineata_Beta_A <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                             EcolPC1,
                           family = betabinomial(link = "logit"),
                           data = PopLev_A)

# Run diagnostics: 
SimulatedResiduals_Lineata_Beta_A <- simulateResiduals(Lineata_Beta_A)
plot(SimulatedResiduals_Lineata_Beta_A)
# OUTPUT: The fit looks reasonable

# See the results table: 
summary(Lineata_Beta_A)

# Prediction grid: 
NewData_Lineata_A <- data.frame(
  EcolPC1 = seq(
    min(PopLev_A$EcolPC1, na.rm = TRUE),
    max(PopLev_A$EcolPC1, na.rm = TRUE),
    length.out = 200))

# Predicted probabilities for regression line:
Prediction_Lineata_A <- predict(Lineata_Beta_A,
                                 newdata = NewData_Lineata_A,
                                 type = "response",
                                 se.fit = TRUE) 

NewData_Lineata_A$Fit <- Prediction_Lineata_A$fit
NewData_Lineata_A$SE <- Prediction_Lineata_A$se.fit

# Confidence intervals:
NewData_Lineata_A <- NewData_Lineata_A |>
  mutate(
    Lower = pmax(0, Fit - 1.96 * SE),
    Upper = pmin(1, Fit + 1.96 * SE))

# Plot with aesthetics:
Fig_4A_6 <- ggplot(PopLev_A, aes(x = EcolPC1, y = Lineata_Frequency)) +
  geom_point(size = 4, shape = 21, fill = "#A09D97", color = "black") +
  geom_ribbon(data = NewData_Lineata_A, aes(x = EcolPC1, ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = NewData_Lineata_A, aes (y = Fit), linewidth = 1.5, color = "black") +
  labs(title = "Arousa",
       x = "Ecological PC1",
       y = expression(paste("Freq. of ", italic("Lineata")))) +
  theme_cd(18) + 
  theme(plot.title = element_text(color = "black", size = 18)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  #scale_x_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Same for the x axis
  scale_x_reverse(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_4A_6 # See the plot

###### 1.3.3.- Lineata in Pontevedra ----

# Filter the Pontevedra ria:
PopLev_P <- PopLev %>% 
  filter(Ria %in% c("Pontevedra"))

# Run the model for Pontevedra:
Lineata_Beta_P <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                            EcolPC1,
                          family = betabinomial(link = "logit"),
                          data = PopLev_P)

# Run diagnostics: 
SimulatedResiduals_Lineata_Beta_P <- simulateResiduals(Lineata_Beta_P)
plot(SimulatedResiduals_Lineata_Beta_P)
# OUTPUT: The fit looks good

# See the results table: 
summary(Lineata_Beta_P)

# Prediction grid: 
NewData_Lineata_P <- data.frame(
  EcolPC1 = seq(
    min(PopLev_P$EcolPC1, na.rm = TRUE),
    max(PopLev_P$EcolPC1, na.rm = TRUE),
    length.out = 200))

# Predicted probabilities for regression line:
Prediction_Lineata_P <- predict(Lineata_Beta_P,
                                newdata = NewData_Lineata_P,
                                type = "response",
                                se.fit = TRUE) 

NewData_Lineata_P$Fit <- Prediction_Lineata_P$fit
NewData_Lineata_P$SE <- Prediction_Lineata_P$se.fit

# Confidence intervals:
NewData_Lineata_P <- NewData_Lineata_P |>
  mutate(
    Lower = pmax(0, Fit - 1.96 * SE),
    Upper = pmin(1, Fit + 1.96 * SE))

# Plot with aesthetics:
Fig_4A_7 <- ggplot(PopLev_P, aes(x = EcolPC1, y = Lineata_Frequency)) +
  geom_point(size = 4, shape = 21, fill = "#A09D97", color = "black") +
  geom_ribbon(data = NewData_Lineata_P, aes(x = EcolPC1, ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = NewData_Lineata_P, aes (y = Fit), linewidth = 1.5, color = "black") +
  labs(title = "Pontevedra",
       x = "Ecological PC1",
       y = expression(paste("Freq. of ", italic("Lineata")))) +
  theme_cd(18) + 
  theme(plot.title = element_text(color = "black", size = 18)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  #scale_x_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Same for the x axis
  scale_x_reverse(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_4A_7 # See the plot

###### 1.3.4.- Lineata in Vigo ----

# Filter the Vigo ria:
PopLev_V <- PopLev %>% 
  filter(Ria %in% c("Vigo"))

# Run the model for Vigo:
Lineata_Beta_V <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                            EcolPC1,
                          family = betabinomial(link = "logit"),
                          data = PopLev_V)

# Run diagnostics: 
SimulatedResiduals_Lineata_Beta_V <- simulateResiduals(Lineata_Beta_V)
plot(SimulatedResiduals_Lineata_Beta_V)
# OUTPUT: The fit looks good, although with a small problem

# See the results table: 
summary(Lineata_Beta_V)

# Prediction grid: 
NewData_Lineata_V <- data.frame(
  EcolPC1 = seq(
    min(PopLev_V$EcolPC1, na.rm = TRUE),
    max(PopLev_V$EcolPC1, na.rm = TRUE),
    length.out = 200))

# Predicted probabilities for regression line:
Prediction_Lineata_V <- predict(Lineata_Beta_V,
                                newdata = NewData_Lineata_V,
                                type = "response",
                                se.fit = TRUE) 

NewData_Lineata_V$Fit <- Prediction_Lineata_V$fit
NewData_Lineata_V$SE <- Prediction_Lineata_V$se.fit

# Confidence intervals:
NewData_Lineata_V <- NewData_Lineata_V |>
  mutate(
    Lower = pmax(0, Fit - 1.96 * SE),
    Upper = pmin(1, Fit + 1.96 * SE))

# Plot with aesthetics:
Fig_4A_8 <- ggplot(PopLev_V, aes(x = EcolPC1, y = Lineata_Frequency)) +
  geom_point(size = 4, shape = 21, fill = "#A09D97", color = "black") +
  geom_ribbon(data = NewData_Lineata_V, aes(x = EcolPC1, ymin = Lower, ymax = Upper), fill = "grey80", alpha = 0.4, inherit.aes = FALSE) +
  geom_line(data = NewData_Lineata_V, aes (y = Fit), linewidth = 1.5, color = "black") +
  labs(title = "Vigo",
       x = "Ecological PC1",
       y = expression(paste("Freq. of ", italic("Lineata")))) +
  theme_cd(18) + 
  theme(plot.title = element_text(color = "black", size = 18)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  #scale_x_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Same for the x axis
  scale_x_reverse(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_4A_8 # See the plot

### 2.- Temporal stability (Vigo) ----

# Subsection aim: Present the temporal stability of fulva and lineata frequencies in figures

# ///

###### 2.1.- Fulva snails ----

# Testing differences (beta-binomial):
Fulva_Year_Beta <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                                       DistanceRiver_Absolute * Year,
                                     family = betabinomial(link = "logit"),
                                     data = TempSer)

# Run diagnostics:
SimulatedResiduals_Fulva_Year_Beta <- simulateResiduals(Fulva_Year_Beta)
plot(SimulatedResiduals_Fulva_Year_Beta)
# OUTPUT: The fit is good

# Inspect the results' table:
summary(Fulva_Year_Beta)

# Likelihood ratio test:
Fulva_Year_Beta_Reduced <- glmmTMB(cbind(Fulva_Counts, n-Fulva_Counts) ~ DistanceRiver_Absolute + Year, 
                         family = betabinomial(link = "logit"),
                         data = TempSer) # Model without interaction

anova(Fulva_Year_Beta, Fulva_Year_Beta_Reduced) # Formally assess the differences between models
# OUTPUT: There are differences between years

# Extract the probabilities and confidence intervals of the beta-binomial model:
Predicted_Fulva_Year <- emmeans(Fulva_Year_Beta, 
                ~ DistanceRiver_Absolute | Year,
                at = list(DistanceRiver_Absolute = unique(TempSer$DistanceRiver_Absolute)),
                type = "response")

Predicted_Fulva_Year_Grid <- as.data.frame(Predicted_Fulva_Year) # Convert to a dataframe to plot it afterwards

# Define custom colors:
Year_Colors <- c(
  "1979" = "#2B31D9",
  "2022" = "#D19215"
)

# Plot with aesthetics:
Fig_4B_1 <- ggplot(TempSer, aes(x = DistanceRiver_Absolute, y = Fulva_Frequency, color = Year)) +
  geom_point(aes(fill = Year), size = 4, alpha = 1, show.legend = TRUE, color = "black", shape = 21) +
  geom_line(data = Predicted_Fulva_Year_Grid, aes(x = DistanceRiver_Absolute, y = prob, color = Year), linewidth = 1.5) +
  geom_ribbon(data = Predicted_Fulva_Year_Grid, aes(x = DistanceRiver_Absolute, y = prob, ymin = asymp.LCL, ymax = asymp.UCL, fill = Year),
              alpha = 0.2, color = NA) +
  labs(title = "",
       x = "Distance from river mouth (m)",
       y = expression(paste("Freq. of ", italic("Fulva")))) +
  scale_color_manual(values = Year_Colors) +
  scale_fill_manual(values = Year_Colors) +
  theme_cd(18) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  scale_x_reverse(breaks = c("0" = 0, "20,000" = 20000, "40,000" = 40000, "60,000" = 60000))
  
Fig_4B_1 # View the plot

#### 2.2.- Lineata snails ----

# Testing differences (beta-binomial):
Lineata_Year_Beta <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                                       DistanceRiver_Absolute * Year,
                                     family = betabinomial(link = "logit"),
                                     data = TempSer)

# Run diagnostics: 
SimulatedResiduals_Lineata_Year_Beta <- simulateResiduals(Lineata_Year_Beta)
plot(SimulatedResiduals_Lineata_Year_Beta)
# OUTPUT: The fit is reasonably good, although with some problems

# Inspect the results' table:
summary(Lineata_Year_Beta)

# Likelihood ratio test:
Lineata_Year_Beta_Reduced <- glmmTMB(cbind(Lineata_Counts, n-Lineata_Counts) ~ DistanceRiver_Absolute + Year, 
                                        family = betabinomial(link = "logit"),
                                        data = TempSer) # Model without interaction

anova(Lineata_Year_Beta, Lineata_Year_Beta_Reduced) # Formally assess the differences between models
# OUTPUT: There are no differences between years

# Extract the probabilities and confidence intervals of the beta-binomial model:
Predicted_Lineata_Year <- emmeans(Lineata_Year_Beta, 
                ~ DistanceRiver_Absolute | Year,
                at = list(DistanceRiver_Absolute = unique(TempSer$DistanceRiver_Absolute)),
                type = "response")

Predicted_Lineata_Year_Grid <- as.data.frame(Predicted_Lineata_Year) # Convert to a dataframe to plot it afterwards

# Plot with aesthetics:
Fig_4B_2 <- ggplot(TempSer, aes(x = DistanceRiver_Absolute, y = Lineata_Frequency, color = Year)) +
  geom_point(aes(fill = Year), size = 4, alpha = 1, show.legend = TRUE, color = "black", shape = 21) +
  geom_line(data = Predicted_Lineata_Year_Grid, aes(x = DistanceRiver_Absolute, y = prob, color = Year), linewidth = 1.5) +
  geom_ribbon(data = Predicted_Lineata_Year_Grid, aes(x = DistanceRiver_Absolute, y = prob, ymin = asymp.LCL, ymax = asymp.UCL, fill = Year),
              alpha = 0.2, color = NA) +
  labs(title = "",
       x = "Distance from river mouth (m)",
       y = expression(paste("Freq. of ", italic("Lineata")))) +
  scale_color_manual(values = Year_Colors) +
  scale_fill_manual(values = Year_Colors) +
  theme_cd(18) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) + # Remove the annoying leading 0s
  scale_x_reverse(breaks = c("0" = 0, "20,000" = 20000, "40,000" = 40000, "60,000" = 60000))

Fig_4B_2 # View the plot

### 3.- Test of isolation by distance ----

# Subsection aim: Are differences in fulva and lineata frequencies across localities explained by geographic or ecological distance?

# ///

#### 3.1.- Morph composition differences between localities (Fulva and Lineata) ----

# Select the variables needed to calculate the index:
PopLev_BC <-PopLev %>% 
  dplyr::select(LocName, Fulva_Counts,
                Lineata_Counts)

# Keep the morphs only:
PopLev_BC_Sp <- PopLev_BC[, c("Fulva_Counts", "Lineata_Counts")]

# Calculate the Bray-Curtis dissimilarity index:
Dist_Morph <- vegdist(PopLev_BC_Sp, method = "bray")

# Inspect the distribution of values (just out of curiosity):
hist(Dist_Morph)

# Turn the data into a matrix format:
Dist_Morph_Matrix <- as.matrix(Dist_Morph)

# Assign row and column names to the matrix:
rownames(Dist_Morph_Matrix) <- PopLev_BC$LocName
colnames(Dist_Morph_Matrix) <- PopLev_BC$LocName

# Turn into a long format:
Dist_Morph_Matrix_Long <- melt(Dist_Morph_Matrix, 
                    varnames = c("Loc1","Loc2"), 
                    value.name = "BrayCurtisDiss")
Dist_Morph_Matrix_Long # Ready to be merged with other distances

#### 3.2.- Geographic distances between localities ----

# Select the variables needed to calculate the differences:
PopLev_Geo <- PopLev %>% 
  dplyr::select(LocName, Latitude, Longitude)

# Extract the coordinates only:
Coordinates <- PopLev_Geo[, c("Longitude", "Latitude")] # Switch the order for geosphere

# Calculate geographic distances (returns a matrix):
Dist_Geo_Matrix <- distm(Coordinates, fun = distHaversine)

# Transform them into kilometers:
#Dist_Geo_Matrix <- Dist_Geo_Matrix / 1000

# Assign row and column names to the matrix:
rownames(Dist_Geo_Matrix) <-PopLev_Geo$LocName
colnames(Dist_Geo_Matrix) <- PopLev_Geo$LocName

# Turn into a long format:
Dist_Geo_Matrix_Long <- melt(Dist_Geo_Matrix,
                      varnames = c("Loc1", "Loc2"),
                      value.name = "GeoDist")
Dist_Geo_Matrix_Long # Ready to be merged with other distances

#### 3.3.- Ecological distances between localities ----

# Select the variables needed to calculate the differences:
PopLev_Ecol <- PopLev %>% 
  dplyr::select(LocName, EcolPC1)

# Calculate Euclidean distances:
Dist_Ecol <- dist(PopLev_Ecol$EcolPC1, method = "euclidean")

# Inspect the distribution of values (just out of curiosity):
hist(Dist_Ecol)

# Turn the data into a matrix format:
Dist_Ecol_Matrix <- as.matrix(Dist_Ecol)

# Assign row and column names to the matrix:
rownames(Dist_Ecol_Matrix) <- PopLev_Ecol$LocName
colnames(Dist_Ecol_Matrix) <- PopLev_Ecol$LocName

# Turn into a long format:
Dist_Ecol_Matrix_Long <- melt(Dist_Ecol_Matrix, 
                       varnames = c("Loc1","Loc2"), 
                       value.name = "EcolDist")
Dist_Ecol_Matrix_Long # Ready to be merged with other distances

#### 3.4.- Fusing all distances into a dataset and creating a plot ----

# Fuse all distances into a dataset:
Distances <- list(
  Dist_Ecol_Matrix_Long,
  Dist_Geo_Matrix_Long,
  Dist_Morph_Matrix_Long
)

view(Distances)

# Remove the duplicate localities:
Distances_Reduced <- reduce(Distances, left_join, by = c("Loc1", "Loc2"))
view(Distances_Reduced)

# Remove duplicate comparisons:
Loc_Order <- rownames(Dist_Ecol_Matrix) # Set an order for the localities to ensure consistency

Distances_Reduced <- Distances_Reduced %>% 
  filter(match(Loc1, Loc_Order) < match(Loc2, Loc_Order))
view(Distances_Reduced)


# Plot of BrayCurtisDiss ~ GeoDist:
Fig_4C_1 <- ggplot(Distances_Reduced, aes(x = GeoDist, y = BrayCurtisDiss)) +
  geom_point(size = 3, alpha = 1, fill = "#91633C", show.legend = TRUE, shape = 21) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.5) +
  labs(title = "",
       x = "Pairwise geographical \n distances (m)",
       y = "Pairwise morph \n composition distances") +
  theme_cd(18) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1)) +  # Remove the annoying leading 0s
  scale_x_continuous(breaks = c("0" = 0, "20,000" = 20000, "40,000" = 40000, "60,000" = 60000))

Fig_4C_1 # View the plot

# Plot of BrayCurtisDiss ~ EcolDist:
Fig_4C_2 <- ggplot(Distances_Reduced, aes(x = EcolDist, y = BrayCurtisDiss)) +
  geom_point(size = 3, alpha = 1, fill = "#69b3a2", show.legend = TRUE, shape = 21) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1.5) +
  labs(title = "",
       x = "Pairwise ecological \n distances (Ecological PC1)",
       y = "Pairwise morph \n composition distances") +
  theme_cd(18) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1))  # Remove the annoying leading 0s

Fig_4C_2 # View the plot

#### 3.5.- Assess differences statistically ----

# Mantel test:
mantel.test(Dist_Morph_Matrix, Dist_Ecol_Matrix, nperm = 999) # OUTPUT: P = 0.001
mantel.test(Dist_Morph_Matrix, Dist_Geo_Matrix, nperm = 999) # OUTPUT: P > 0.05

# Multiple regression on distance matrices:
MRM_Dist <- MRM(BrayCurtisDiss ~ EcolDist + GeoDist, data = Distances_Reduced, nperm = 999)

MRM_Dist
# OUTPUT: Ecological distances are significant, geographical distances aren't

### 4.- SLiM simulations ----

# NOTE: This code corresponds only to that needed to generate the figures for the manuscript and supplementary materials. All raw data as extracted directly from SLiM is available at https://gitlab.com/elcortegano/rias_gallegas.
# Subsection aim: Can mechanisms other than selection explain the existence of the cline?

# ///

#### 4.1.- Replicate slopes ----

# Set generation categories in the simulations (for dashed lines):
burnin <- 10 # Burn-in generations
t_expansion <- 610 # Where expansion finishes

# Figure of median slopes and 95% CIs:
Fig_4D <- RepDef %>%
  ggplot(aes(x = generation)) +
  geom_vline(xintercept = burnin, linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = t_expansion, linetype = "dashed", linewidth = 0.5) +
  geom_line(aes(y = b, color = factor(transect))) +
  geom_line(aes(y = bl, color = factor(transect)), linetype = "dotted") +
  geom_line(aes(y = bh, color = factor(transect)), linetype = "dotted") +
  labs(x = "Generation",
       y = "Slope") +
  scale_color_brewer(palette = "Dark2") +
  annotate("text", x = burnin + 25, y = -4.2, label = "Start expansion", color = "#000000", size = 6, hjust = 0) +
  annotate("text", x = t_expansion + 25, y = -4.2, label = "End expansion", color = "#000000", size = 6, hjust = 0) +
  theme_cd(18) +
  scale_y_continuous(breaks = c("-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Remove the annoying leading 0s
  scale_x_continuous(breaks = c("0" = 0, "250" = 250, "500" = 500, "750" = 750, "1,000" = 1000, "1,250" = 1250))

Fig_4D

### 5.- Composite figure (Figure 4) ----

#### 5.1.- Fig. 4A ----

# Set margins:
Margin_4A <- theme(plot.margin = margin(t = 0, r = 2, b = 2, l = 0, "mm"))

Fig_4A_1 <- Fig_4A_1 + Margin_4A 
Fig_4A_2 <- Fig_4A_2 + Margin_4A
Fig_4A_3 <- Fig_4A_3 + Margin_4A
Fig_4A_4 <- Fig_4A_4 + Margin_4A
Fig_4A_5 <- Fig_4A_5 + Margin_4A
Fig_4A_6 <- Fig_4A_6 + Margin_4A
Fig_4A_7 <- Fig_4A_7 + Margin_4A
Fig_4A_8 <- Fig_4A_8 + Margin_4A

# Remove axis texts:
Fig_4A_1 <- Fig_4A_1 + theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank())

Fig_4A_2 <- Fig_4A_2 + theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.title.x = element_blank())

Fig_4A_3 <- Fig_4A_3 + theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.title.x = element_blank())

Fig_4A_4 <- Fig_4A_4 + theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.title.x = element_blank())

Fig_4A_5 <- Fig_4A_5 + theme(
  plot.title = element_blank())

Fig_4A_6 <- Fig_4A_6 + theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  plot.title = element_blank())

Fig_4A_7 <- Fig_4A_7 + theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  plot.title = element_blank())

Fig_4A_8 <- Fig_4A_8 + theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank(),
  plot.title = element_blank())

# Figure 4A - top:
Fig_4A_top <- plot_grid(Fig_4A_1, Fig_4A_2, Fig_4A_3, Fig_4A_4,
                           ncol = 4,
                           align = "h",
                           axis = "tb", rel_widths = c(1.2, 0.9, 0.9, 0.9))

# Figure 4A - bottom:
Fig_4A_bottom <- plot_grid(Fig_4A_5, Fig_4A_6, Fig_4A_7, Fig_4A_8,
                        ncol = 4,
                        align = "h",
                        axis = "tb", rel_widths = c(1.2, 0.9, 0.9, 0.9))

# Figure 4A:
Fig_4A <- plot_grid(Fig_4A_top, Fig_4A_bottom, 
                    ncol = 1, align = "v", axis = "lr",
                    rel_heights = c(1, 1))

# Add the tag to Figure 4A:
Fig_4A_Tag <- tags_cd(Fig_4A, "A", size = 22, x = 0.01, y = 0.99)

# Remove axis texts and legend:
Fig_4B_1 <- Fig_4B_1 + theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  legend.position = "none")

#### 5.2.- Fig. 4B ----

# Set margins:
Margin_4B <- theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 4, "mm"))

Fig_4B_1 <- Fig_4B_1 + Margin_4B
Fig_4B_2 <- Fig_4B_2 + Margin_4B

# Remove legend:
Fig_4B_2 <- Fig_4B_2 + theme(
  legend.position = "none")

# Merge them:
Fig_4B <- plot_grid(Fig_4B_1, Fig_4B_2,
                    ncol = 1, 
                    align = "v",
                    axis = "lr",
                    rel_heights = c(1, 1))

# Add tag:
Fig_4B_Tag <- tags_cd(Fig_4B, "B", size = 22, x = 0.02, y = 0.99)

#### 5.3.- Fig. 4C ----

# Set margins:
Margin_4C <- theme(plot.margin = margin(t = 2, r = 2, b = 0, l = 0, "mm"))

Fig_4C_1 <- Fig_4C_1 + Margin_4C
Fig_4C_2 <- Fig_4C_2 + Margin_4C

# Remove axis texts:
Fig_4C_2 <- Fig_4C_2 + theme(
  axis.text.y = element_blank(),
  axis.title.y = element_blank())

# Add tag:
Fig_4C_1_Tag <- tags_cd(Fig_4C_1, "C", size = 22, x = 0.04, y = 0.98)

#### 5.4.- Fig. 4D ----

# Remove legend:
Fig_4D <- Fig_4D + theme(
  legend.position = "none")

# Set margins:
Margin_4D <- theme(plot.margin = margin(t = 0, r = 8, b = 0, l = 0, "mm"))

Fig_4D <- Fig_4D + Margin_4D

# Add tag: 
Fig_4D_Tag <- tags_cd(Fig_4D, "D", size = 22, x = 0.04, y = 0.98)

#### 5.5.- Top and bottom rows ----

# Top row: 
Fig_4_TopRow <- plot_grid(Fig_4A_Tag, Fig_4B_Tag,
                           ncol = 2, align = "h", axis = "tb",
                           rel_widths = c(4, 2))

# Bottom row: 
Fig_4_BottomRow <- plot_grid(Fig_4C_1_Tag, Fig_4C_2, Fig_4D_Tag,
                             ncol = 3, align = "h", axis = "tb",
                             rel_widths = c(1.3, 1, 1.5))

# Set margin:
Margin_BottomRow <- theme(plot.margin = margin(t = 4, r = 0, b = 0, l = 0, "mm"))

Fig_4_BottomRow <- Fig_4_BottomRow + Margin_BottomRow

# Full figure:
Figure4 <- plot_grid(Fig_4_TopRow, Fig_4_BottomRow,
                     ncol = 1, align = "h", axis = "tb",
                     rel_widths = c(2, 1))

# Save in .svg to export to InkScape:
ggsave("Figure4.svg",
       Figure4,
       width = 460, height = 255, units = "mm")

## Supplementary analyses ----

### 1.- Map figure with finer-grained classes ----

# Preparing the dataset for the fulva, lineata, and rest pie charts
PopLev_PieCharts <- PopLev %>% # Selecting the variables that I will use
  dplyr::select(Ria, LocNumCumulative, LocName, Fulva_Frequency, 
                Lutea_Frequency, Fusca_Frequency, Lineata_Frequency, 
                Minority_Frequency)

PopLev_PieCharts_long <- PopLev_PieCharts %>% 
  pivot_longer(
    cols = ends_with("_Frequency"),
    names_to = "Morph",
    values_to = "Frequency"
  )

# Plot
PieLFFLM <- ggplot(PopLev_PieCharts_long, aes(x = "", y = Frequency, fill = Morph)) +
  geom_bar(stat = "identity", width = 1, color = "black",linewidth = 0.2) +
  coord_polar(theta = "y") +
  facet_wrap(~ LocNumCumulative) +
  scale_fill_manual(
    values = c(
      "Fulva_Frequency" = "#D6D2A9",
      "Fusca_Frequency" = "#856E42",
      "Lineata_Frequency" = "#A09D97",
      "Lutea_Frequency" = "#F0DC75",
      "Minority_Frequency" = "#95B877"
    )
  ) +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold")
  )

PieLFFLM

# Save the plot
ggsave("CD_DA_Out_PL_SCE_PieLFFLM.svg", # Choose .svg format for later aesthetic additions
       plot = PieLFFLM,
       width = 7,
       height = 6)

### 2.- Supplementary SLiM figures ----

# Set generation categories in the simulations (for dashed lines):
burnin <- 10 # Burn-in generations
t_expansion <- 610 # Where expansion finishes

#### 1.1.- Slope calculated over the mean frequency across replicates (i.e, over the expected frequencies; variation is shown only among transects) ----
Fig_SM_Expected <- Expected %>%
  ggplot(aes(x = generation)) +
  geom_vline(xintercept = burnin, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = t_expansion, linetype = "dashed", linewidth = 0.3) +
  geom_line(aes(y = b, colour = factor(transect))) +
  labs(x = "Generation",
       y = "Slope") +
  scale_color_brewer(palette = "Dark2") +
  annotate("text", x = burnin + 20, y = 0.9, label = "Start expansion", color = "#000000", size = 4, hjust = 0) +
  annotate("text", x = t_expansion + 20, y = 0.9, label = "End expansion", color = "#000000", size = 4, hjust = 0) +
  theme_cd(18) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = c("-1" = -1, "-0.5" = -0.5, "0" = 0, "0.5" = 0.5, "1" = 1)) + # Remove the annoying leading 0s
  scale_x_continuous(breaks = c("0" = 0, "250" = 250, "500" = 500, "750" = 750, "1,000" = 1000, "1,250" = 1250))


Fig_SM_Expected # View the plot

# Save the plot:
ggsave("CD_MS_SM_Fig_S11.svg",
       plot = Fig_SM_Expected,
       width = 6,
       height = 5)

#### 1.2.- Maintenance of significant slopes ----
Fig_SM_Cons <- ConsGen %>%
  ggplot(aes(x = run_length)) +
  geom_point(aes(y = score/1210, colour = factor(T_type)), position = position_dodge2(width = 0.95)) +
  geom_errorbar(aes(ymin = sl/1210, ymax = sh/1210, colour = factor(T_type)), position = position_dodge2(width = 0.95)) +
  scale_y_continuous(name = "% of time", limits = c(0, 0.23), breaks = seq(0, 0.2, 0.05), labels = scales::label_percent(accuracy = 1)) +
  scale_x_continuous(name = expression(atop("Number of consecutive",
                                            "generations where " * italic(P) < 0.05)),
                     limits = c(0.5, 11.5), breaks = 1:11, labels = c("1","2","3","4","5","6","7","8","9","10",">10")) +
  scale_color_brewer(palette = "Dark2") +
  theme_cd(18)

Fig_SM_Cons # View the plot

# Save the plot:
ggsave("CD_MS_SM_Fig_S12.svg",
       plot = Fig_SM_Cons,
       width = 7,
       height = 5)

#### 1.3.- Among-replicates variation (RepDef) but for the 'drift' and 'migration' cases (alternative scenarios) ----
Fig_SM_AltScen <- AltScen %>%
  ggplot(aes(x = generation)) +
  geom_vline(xintercept = burnin, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = t_expansion, linetype = "dashed", linewidth = 0.3) +
  geom_line(aes(y = b, color = factor(transect))) +
  geom_line(aes(y = bl, color = factor(transect)), linetype = "dotted") +
  geom_line(aes(y = bh, color = factor(transect)), linetype = "dotted") +
  labs(x = "Generation",
       y = "Slope") +
  scale_color_brewer(palette = "Dark2") +
  annotate("text", x = burnin + 20, y = -4.2, label = "Start expansion", color = "#000000", size = 5, hjust = 0) +
  annotate("text", x = t_expansion + 20, y = -4.2, label = "End expansion", color = "#000000", size = 5, hjust = 0) +
  facet_grid(. ~ scenario) +
  theme_cd(18) +
  theme(legend.position = "none",
        plot.margin = margin(t = 10, b =10, l = 10, r = 30)) +
  scale_y_continuous(limits = c(-4.4, 4.4),
                     breaks = c("-2.5" = -2.5, "0" = 0, "2.5" = 2.5)) + # Remove the annoying leading 0s
  scale_x_continuous(breaks = c("0" = 0, "250" = 250, "500" = 500, "750" = 750, "1,000" = 1000, "1,250" = 1250))

Fig_SM_AltScen # View the plot

# Save the plot:
ggsave("CD_MS_SM_Fig_S13.svg",
       plot = Fig_SM_AltScen,
       width = 11,
       height = 6)
