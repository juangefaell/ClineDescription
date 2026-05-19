############### ClineDescription - Selective Pressures #################

# BASIC INFO ----
 
# Researcher: [AuthorName]
# Section aim: What are the selective pressures that drive the evolution of different traits? 
# Last update: 2026-05-19

# SETUP SECTION ----

## 1.- Packages ----
library(tidyverse) # Basic plots and statistical analyses
library(lme4) # LME models
library(pavo) # Reflectance spectrometry analyses
library(DHARMa) # Model diagnostics
library(heplots) # Partial ffect size
library(cowplot) # Composite figures
library(grid) # Composite figures 

## 2.- Dataset loading ----
setwd("") # Set the directory of the dataset file
read.csv("CD_DA_RD_ReflectanceSpectrometry.csv") 
read.csv("CD_DA_RD_PopulationLevel.csv") 
read.csv("CD_DA_RD_IndividualLevel.csv")

## 3.- Data handling and overview ----

# ReflectanceSpectrometry dataset:
RefSpec <- read.csv("CD_DA_RD_ReflectanceSpectrometry.csv") %>% # Set shortcut name and convert the `ReflectanceSpectrometry` dataset to the appropriate format
  drop_na() %>% # Remove NAs
  as.rspec(lim = c(300, 700))  # Convert to 'pavo' format and select the visible range only

# PopulationLevel dataset:
PopLev <- read.csv("CD_DA_RD_PopulationLevel.csv") %>%
  mutate(
    Ria = factor(Ria, levels = c("MurosNoia", "Arousa","Pontevedra", "Vigo")), # Set an order from north to south
    Coast = factor(Coast, levels = c("South", "North")),
    Cluster = factor(Cluster, levels = c("WaveSheltered", "WaveIntermediate", "WaveExposed"))) # Convert variables into adequate formats
str(PopLev) # Check the structure of the data

# IndividualLevel dataset:
IndLev <- read.csv("CD_DA_RD_IndividualLevel.csv") %>%
  mutate(
    Ria = factor(Ria, levels = c("MurosNoia", "Arousa", "Pontevedra", "Vigo")),
    Coast = factor(Coast, levels = c("South", "North")),
    Cluster = factor(Cluster, levels = c("WaveSheltered", "WaveIntermediate", "WaveExposed"), ordered = TRUE),
    Lineata = if_else(ColorMorph == "Lineata", "1", "0"), # Create a Lineata variable
    Fulva = if_else(ColorMorph == "Fulva", "1", "0"), # Create a Fulva variable
    SculpturedShell = if_else(ShellSculpture == "Jugosa", "1", "0"), # Create the SculpturedShell variable (Sculptured snail or not)
    Scars = factor(Scars, levels = c("Absent", "Present")),
    RW1_Z = scale(RW1), # Scale the RW1 variable
    ShellLength_Z = scale(ShellLength), # Scale the ShellLength variable
    ShellThickness_Rel_Z = scale(ShellThickness_Rel)) # Scale the ShellThickness_Rel variable
str(IndLev) # Check the structure of the data

## 4.- Themes for figures ----

# Theme figures, v. 2:
theme_cd2 <- function(BaseSize = 24, BaseFamily = "sans") {
  theme_bw(base_size = BaseSize, base_family = BaseFamily) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(), 
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

### 1.- Color as an adaptation for camouflage ----

# Subsection aim: Are shell colors (lineata and fulva) adaptations for camouflage in their respective environments?

# ///

#### 1.1.- Data curation ----

# Identify the position of negative values, as these can distort the overall spectra:
RefSpec %>% 
  plot(type = "o") %>% 
  abline(h = 0, lty = 3) 
# OUTPUT: No negative values within the visible range.

# Correct negative values by converting them to 0 and apply smoothing (span = 0.6):
RefSpec_Zero <- RefSpec %>%
  procspec(., opt = c("smooth"), # Apply LOESS smoothing
           span = 0.6, # Neighboring size with span = 0.6 for the LOESS smoothing
           fixneg = "zero") # Turn negative values into to zero

# Confirm that the selected span works well:
RefSpec_Zero %>%
  subset("Fulva") %>% # Change the object at will
  plotsmooth(minsmooth = 0.05, 
             maxsmooth = 1, 
             curves = 4, 
             ask = FALSE) 
# OUTPUT: It works well, since it distinguishes between morphs, yet softens the noise in the data (either way, there is not too much noise in the data) 

# Visually inspect the selected parameters on different morphs:
RefSpec_Zero %>% 
  subset("Fulva") %>% # Change between objects
  plot(type = "o", # Select the type of graph
       ylim = c(0, 40)) # Set plot limits at will
  

#### 1.2.- Visual models of a fish predator -----

##### 1.2.1.- Lineata, Fulva, and barnacles ----

# Extract the names of the objects from the column names:
RefSpec_Col <- colnames(RefSpec_Zero) %>% # Extract the column names
  strsplit("_") %>% # Split the column names based on '_'
  sapply(function(x) x[2]) %>% # Select the second element of the column names (which contains the names of color morphs)
  .[-1] # Remove the first column, which corresponds to 'wl'

# Select the objects to plot:
Exposed <- RefSpec_Zero %>%
  dplyr::select(wl, which(RefSpec_Col %in% c("Lineata", "Fulva", "Barnacles")) + 1)  # +1 because column 1 is 'wl'

# Create a visual model for a fish ("Rhinecanthus")
VisModel_Exposed <- vismodel(Exposed, visual ="rhinecanthus", achromatic= "all", relative = FALSE)

# Check the summary of the visual model:
summary(VisModel_Exposed)

# Calculate the chromatic distances between the morphs and the barnacles:
coldist(VisModel_Exposed, 
        noise = "neural", achromatic = TRUE, n = c(1, 2, 4),
        weber = 0.05, weber.achro = 0.05, subset = "Barnacles")

# Set the names to identify objects:
names_vec_exp <- rownames(VisModel_Exposed)

# Establish which are the colors to calculate the chromatic distances:
colors_exp <- sapply(strsplit(names_vec_exp, "_"), function(z) z[2])

# Compare chromatic distances by bootstrapping:
ChromaticDistances_Exposed <- bootcoldist(VisModel_Exposed, 
                                          by = colors_exp,
                                          n = c(1, 2, 4),
                                          weber = 0.05,
                                          weber.achro = 0.05)

ChromaticDistances_Exposed # View the result

# Create a data frame of the chromatic distances:
ChromaticDistances_Exposed_df <- data.frame(
  x = 1:nrow(ChromaticDistances_Exposed),
  Mean = ChromaticDistances_Exposed[, 1],
  Lower = ChromaticDistances_Exposed[, 2],
  Upper = ChromaticDistances_Exposed[, 3],
  Label = rownames(ChromaticDistances_Exposed))

ChromaticDistances_Exposed_df # View the data frame

# Remove the last row (not very interesting):
ChromaticDistances_Exposed_df <- ChromaticDistances_Exposed_df %>% 
  filter(Label != "Fulva-Lineata")

# Set a color to the classes in the plot below:
Colors_ChromaticDistances_Exposed <- c(
  "Barnacles-Fulva" = "#D6D2A9",
  "Barnacles-Lineata" = "#A09D97")

# Chromatic contrast (exposed environment) (Figure 6A top):
Fig_6A_T <- ggplot(ChromaticDistances_Exposed_df, aes(x = factor(x), y = Mean, fill = Label)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "red") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, linewidth = 0.9) +
  geom_point(shape = 21, size = 5, stroke = 1.2) +
  scale_fill_manual(values = Colors_ChromaticDistances_Exposed) +
  scale_x_discrete(breaks = ChromaticDistances_Exposed_df$x, labels = c(
    expression(vs.~italic(Ful.)),
    expression(vs.~italic(Lin.))), 
    expand = expansion(mult = 0.3)) +
  labs(
    x = "",
    y = "Chrom. contrast (dS)",
    title = "Barnacle-covered rock"
  ) +
  coord_cartesian(ylim = c(0, 4)) +
  theme_cd2(18) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(colour = "black", margin = margin(b = 6), size = 16))

Fig_6A_T # View the plot

# Achromatic contrast:
AchromaticDistances_Exposed_df <- data.frame(
  x = 1:nrow(ChromaticDistances_Exposed),
  Mean = ChromaticDistances_Exposed[, 4],
  Lower = ChromaticDistances_Exposed[, 5],
  Upper = ChromaticDistances_Exposed[, 6],
  Label = rownames(ChromaticDistances_Exposed))

AchromaticDistances_Exposed_df # View the data frame

# Remove the last row (not very interesting):
AchromaticDistances_Exposed_df <- AchromaticDistances_Exposed_df %>% 
  filter(Label != "Fulva-Lineata")

# Chromatic contrast (exposed environment) (Figure 6A top):
ggplot(AchromaticDistances_Exposed_df, aes(x = factor(x), y = Mean, fill = Label)) +
  #geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "red") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, linewidth = 0.9) +
  geom_point(shape = 21, size = 5, stroke = 1.2) +
  scale_fill_manual(values = Colors_ChromaticDistances_Exposed) +
  scale_x_discrete(breaks = AchromaticDistances_Exposed_df$x, labels = c(
    expression(vs.~italic(Ful.)),
    expression(vs.~italic(Lin.))), 
    expand = expansion(mult = 0.3)) +
  labs(
    x = "",
    y = "Achromatic contrast",
    title = "Barnacle-covered rocks"
  ) +
  #coord_cartesian(ylim = c(0, 4)) +
  theme_cd(18) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(colour = "black", margin = margin(b = 6), size = 16))

##### 1.2.2.- Lineata, Fulva, and barnacle-free rocks (Which morph is more cryptic against a wave-sheltered background? Background 1) ----

# Extract the names of the objects from the column names:
RefSpec_Col <- colnames(RefSpec_Zero) %>% # Extract the column names
  strsplit("_") %>% # Split the column names based on '_'
  sapply(function(x) x[2]) %>% # Select the fifth element of the column names (which contains the names of color morphs)
  .[-1] # Remove the first column, which corresponds to 'wl'

# Select the objects to plot:
Sheltered <- RefSpec_Zero %>%
  dplyr::select(wl, which(RefSpec_Col %in% c("Lineata", "Fulva", "Rock")) + 1)  # +1 because column 1 is 'wl'

# Create a visual model for a fish ("Rhinecanthus")
VisModel_Sheltered <- vismodel(Sheltered, visual ="rhinecanthus", achromatic= "all", relative = FALSE)

# Check the summary of the visual model:
summary(VisModel_Sheltered)

# Calculate the chromatic distances between the morphs and rocks:
coldist(VisModel_Sheltered, 
        noise = "neural", achromatic = TRUE, n = c(1, 2, 4),
        weber = 0.05, weber.achro = 0.05, subset = "Rock")

# Set the names to identify objects:
names_vec_shelt <- rownames(VisModel_Sheltered)

# Establish which are the colors to calculate the chromatic distances:
colors_shelt <- sapply(strsplit(names_vec_shelt, "_"), function(z) z[2])

# Compare chromatic distances by bootstrapping:
ChromaticDistances_Sheltered <- bootcoldist(VisModel_Sheltered, 
                       by = colors_shelt,
                       n = c(1, 2, 4),
                       weber = 0.05,
                       weber.achro = 0.05)

ChromaticDistances_Sheltered # View the result

# Create a data frame of the chromatic distances:
ChromaticDistances_Sheltered_df <- data.frame(
  x = 1:nrow(ChromaticDistances_Sheltered),
  Mean = ChromaticDistances_Sheltered[, 1],
  Lower = ChromaticDistances_Sheltered[, 2],
  Upper = ChromaticDistances_Sheltered[, 3],
  Label = rownames(ChromaticDistances_Sheltered))

ChromaticDistances_Sheltered_df # View the data frame

# Remove the last row (not very interesting):
ChromaticDistances_Sheltered_df <- ChromaticDistances_Sheltered_df %>% 
  filter(Label != "Fulva-Lineata")

# Set a color to the classes in the plot below:
Colors_ChromaticDistances_Sheltered <- c(
  "Fulva-Rock" = "#D6D2A9",
  "Lineata-Rock" = "#A09D97")

# Chromatic contrast (sheltered environment) (Figure 6A bottom):
Fig_6A_B <- ggplot(ChromaticDistances_Sheltered_df, aes(x = factor(x), y = Mean, fill = Label)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "red") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, linewidth = 0.9) +
  geom_point(shape = 21, size = 5, stroke = 1.2) +
  scale_fill_manual(values = Colors_ChromaticDistances_Sheltered) +
  scale_x_discrete(breaks = ChromaticDistances_Sheltered_df$x, labels = c(
    expression(vs.~italic(Ful.)),
    expression(vs.~italic(Lin.))), 
    expand = expansion(mult = 0.3)) +
  labs(
    x = "",
    y = "Chrom. contrast (dS)",
    title = "Barnacle-free rock"
  ) +
  coord_cartesian(ylim = c(0, 4)) +
  theme_cd2(18) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(colour = "black", margin = margin(b = 6), size = 16))

Fig_6A_B # View the plot

# Achromatic contrast:
AchromaticDistances_Sheltered_df <- data.frame(
  x = 1:nrow(ChromaticDistances_Sheltered),
  Mean = ChromaticDistances_Sheltered[, 4],
  Lower = ChromaticDistances_Sheltered[, 5],
  Upper = ChromaticDistances_Sheltered[, 6],
  Label = rownames(ChromaticDistances_Sheltered))

AchromaticDistances_Sheltered_df # View the data frame

# Remove the last row (not very interesting):
AchromaticDistances_Sheltered_df <- AchromaticDistances_Sheltered_df %>% 
  filter(Label != "Fulva-Lineata")

# Chromatic contrast (exposed environment) (Figure 6A top):
ggplot(AchromaticDistances_Sheltered_df, aes(x = factor(x), y = Mean, fill = Label)) +
  #geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "red") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.1, linewidth = 0.9) +
  geom_point(shape = 21, size = 5, stroke = 1.2) +
  scale_fill_manual(values = Colors_ChromaticDistances_Sheltered) +
  scale_x_discrete(breaks = AchromaticDistances_Sheltered_df$x, labels = c(
    expression(vs.~italic(Ful.)),
    expression(vs.~italic(Lin.))), 
    expand = expansion(mult = 0.3)) +
  labs(
    x = "",
    y = "Achromatic contrast (dL)",
    title = "Barnacle-free rocks"
  ) +
  #coord_cartesian(ylim = c(0, 4)) +
  theme_cd(18) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(colour = "black", margin = margin(b = 6), size = 16))


### 2.- Size and shape as adaptations to crabs and waves ----

# Subsection aim: Are size and shape adaptations to wave action and crab predation (in exposed and sheltered environments, respectively)?

# ///

#### 2.1.- MANOVA on Shell length and shell shape as a function of Wave action and crab predation ----

# Run the MANOVA test: 
MANOVA_LengthShape <- manova(cbind(RW1_Mean, ShellLength_Mean) ~ Wave_Mean + ScarsPresent_Frequency,
                             data = PopLev)
summary(MANOVA_LengthShape)

# Calculate the eta squared of the predictors:
etasq(MANOVA_LengthShape, test="Wilks")

#### 2.2.- Shape and size as adaptations to crabs ----

# Models: 
ScarsRW1_model <- lm(RW1_Mean ~ ScarsPresent_Frequency * Ria, data = PopLev) # Shape

ScarsShellLength_model <- lm(ShellLength_Mean ~ ScarsPresent_Frequency * Ria, data = PopLev) # Size

# Run diagnostics: 
ScarsRW1_res <- simulateResiduals(ScarsRW1_model) # Shape
plot(ScarsRW1_res)

ScarsShellLength_res <- simulateResiduals(ScarsShellLength_model) # Size
plot(ScarsShellLength_model)
# OUTPUT: The fits looks reasonably good

# See the summary of the models:
summary(ScarsRW1_model)

summary(ScarsShellLength_model)

#### 2.3.- Shape and size as adaptations to wave action ----

# Models: 
WaveRW1_model <- lm(RW1_Mean ~ Wave_Mean * Ria, data = PopLev) # Shape

WaveShellLength_model <- lm(ShellLength_Mean ~ Wave_Mean * Ria, data = PopLev) # Size

# Run diagnostics: 
WaveRW1_res <- simulateResiduals(WaveRW1_model) # Shape
plot(WaveRW1_res)

WaveShellLength_res <- simulateResiduals(WaveShellLength_model) # Size
plot(WaveShellLength_res)
# OUTPUT: The fits look reasonable

# See the summary of the model:
summary(WaveRW1_model)

summary(WaveShellLength_model)

#### 2.4.- Plots with linear models ----

# Plot of wave action and rw1 with linear model (Figure 6B top):
Fig_6B_t <- ggplot(PopLev, aes(x = Wave_Max, y = RW1_Mean)) +
  geom_point(aes(fill = Cluster), size = 3, alpha = 1,stroke = 1, show.legend = FALSE, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1), linewidth = 1) +
  scale_fill_manual(values = c("#63069D", "#C73890", "#EBD809")) +
  labs(title = "",
       x = "",
       y = "Relative warp 1") +
  scale_y_continuous(limits = c(-0.06, 0.06),
                     breaks = c("-0.05" = -0.05,"0" = 0, "0.05" = 0.05)) +
  scale_x_continuous(breaks = c("1" = 1, "3" = 3, "5" = 5)) +
  theme_cd2(18) +
  theme(axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    legend.text = element_blank(),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", linewidth = 1))

Fig_6B_t # View the plot 

# Plot of wave action and shell length with linear model (Figure 6B bottom):
Fig_6B_b <- ggplot(PopLev, aes(x = Wave_Max, y = ShellLength_Mean)) +
  geom_point(aes(fill = Cluster), size = 3, alpha = 1, stroke = 1, show.legend = FALSE, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1), linewidth = 1) +
  scale_fill_manual(values = c("#63069D", "#C73890", "#EBD809")) +
  labs(title = "",
       x = "Wave action",
       y = "Shell length (mm)") +
  scale_y_continuous(breaks = c("  5" = 5,"  7.5" = 7.5, "  10" = 10)) +
  scale_x_continuous(breaks = c("1" = 1, "3" = 3, "5" = 5)) +
  theme_cd2(18) +
  theme(panel.border = element_rect(colour = "black", linewidth = 1))

Fig_6B_b # View the plot 

# Plot of crab predation and rw1 with linear model (Figure 6C top):
Fig_6C_t <- ggplot(PopLev, aes(x = ScarsPresent_Frequency, y = RW1_Mean)) +
  geom_point(aes(fill = Cluster), size = 3, alpha = 1, stroke = 1, show.legend = FALSE, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1), linewidth = 1) +
  scale_fill_manual(values = c("#63069D", "#C73890", "#EBD809")) +
  labs(title = "",
       x = "",
       y = "Relative warp 1") +
  scale_y_continuous(limits = c(-0.06, 0.06),
                     breaks = c("-0.05" = -0.05,"0" = 0, "0.05" = 0.05)) +
  scale_x_continuous(breaks = c("0.1" = 0.1, "0.3" = 0.3, "0.5" = 0.5)) +
  theme_cd2(18) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        legend.text = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1)) 

Fig_6C_t # View the plot

# Plot of crab predation and shell length with linear model (Figure 6C bottom)
Fig_6C_b <- ggplot(PopLev, aes(x = ScarsPresent_Frequency, y = ShellLength_Mean)) +
  geom_point(aes(fill = Cluster), size = 3, alpha = 1, stroke = 1, show.legend = FALSE, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1), linewidth = 1) +
  scale_fill_manual(values = c("#63069D", "#C73890", "#EBD809")) +
  labs(title = "",
       x = "Crab predation",
       y = "Shell length (mm)") +
  scale_y_continuous(breaks = c("  5" = 5,"  7.5" = 7.5, "  10" = 10)) + # Remove the annoying leading 0s
  scale_x_continuous(breaks = c("0.1" = 0.1, "0.3" = 0.3, "0.5" = 0.5)) +
  theme_cd2(18) +
  theme(legend.text = element_blank(),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1))

Fig_6C_b # View the plot

### 3.- Composite figure (Figure 6) ----

# Set margins:
Margin <- theme(plot.margin = margin(2, 2, 2, 6, "mm"))

Fig_6A_T <- Fig_6A_T + Margin 
Fig_6A_B <- Fig_6A_B + Margin 

Fig_6B_t <- Fig_6B_t + Margin 
Fig_6B_b <- Fig_6B_b + Margin 

Fig_6C_t <- Fig_6C_t + Margin 
Fig_6C_b <- Fig_6C_b + Margin

# Add tags:
Fig_6A_T <- tags_cd(Fig_6A_T, "A", size = 22, x = 0.02, y = 0.98)
Fig_6B_t <- tags_cd(Fig_6B_t, "B", size = 22, x = 0.02, y = 0.98)
Fig_6C_t <- tags_cd(Fig_6C_t, "C", size = 22, x = 0.02, y = 0.98)

# Left column (Chromatic contrasts):
Fig_6_L <- plot_grid(
  Fig_6A_T, Fig_6A_B,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

# Center column (Shell length or Relative warp 1 ~ Wave action):
Fig_6_C <- plot_grid(
  Fig_6B_t, Fig_6B_b,
  ncol = 1,
  align = "v",  
  axis = "lr",
  rel_heights = c(0.8, 0.9)
)

# Right column (Shell length or Relative warp 1 ~ Crab predation):
Fig_6_R <- plot_grid(
  Fig_6C_t, Fig_6C_b,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(0.8, 0.9)
)

# Full figure:
Figure6 <- plot_grid(
  Fig_6_L, Fig_6_C, Fig_6_R,
  ncol = 3,
  align = "h",
  axis = "tb",
  rel_heights = c(1.0, 0.8, 0.8)
)

# Save in .svg to export to InkScape:
ggsave("Figure6.svg",
       Figure6,
       width = 330, height = 200, units = "mm")

## Supplementary analyses ----

### 1.- Spectra of color morphs and backgrounds ----

#### 1.1. Spectra of wave-sheltered ecosystems: Fulva morph, granite rocks, and Fucus ----

# Extract the names of the objects from the column names:
RefSpec_Col <- colnames(RefSpec_Zero) %>% # Extract the column names
  strsplit("_") %>% # Split the column names based on '_'
  sapply(function(x) x[2]) %>% # Select the fifth element of the column names (which contains the names of color morphs)
  .[-1] # Remove the first column, which corresponds to 'wl'

# Select the objects to plot:
RefSpec_Sheltered <- RefSpec_Zero %>%
  dplyr::select(wl, which(RefSpec_Col %in% c("Fulva", "Rock", "Fucus")) + 1)  # +1 because column 1 is 'wl'

# Now extract the names of the objects that you want to plot:
RefSpec_Sheltered_Names <- colnames(RefSpec_Sheltered) %>% 
  strsplit("_") %>%
  sapply(function(x) x[2]) %>% 
  .[-1] 

# Assign colors to the objects:
Colors_Sheltered <- c(
  Fulva = "#D6D2A9",
  Rock = "#856E42")

# Map the morphs in the reflectance spectra dataset to the colors:
Color_Unique <- unique(RefSpec_Sheltered_Names) #Set a unique color to each of the selected objects
Colors_Sheltered_ForPlot <- Colors_Sheltered[Color_Unique] # Map the colors 

# Plot with aesthetics:
aggplot(RefSpec_Sheltered,
                             by = RefSpec_Sheltered_Names,
                             FUN.center = mean, # We use the mean to calculate the line
                             FUN.error = function(x) 1.96 * (sd(x) / sqrt(length(x))), # We calculate the error with the 95% CI of the mean
                             legend = TRUE,
                             lcol = Colors_Sheltered_ForPlot, 
                             shadecol = Colors_Sheltered_ForPlot,
                             alpha = 0.6,
                             lty = c("solid"),
                             lwd = 1)

Fig_SM_SheltSpec # View the plot

#### 1.2. Spectra of wave-exposed ecosystems: Lineata morph, granite rock, and Lichina ----

# Extract the names of the objects from the column names:
RefSpec_Col <- colnames(RefSpec_Zero) %>% # Extract the column names
  strsplit("_") %>% # Split the column names based on '_'
  sapply(function(x) x[2]) %>% # Select the fifth element of the column names (which contains the names of color morphs)
  .[-1] # Remove the first column, which corresponds to 'wl'

# Select the objects to plot:
RefSpec_Exposed <- RefSpec_Zero %>%
  dplyr::select(wl, which(RefSpec_Col %in% c("Lineata", "Barnacles")) + 1)  # +1 because column 1 is 'wl'

# Now extract the names of the objects that you want to plot:
RefSpec_Exposed_Names <- colnames(RefSpec_Exposed) %>% 
  strsplit("_") %>%
  sapply(function(x) x[2]) %>% 
  .[-1] 

# Assign colors to the objects:
Colors_Exposed <- c(
  Lineata = "#A09D97",
  Barnacles = "black")

# Map the morphs in the reflectance spectra dataset to the colors:
Color_Unique <- unique(RefSpec_Exposed_Names) #Set a unique color to each of the selected objects
Colors_Exposed_ForPlot <- Colors_Exposed[Color_Unique] # Map the colors 

# Plot with aesthetics:
aggplot(RefSpec_Exposed,
                           by = RefSpec_Exposed_Names,
                           FUN.center = mean, # We use the mean to calculate the line
                           FUN.error = function(x) 1.96 * (sd(x) / sqrt(length(x))), # We calculate the error with the 95% CI of the mean
                           legend = TRUE,
                           #ylim = c(0, 25),
                           lcol = Colors_Exposed_ForPlot, 
                           shadecol = Colors_Exposed_ForPlot,
                           alpha = 0.6,
                           lty = c("solid"),
                           lwd = 1)

### 2.- Crab attack-induced plasticity ----

# GLMM to rule out the hypothesis of shape change as a plastic response to failed predation attempts:
ScarsPlasticity <- glmer(Scars ~ ShellThickness_Rel_Z + ShellLength_Z + RW1_Z + (1 | LocName),
                         data = IndLev,
                         family = binomial)
# Rationale: If predation attempts (i.e., scars) cause changes in shell shape, we should see differences in shape between individuals with and without scars in each locality

# Run diagnostics: 
ScarsPlasticity_SimulatedResiduals <- simulateResiduals(ScarsPlasticity) 
plot(ScarsPlasticity_SimulatedResiduals)

# Summary of the model:
summary(ScarsPlasticity)
# OUTPUT: RW1 is not significant, therefore the plasticity hypothesis is ruled out. Shell length is significant, but that might be due to larger snails having more time to be attacked
