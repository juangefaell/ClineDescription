############### ClineDescription - Covariation between Shell Traits #################

# BASIC INFO ----

# Researcher: [AuthorName] 
# Section aim: Is there covariation between color and other shell traits? Is such covariation equivalent across evolutionary units?
# Last update: 2026-05-19

# SETUP SECTION ----

## 1.- Packages ----
library(tidyverse) # Basic plots and statistical analyses
library(BSDA) # Vote counting
library(FactoMineR) # Dimension Reduction Analyses (PCA)
library(factoextra) # Dimension Reduction Analyses (PCA)
library(car) # Assumptions and ANOVA
library(effectsize) # Effect size of ANOVA (Eta-squared)
library(heplots) # Partial effect size
library(performance) # To calculate Tjur's R^2
library(DHARMa) # Model diagnostics
library(vegan) # Euclidean distances between morph centroids (PCA)
library(spdep) # Moran's I
library(ggrepel) # To repel labels in the partial effect size plot
library(cowplot) # Composite figures
library(grid) # Composite figures 

## 2.- Dataset loading ----
setwd("") # Set the directory of the dataset file
read.csv("CD_DA_RD_IndividualLevel.csv") 
read.csv("CD_DA_RD_PopulationLevel.csv")

## 3.- Data handling and overview ----

# IndividualLevel dataset:
IndLev <- read.csv("CD_DA_RD_IndividualLevel.csv")

  # Create the ShellThickness_Rel (relativized to shell size) to account for allometric effects: 
    # ST_model <- lm(log(ShellThickness) ~ log(ShellLength), data = IndLev)
  
  # Incorporate the residuals of the model as "Shell thickness relativized" to the dataset:
    # IndLev <- IndLev %>% 
      # add_residuals(ST_model, var = "ShellThickness_Rel")

# Format the variables:
IndLev <- IndLev %>%
  mutate(
    Ria = factor(Ria, levels = c("MurosNoia", "Arousa", "Pontevedra", "Vigo")),
    RiaCoast = factor(Ria, levels = c("MurosNoiaNorth", "MurosNoiaSouth", "ArousaNorth", "ArousaSouth", "PontevedraNorth", "PontevedraSouth", "VigoNorth", "VigoSouth")),
    Coast = factor(Coast, levels = c("South", "North")),
    Cluster = factor(Cluster, levels = c("WaveSheltered", "WaveIntermediate", "WaveExposed"), ordered = TRUE),
    Lineata = if_else(ColorMorph == "Lineata", "1", "0"), # Create a Lineata variable
    Fulva = if_else(ColorMorph == "Fulva", "1", "0"), # Create a Fulva variable
    SculpturedShell = if_else(ShellSculpture == "Jugosa", "1", "0"), # Create the SculpturedShell variable (Sculptured snail or not)
    RW1_Z = scale(RW1), # Scale the RW1 variable
    ShellLength_Z = scale(ShellLength), # Scale the ShellLength variable
    ShellThickness_Rel_Z = scale(ShellThickness_Rel)) # Scale the ShellThickness_Rel variable
str(IndLev) # Check the structure of the data

# PopulationLevel dataset: 
PopLev <- read.csv("CD_DA_RD_PopulationLevel.csv")

# Create the ShellThickness_Rel_Mean variable and incorporate it on the PopLev dataset:
  # ShellThickness_Rel_MeanSD <- read.csv("CD_DA_RD_IndividualLevel.csv") %>% # Create a new object
  # group_by(LocNumCumulative) %>% # Group by locality 
  # summarise(ShellThickness_Rel_Mean = mean(ShellThickness_Rel, na.rm = TRUE),
  # ShellThickness_Rel_SD = sd(ShellThickness_Rel, na.rm = TRUE))
  # write.csv(ShellThickness_Rel_MeanSD, "ShellThickness_Rel_MeanSD_ForPopulationLevel.csv", row.names = FALSE) # Save the new variable to a CSV file (then incorporate it into the PopLev dataset)

# Format the variables:
PopLev <- PopLev %>%
  mutate(
    Ria = factor(Ria, levels = c("MurosNoia", "Arousa","Pontevedra", "Vigo"))) # Set an order from north to south
str(PopLev) # Check the structure of the data

## 4.- Themes for figures ----

# Theme for figures, v. 1:
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

### 1.- Magnitude of covariation between color and other shell traits by color morph ----

# Subsection aim: Reduce the dimensionality of the shell traits and see whether they are different across color morphs in the morphospace

#///

# Transform the SculpturedShell variable to numeric format to be able to include it on the PCA:
IndLev <- IndLev %>%
  mutate(
    SculpturedShell = as.numeric(SculpturedShell))

# Select only the Fulva and Lineata morphs:
IndLev_LinFul <- IndLev %>%
  filter(ColorMorph %in% c("Fulva", "Lineata")) %>% 
  dplyr::select(ColorMorph, LocName, RiaCoast, RW1, ShellLength, ShellThickness_Rel, SculpturedShell) %>% 
  drop_na()

# Selecting the variables:
ShellPCA_Vars <- IndLev_LinFul[, c("ShellLength", "RW1", 
                                   "ShellThickness_Rel", "SculpturedShell")] 
# Check the format of the variables (IMPORTANT!): 
str(ShellPCA_Vars) 

# Perform PCA:
ShellPCA_Result <- prcomp(ShellPCA_Vars, center = TRUE, scale = TRUE) # Perform PCA, with variables scaled
summary(ShellPCA_Result) # Display the % of variance explain by each PC
# OUTCOME: PC1 explains 37.87% of variance; PC2: 27.66%; PC3: 21.27%

# Display the SD and component loadings of the individual variables:
print(ShellPCA_Result) 

# Visualization of PCA results (guide: https://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining):
fviz_pca_biplot(ShellPCA_Result, label = "var", addEllipses = TRUE, ellipse.level = 0.95, repel = TRUE) # Biplot of variables to visualize the directions of variables (individuals can be added)

# Extract the scores:
ShellPCA_Scores <- as.data.frame(ShellPCA_Result$x)

# Set ColorMorph as grouping factor:
ShellPCA_Scores$ColorMorph <- IndLev_LinFul$ColorMorph

# Assess differences in centroids between morphs:
EuclideanPCs_Morphs <- adonis2(ShellPCA_Scores[, c("PC1", "PC2")] ~ ColorMorph, data = ShellPCA_Scores, method = "euclidean")

EuclideanPCs_Morphs # See the results
# OUTPUT: The color morph centroids are different, with morph explaining 15.63% of variation in PC1 and PC2 scores

# Checking dispersion:
anova(betadisper(dist(ShellPCA_Scores[, c("PC1", "PC2")]), ShellPCA_Scores$ColorMorph, type = "centroid"))
# OUTPUT: Dispersion exists, interpret cautiously

# Aesthetics for ColorMorph:
Color_FulvaLineata <- c(
  Fulva = "#D6D2A9",
  Lineata = "#A09D97"
)

# Extract the centroids:
ShellPCA_Centroids <- ShellPCA_Scores %>% 
  group_by(ColorMorph) %>% 
  summarise(across(c(PC1, PC2), mean))

# Plot of PCA centroids:
Fig_5A <- ggplot(ShellPCA_Scores, aes(x = PC1, y = PC2, fill = ColorMorph)) +
  geom_point(size = 3, shape = 21, color = "black", alpha = 0.8) +
  stat_ellipse(type = "t", level = 0.95, alpha = 0.3, geom = "polygon", color = "black") +
  geom_point(data = ShellPCA_Centroids, size = 4, shape = 23, stroke = 1) +
  scale_fill_manual(values = Color_FulvaLineata) +
  theme_cd(18) +
  labs(
    x = "PC1 (37.9%)",
    y = "PC2 (27.7%)") +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c("-2.5" = -2.5, "0" = 0, "2.5" = 2.5, "5" = 5))

Fig_5A # View the plot

### 2.- Differences in covariation across coasts (Based on James et al. 2021: https://github.com/OB-lab/James_et_al._2021-Evolution/blob/main/phenotype/R_code/phenotype.R) ----

# Subsection aim: Is shell color covariation equivalent in direction across coasts (i.e., lineages)?

# /// 

#### 2.1.- Direction of variation plots ----

# Filter the Fulva and Lineata cases:
IndLev_LinFul <- IndLev %>% 
  filter(ColorMorph %in% c("Fulva", "Lineata")) # Filter only the localities with Wave-sheltered and Wave-exposed ecosystems

# Calculate means of quantitative traits by morph and coast: 
Coast_Means <- IndLev_LinFul %>% 
  group_by(RiaCoast, ColorMorph) %>% 
  summarise(Mean_ShellLength = mean(ShellLength, na.rm = TRUE), # Do not z-scale to retain original units
            Mean_ShellThickness = mean(ShellThickness_Rel_Z, na.rm = TRUE),
            Mean_RelativeWarp1 = mean(RW1_Z, na.rm = TRUE))

# Calculate the proportion of sculptured individuals per morph and coast:
Coast_Prop <- IndLev_LinFul %>% 
  group_by(RiaCoast, ColorMorph) %>% 
  summarise(
    n = n(),
    n_Sculptured = sum(SculpturedShell == 1, na.rm = TRUE),
    Proportion_Sculptured = n_Sculptured / n)

Color_FulvaLineata <- c(
  Fulva = "#D6D2A9",
  Lineata = "#A09D97"
)

# Shell length plot with aesthetics (Figure 5B, part 1): 
Fig_5B_1 <- ggplot(Coast_Means, aes(x = ColorMorph, y = Mean_ShellLength, group = RiaCoast)) +
  geom_line(linewidth = 0.5, color = "black") +
  geom_point(aes(fill = ColorMorph), color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = Color_FulvaLineata) +
  labs(x = "",
       y = "Shell length (mm)",
       color = "") +
  theme_cd(18) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic"))

Fig_5B_1 # View the plot

# Shell thickness plot with aesthetics (Figure 5B, part 2):
Fig_5B_2 <- ggplot(Coast_Means, aes(x = ColorMorph, y = Mean_ShellThickness, group = RiaCoast)) +
  geom_line(linewidth = 0.5, color = "black") +
  geom_point(aes(fill = ColorMorph), color = "black", size = 4, shape = 21) +
  # geom_text(aes(label = RiaCoast), hjust = -0.1, vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = Color_FulvaLineata) +
  labs(x = "",
       y = "Shell thickness (rel.)",
       color = "") +
  theme_cd(18) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic")) +
  scale_y_continuous(breaks = c("-0.5" = -0.5, "0" = 0, "0.5" = 0.5)) # Remove the annoying leading 0s

Fig_5B_2 # View the plot

# Relative warp 1 plot with aesthetics (Figure 5B, part 3):
Fig_5B_3 <- ggplot(Coast_Means, aes(x = ColorMorph, y = Mean_RelativeWarp1, group = RiaCoast)) +
  geom_line(linewidth = 0.5, color = "black") +
  geom_point(aes(fill = ColorMorph), color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = Color_FulvaLineata) +
  labs(x = "",
       y = "Relative warp 1",
       color = "") +
  theme_cd(18) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic")) +
  scale_y_continuous(breaks = c("-0.5" = -0.5, "-0.25" = -0.25, "0" = 0, "0.25" = 0.25, "0.5" = 0.5)) # Remove the annoying leading 0s

Fig_5B_3 # View the plot

# Shell sculpture (proportions) plot with aesthetics (Figure 5B, part 4):
Fig_5B_4 <- ggplot(Coast_Prop, aes(x = ColorMorph, y = Proportion_Sculptured, group = RiaCoast)) +
  geom_line(linewidth = 0.5, color = "black") +
  geom_point(aes(fill = ColorMorph), color = "black", size = 4, shape = 21) +
  scale_fill_manual(values = Color_FulvaLineata) +
  labs(x = "",
       y = "Sculptured shells (prop.)",
       color = "") +
  theme_cd(18) +
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic")) +
  scale_y_continuous(breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75))

Fig_5B_4 # View the plot

#### 2.2.- Vote counting ----

# Shell length: 
Coast_SL_Wide <- Coast_Means %>% # Turn into a wide format
  select(RiaCoast, ColorMorph, Mean_ShellLength) %>% 
  pivot_wider(names_from = ColorMorph, values_from = Mean_ShellLength)

VoteCounting_SL <- SIGN.test( # Calculate the vote counting
  x = Coast_SL_Wide$Fulva,
  y = Coast_SL_Wide$Lineata,
  alternative = "two.sided"
) 
VoteCounting_SL$p.value 
# OUTPUT: 0.0078125

# Shell thickness: 
Coast_ST_Wide <- Coast_Means %>% # Turn into a wide format
  select(RiaCoast, ColorMorph, Mean_ShellThickness) %>% 
  pivot_wider(names_from = ColorMorph, values_from = Mean_ShellThickness)

VoteCounting_ST <- SIGN.test( # Calculate the vote counting
  x = Coast_ST_Wide$Fulva,
  y = Coast_ST_Wide$Lineata,
  alternative = "two.sided"
) 
VoteCounting_ST$p.value 
# OUTPUT: 0.2890625

# Relative warp 1: 
Coast_RW1_Wide <- Coast_Means %>% # Turn into a wide format
  select(RiaCoast, ColorMorph, Mean_RelativeWarp1) %>% 
  pivot_wider(names_from = ColorMorph, values_from = Mean_RelativeWarp1)

VoteCounting_RW1 <- SIGN.test( # Calculate the vote counting
  x = Coast_RW1_Wide$Fulva,
  y = Coast_RW1_Wide$Lineata,
  alternative = "two.sided"
) 
VoteCounting_RW1$p.value 
# OUTPUT: 0.0078125

# Shell sculpture: 
Coast_SS_Wide <- Coast_Prop %>% # Turn into a wide format
  select(RiaCoast, ColorMorph, Proportion_Sculptured) %>% 
  pivot_wider(names_from = ColorMorph, values_from = Proportion_Sculptured)

VoteCounting_SS <- SIGN.test( # Calculate the vote counting
  x = Coast_SS_Wide$Fulva,
  y = Coast_SS_Wide$Lineata,
  alternative = "two.sided"
) 
VoteCounting_SS$p.value 
# OUTPUT: 0.0078125

### 3.- Relative importance of lineage and repeated selection (Based on James et al. 2021: https://github.com/OB-lab/James_et_al._2021-Evolution/blob/main/phenotype/R_code/phenotype.R) ----

# Subsection aim: What is more important for covariation, lineage or selection?

#/// 

# Prepare the data (selecting only Fulva and Lineata and the variables of interest):
IndLev_MANOVA <- IndLev %>%
  filter(ColorMorph %in% c("Fulva", "Lineata")) %>% 
  dplyr::select(ColorMorph, LocName, RiaCoast, RW1_Z, ShellLength_Z, ShellThickness_Rel_Z, SculpturedShell) %>% 
  drop_na()

# Convert "SculpturedShell" into a factor:
IndLev_MANOVA <- IndLev_MANOVA %>% 
  mutate(SculpturedShell = as.factor(SculpturedShell))

# Set a reference level: 
IndLev_MANOVA <- IndLev_MANOVA %>% 
  mutate(ColorMorph = as.factor(ColorMorph)) %>% 
  mutate(ColorMorph = relevel(ColorMorph, ref = "Fulva"))

# Check the dataset:
str(IndLev_MANOVA)
head(IndLev_MANOVA)

# Run the MANOVA test: 
MANOVA_Traits_RiaMorph <- manova(cbind(RW1_Z, ShellLength_Z, ShellThickness_Rel_Z) ~ ColorMorph * RiaCoast,
                                 data = IndLev_MANOVA)
summary(MANOVA_Traits_RiaMorph)

# Calculate the eta squared of the predictors:
etasq(MANOVA_Traits_RiaMorph, test = "Wilks")

# Linear models for individual traits (Shell length):
ShellLength_lm <- lm(ShellLength_Z ~ ColorMorph * RiaCoast, data = IndLev_MANOVA)
summary(ShellLength_lm)

etasq(ShellLength_lm) # Effect size

## Linear models for individual traits (Shell thickness):
ShellThickness_lm <- lm(ShellThickness_Rel_Z ~ ColorMorph * RiaCoast, data = IndLev_MANOVA)
summary(ShellThickness_lm)

etasq(ShellThickness_lm)

# Linear models for individual traits (RW1):
RW1_lm <- lm(RW1_Z ~ ColorMorph * RiaCoast, data = IndLev_MANOVA)
summary(RW1_lm)

etasq(RW1_lm) # Effect size

# Logistic model for shell sculpture:
SculpturedShell_glm <- glm(SculpturedShell ~ ColorMorph * RiaCoast, family = binomial(link = "logit"), data = IndLev_MANOVA)
summary(SculpturedShell_glm)

# Run diagnostics: 
SimulatedResiduals_SculpturedShell_glm <- simulateResiduals(SculpturedShell_glm)
plot(SimulatedResiduals_SculpturedShell_glm)
# OUTPUT: Very good fit

# Run models without ColorMorph and the interaction:
SS_WTColor <- glm(SculpturedShell ~ RiaCoast, family = binomial, data = IndLev_MANOVA) # Logistic model without ColorMorph
SS_WTInteraction <- glm(SculpturedShell ~ ColorMorph + RiaCoast, family = binomial, data = IndLev_MANOVA) # Logistic model without the interaction

# Calculate Tjur's R2s of each model:
SS_R2_fullmodel <- r2_tjur(SculpturedShell_glm) 
SS_R2_WTColor <- r2_tjur(SS_WTColor)
SS_R2_WTInteraction <- r2_tjur(SS_WTInteraction)

# Get the partial effect sizes subtracting Tjur's R2s values:
partial_R2Tjur_ColorMorph <- SS_R2_fullmodel - SS_R2_WTColor # OUTPUT for ColorMorph: 0.1229063
partial_R2Tjur_Interaction <- SS_R2_fullmodel - SS_R2_WTInteraction # OUTPUT for interaction: 0.005088147

## Create variables for a data frame of partial eta-squares:
Trait <- c("All traits", "Shell length", "Shell thickness", "Relative warp 1", "Shell sculpture")
Morph_EtaSq <- c("0.28649286", "0.22806369", "0.03348837", "0.04977246", "0.1229063")
MorphCoast_EtaSq <- c("0.03970950", "0.08486516", "0.02274045", "0.01123074", "0.005088147")

# Create a data frame with the values: 
EtaSq_df <- data.frame(Trait, Morph_EtaSq, MorphCoast_EtaSq)

# Set the correct variable format:
EtaSq_df <- EtaSq_df %>% 
  mutate(
    Morph_EtaSq = as.numeric(Morph_EtaSq),
    MorphCoast_EtaSq = as.numeric(MorphCoast_EtaSq)
  )

# Plot of partial effect sizes (Figure 5C)
Fig_5C <- ggplot(EtaSq_df, aes(x = MorphCoast_EtaSq, y = Morph_EtaSq)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linewidth = 1, linetype = "dashed") +
  geom_point(size = 5, alpha = 1, shape = 21, fill = "#69b3a2", color = "black") +
  geom_point(data = subset(EtaSq_df, Trait == "All traits"), size = 5, shape = 21, fill = "#CC1C00", color = "black") +
  geom_text_repel(aes(label = Trait), size = 6, hjust = -0.2, vjust = -0.1) +
  labs(title = "",
       x = "Morph * Coast",
       y = "Morph") +
  theme_cd(18) +
  scale_y_continuous(limits = c(0, 0.4),
    breaks = c("0" = 0.0, "0.1" = 0.1, "0.2" = 0.2, "0.3" = 0.3, "0.4" = 0.4))+
  scale_x_continuous(limits = c(0, 0.4),
                     breaks = c("0" = 0.0, "0.1" = 0.1, "0.2" = 0.2, "0.3" = 0.3, "0.4" = 0.4))

Fig_5C # View the plot

### 4.- Phenotypic Change Vector Analysis (PCVA) ----

# Subsection aim: Test if, when all traits are considered, the direction and magnitude of phenotypic change between Fulva and Lineata is equivalent across coasts (i.e., lineages)

# ///

#### 4.1.- Preparing the data ----

# Transform the SculpturedShell variable to be able to include it on the PCA:
IndLev <- IndLev %>%
  mutate(
    SculpturedShell = as.numeric(SculpturedShell))

# Select only the Fulva and Lineata morphs:
IndLev_LinFul <- IndLev %>%
  filter(ColorMorph %in% c("Fulva", "Lineata")) %>% 
  dplyr::select(ColorMorph, LocName, RiaCoast, RW1, ShellLength, ShellThickness_Rel, SculpturedShell) %>% 
  drop_na()

# Selecting the variables:
ShellPCA_Vars <- IndLev_LinFul[, c("ShellLength", "RW1", 
                                   "ShellThickness_Rel", "SculpturedShell")] 

# Check the format of the variables (IMPORTANT!): 
str(ShellPCA_Vars) 

# Perform PCA:
ShellPCA_Result <- prcomp(ShellPCA_Vars, center = TRUE, scale = TRUE) # Perform PCA, with variables scaled
summary(ShellPCA_Result) # Display the % of variance explain by each PC
print(ShellPCA_Result) # Display the SD and component loadings of the individual variables
# OUTCOME: PC1 explains 37.87% of variance; PC2: 27.66%; PC3: 21.27%

# Extract the scores:
ShellPCA_Scores <- as.data.frame(ShellPCA_Result$x)

# Add the ColorMorph variable:
ShellPCA_Scores$ColorMorph <- IndLev_LinFul$ColorMorph

# Add the RiaCoast variable:
ShellPCA_Scores$RiaCoast <- IndLev_LinFul$RiaCoast

# Extract the centroids:
ShellPCA_Centroids_Coast <- ShellPCA_Scores %>% 
  group_by(RiaCoast, ColorMorph) %>% 
  summarise(
    Mean_PC1 = mean(PC1),
    Mean_PC2 = mean(PC2),
    .groups = "drop")

#### 4.2.- Vectors of phenotypic change ----

# Create centroid pairs:
CentroidPairs <- ShellPCA_Centroids_Coast %>% 
  pivot_wider(
    names_from = ColorMorph,
    values_from = c(Mean_PC1, Mean_PC2))

# Compute distances to get length and angles of vectors:
CentroidVectors <- CentroidPairs %>% 
  mutate(
    dPC1 = Mean_PC1_Lineata - Mean_PC1_Fulva,
    dPC2 = Mean_PC2_Lineata - Mean_PC2_Fulva)

# Add the lengths:
CentroidVectors <- CentroidVectors %>% 
  mutate(Length = sqrt(dPC1^2 + dPC2^2))

# Add the angles:
CentroidVectors <- CentroidVectors %>% 
  mutate(
    Angle_rad = atan2(dPC2, dPC1),
    Angle_deg = Angle_rad * 180 / pi)

# Check if it worked:  
head(CentroidVectors)

# Assign colors:
Colors_PhenoTraj <- c(Fulva = "#D6D2A9",
                      Lineata = "#A09D97")

# Plot of phenotypic trajectories (figure 5D): 
Fig_5D <- ggplot(ShellPCA_Centroids_Coast, aes(x = Mean_PC1, y = Mean_PC2)) +
  geom_segment(data = CentroidVectors, aes(x = Mean_PC1_Fulva, y = Mean_PC2_Fulva, 
                                           xend = Mean_PC1_Lineata, yend = Mean_PC2_Lineata), color = "black", linewidth = 0.5) +
  labs(x = "PC1 (37.9%)",
       y = "PC2 (27.7%)") +
  geom_point(aes(fill = ColorMorph), size = 4, shape = 21) +
  #geom_text(aes(label = RiaCoast), hjust = -0.1, vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = Colors_PhenoTraj) +
  theme_cd(18) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c("-0.75" = -0.75, "-0.5" = -0.5, "-0.25" = -0.25, "0" = 0, "0.25" = 0.25, "0.5" = 0.5)) # Remove the annoying leading 0s

Fig_5D # View the plot 

#### 4.3.- Magnitude of change across replicates ----

# Calculate Euclidean distances in vector length (i.e., magnitude of change):
Dist_length <- dist(CentroidVectors$Length, method = "euclidean")
Dist_length <- as.numeric(Dist_length)
Dist_length_df <- data.frame(Distance = Dist_length)

# Extract the mean and sd:
median(Dist_length_df$Distance)
IQR(Dist_length_df$Distance)
# OUTPUT: Median = 0.5405983, IQR = 0.5531487

# Plot histogram of differences in length (Figure 5E, part 1):
Fig_5E_1 <- ggplot(Dist_length_df, aes(x = Distance)) + # Remove the NAs
  geom_histogram(binwidth = 0.2, fill = "#69b3a2", color = "black", alpha = 0.7) +
  labs(
    x = "Pairwise euclidean  \n distances in vector length",
    y = "Count") + 
  theme_cd(18) +
  theme(plot.margin = margin(4, 4, 2, 4, "mm")) +
  scale_x_continuous(breaks = c("0" = 0, "0.5" = 0.5, "1" = 1, "1.5" = 1.5)) # Remove the annoying leading 0s

Fig_5E_1 # View the plot

#### 4.4.- Contribution of traits to divergence ----

# Extract the angles (in degrees):
Angles <- CentroidVectors$Angle_deg

# Calculate pairwise absolute angle differences:
Angle_Dif_Matrix <- outer(Angles, Angles, function(a, b) abs(a - b))

# Convert to circular angle distances:
Angle_Dif_Matrix <- pmin(Angle_Dif_Matrix, 360 - Angle_Dif_Matrix)

# Convert to a "Dist" object:
Dist_Angle <- as.dist(Angle_Dif_Matrix)

# Create a numeric object:
Angle_Differences <- as.numeric(Dist_Angle)

# Convert to a dataframe:
Angle_Differences_df <- data.frame(Angle_Dif = Angle_Differences)

# Extract the mean and sd:
median(Angle_Differences_df$Angle_Dif)
IQR(Angle_Differences_df$Angle_Dif)
# OUTPUT: Median = 14.03926, IQR = 14.57729

# Plot histogram of angle differences (Figure 5E, part 2):
Fig_5E_2 <- ggplot(data.frame(Angle_Diff = Angle_Differences), aes(x = Angle_Diff)) + # Remove the NAs
  geom_histogram(binwidth = 5, fill = "#69b3a2", color = "black", alpha = 0.7) +
  labs(
    x = "Pairwise differences in angles",
    y = "Count") + 
  theme_cd(18) +
  theme(plot.margin = margin(2, 4, 4, 4, "mm"))

Fig_5E_2 # View the plot

### 5.- Composite figure (Figure 5) ----

# Add tags: 
Fig_5A_Tag <- tags_cd(Fig_5A, "A", size = 22, x = 0.02, y = 0.98)
Fig_5B_1_Tag <- tags_cd(Fig_5B_1, "B", size = 22, x = 0.02, y = 0.98)
Fig_5C_Tag <- tags_cd(Fig_5C, "C", size = 22, x = 0.02, y = 0.98)
Fig_5D_Tag <- tags_cd(Fig_5D, "D", size = 22, x = 0.02, y = 0.98)
Fig_5E_1 <- tags_cd(Fig_5E_1, "", size = 22, x = 0.02, y = 0.98)
Fig_5E_2 <- tags_cd(Fig_5E_2, "", size = 22, x = 0.02, y = 0.98)

# Create top row:
Fig_5_Top <- plot_grid(Fig_5A_Tag, Fig_5B_1_Tag, Fig_5B_2, Fig_5B_3, Fig_5B_4,
                       ncol = 5,
                       align = "h",
                       axis = "tb",
                       rel_widths = c(2, 0.9, 1, 1.1, 1))

# Create bottom row:
Fig_5_Bot <- ggdraw() +
  draw_plot(Fig_5C_Tag, x = 0.00, y = 0.00, width = 0.33, height = 1.00) +  # Panel C
  draw_plot(Fig_5D_Tag, x = 0.33, y = 0.00, width = 0.33, height = 1.00) +  # Panel D
  draw_plot(Fig_5E_1, x = 0.69, y = 0.51, width = 0.31, height = 0.51) +  # Panel E top histogram
  draw_plot(Fig_5E_2, x = 0.69, y = 0.00, width = 0.31, height = 0.43) +  # Panel E bottom histogram
  draw_plot_label("E", x = 0.68, y = 0.98, hjust = 0, vjust = 1,  # Label for whole E block
                  size = 22, fontface = "bold")

# Entire figure:
Figure5 <- plot_grid(Fig_5_Top, Fig_5_Bot,
                   ncol = 1, 
                   align = "v",
                   axis = "lr",
                   rel_heights = c(1.0, 1.0))

# Save in .svg to export to InkScape:
ggsave("Figure5.svg",
       Figure5,
       width = 430, height = 210, units = "mm")


## Supplementary analyses ----

### 1.- Visualizing the allometric relationship between shell thickness and shell length ----
Fig_SM_ShellThick <- ggplot(IndLev, aes(x = ShellLength, y = ShellThickness)) +
  geom_point(fill = "grey", size = 2, alpha = 1, stroke = 0.5, show.legend = TRUE, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1, color = "black") +
  labs(title = "",
       x = "Shell length (mm)",
       y = "Shell thickness (mm)") +
  theme_cd(18) +
  scale_x_continuous(breaks = c("5" = 5, "7.5" = 7.5, "10" = 10, "12.5" = 12.5)) # Remove the annoying leading 0s

Fig_SM_ShellThick # View the plot

# Save the plot:
ggsave("CD_MS_SM_Fig_S2.svg",
       plot = Fig_SM_ShellThick,
       width = 6,
       height = 5)

### 2.- Strength of covariation and ecosystem clustering ---- 

#### 2.1.- Population-level covariation strength ----

#####  2.1.1.- Fulva ----

###### 2.1.1.1.- Fulva in Muros-Noia ----

# Fix the "MurosNoia" name by adding a "-" between them:
PopLev <- PopLev %>% 
  mutate(Ria = fct_recode(Ria, 
                          "Muros-Noia" = "MurosNoia"))

# Filter the ria:
PopLev_MN <- PopLev %>% 
  filter(Ria %in% c("Muros-Noia"))

# Run the models:
Fulva_MN_full <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                           ShellLength_Mean + RW1_Mean + Jugosa_Frequency + ShellThickness_Rel_Mean,
                         family = betabinomial(link = "logit"),
                         data = PopLev_MN)

Fulva_MN_null <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                           1,
                         family = betabinomial(link = "logit"),
                         data = PopLev_MN)

# Run diagnostics: 
SimulatedResiduals_Fulva_MN <- simulateResiduals(Fulva_MN_full)
plot(SimulatedResiduals_Fulva_MN)

# McFadden R^2:
1 - (logLik(Fulva_MN_full) / logLik(Fulva_MN_null))
# OUTPUT: 'log Lik.' 0.2362244 (df=6)

###### 2.1.1.2.- Fulva in Arousa ----

# Filter the ria:
PopLev_A <- PopLev %>% 
  filter(Ria %in% c("Arousa"))

# Run the models:
Fulva_A_full <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                          ShellLength_Mean + RW1_Mean + Jugosa_Frequency + ShellThickness_Rel_Mean,
                        family = betabinomial(link = "logit"),
                        data = PopLev_A)

Fulva_A_null <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                          1,
                        family = betabinomial(link = "logit"),
                        data = PopLev_A)

# Run diagnostics: 
SimulatedResiduals_Fulva_A <- simulateResiduals(Fulva_A_full)
plot(SimulatedResiduals_Fulva_A)

# McFadden R^2:
1 - (logLik(Fulva_A_full) / logLik(Fulva_A_null))
# OUTPUT: 'log Lik.' 0.09004023 (df=6)


###### 2.1.1.3.- Fulva in Pontevedra ----

# Filter the ria:
PopLev_P <- PopLev %>% 
  filter(Ria %in% c("Pontevedra"))

# Run the models:
Fulva_P_full <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                          ShellLength_Mean + RW1_Mean + Jugosa_Frequency + ShellThickness_Rel_Mean,
                        family = betabinomial(link = "logit"),
                        data = PopLev_P)

Fulva_P_null <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                          1,
                        family = betabinomial(link = "logit"),
                        data = PopLev_P)

# Run diagnostics: 
SimulatedResiduals_Fulva_P <- simulateResiduals(Fulva_P_full)
plot(SimulatedResiduals_Fulva_P)

# McFadden R^2:
1 - (logLik(Fulva_P_full) / logLik(Fulva_P_null))
# OUTPUT: 'log Lik.' 0.2278874 (df=6)

###### 2.1.1.4.- Fulva in Vigo ----

# Filter the ria:
PopLev_V <- PopLev %>% 
  filter(Ria %in% c("Vigo"))

# Run the models:
Fulva_V_full <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                          ShellLength_Mean + RW1_Mean + Jugosa_Frequency + ShellThickness_Rel_Mean,
                        family = betabinomial(link = "logit"),
                        data = PopLev_V)

Fulva_V_null <- glmmTMB(cbind(Fulva_Counts, n - Fulva_Counts) ~
                          1,
                        family = betabinomial(link = "logit"),
                        data = PopLev_V)

# Run diagnostics: 
SimulatedResiduals_Fulva_V <- simulateResiduals(Fulva_V_full)
plot(SimulatedResiduals_Fulva_V)

# McFadden R^2:
1 - (logLik(Fulva_V_full) / logLik(Fulva_V_null))
# OUTPUT: 'log Lik.' 0.4232778 (df=6)


#####  2.1.2.- Lineata ----

###### 2.1.2.1.- Lineata in Muros-Noia ----

# Fix the "MurosNoia" name by adding a "-" between them:
PopLev <- PopLev %>% 
  mutate(Ria = fct_recode(Ria, 
                          "Muros-Noia" = "MurosNoia"))

# Filter the ria:
PopLev_MN <- PopLev %>% 
  filter(Ria %in% c("Muros-Noia"))

# Run the models:
Lineata_MN_full <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                             ShellLength_Mean + RW1_Mean + Jugosa_Frequency + ShellThickness_Rel_Mean,
                           family = betabinomial(link = "logit"),
                           data = PopLev_MN)

Lineata_MN_null <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                             1,
                           family = betabinomial(link = "logit"),
                           data = PopLev_MN)

# Run diagnostics: 
SimulatedResiduals_Lineata_MN <- simulateResiduals(Lineata_MN_full)
plot(SimulatedResiduals_Lineata_MN)

# McFadden R^2:
1 - (logLik(Lineata_MN_full) / logLik(Lineata_MN_null))
# OUTPUT: 'log Lik.' 0.2996492 (df=6)

###### 2.1.2.2.- Lineata in Arousa ----

# Filter the ria:
PopLev_A <- PopLev %>% 
  filter(Ria %in% c("Arousa"))

# Run the models:
Lineata_A_full <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                            ShellLength_Mean + RW1_Mean + Jugosa_Frequency + ShellThickness_Rel_Mean,
                          family = betabinomial(link = "logit"),
                          data = PopLev_A)

Lineata_A_null <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                            1,
                          family = betabinomial(link = "logit"),
                          data = PopLev_A)

# Run diagnostics: 
SimulatedResiduals_Lineata_A <- simulateResiduals(Lineata_A_full)
plot(SimulatedResiduals_Lineata_A)

# McFadden R^2:
1 - (logLik(Lineata_A_full) / logLik(Lineata_A_null))
# OUTPUT: log Lik.' 0.1598871 (df=6)

###### 2.1.2.3.- Lineata in Pontevedra ----

# Filter the ria:
PopLev_P <- PopLev %>% 
  filter(Ria %in% c("Pontevedra"))

# Run the models:
Lineata_P_full <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                            ShellLength_Mean + RW1_Mean + Jugosa_Frequency + ShellThickness_Rel_Mean,
                          family = betabinomial(link = "logit"),
                          data = PopLev_P)

Lineata_P_null <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                            1,
                          family = betabinomial(link = "logit"),
                          data = PopLev_P)

# Run diagnostics: 
SimulatedResiduals_Lineata_P <- simulateResiduals(Lineata_P_full)
plot(SimulatedResiduals_Lineata_P)

# McFadden R^2:
1 - (logLik(Lineata_P_full) / logLik(Lineata_P_null))
# OUTPUT: 'log Lik.' 0.1042499 (df=6)

###### 2.1.2.4.- Lineata in Vigo ----

# Filter the ria:
PopLev_V <- PopLev %>% 
  filter(Ria %in% c("Vigo"))

# Run the models:
Lineata_V_full <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                            ShellLength_Mean + RW1_Mean + Jugosa_Frequency + ShellThickness_Rel_Mean,
                          family = betabinomial(link = "logit"),
                          data = PopLev_V)

Lineata_V_null <- glmmTMB(cbind(Lineata_Counts, n - Lineata_Counts) ~
                            1,
                          family = betabinomial(link = "logit"),
                          data = PopLev_V)

# Run diagnostics: 
SimulatedResiduals_Lineata_V <- simulateResiduals(Lineata_V_full)
plot(SimulatedResiduals_Lineata_V)

# McFadden R^2:
1 - (logLik(Lineata_V_full) / logLik(Lineata_V_null))
# OUTPUT: 'log Lik.' 0.2683315 (df=6)

#### 2.2.- Moran's I (To asses clustering) ----



##### 2.2.1.- Muros-Noia ----

# Select the ria:
PopLev_MN <- subset(PopLev, Ria == "MurosNoia") 

# Extract the coordinates:
Coordinates_MN <- cbind(PopLev_MN$Latitude, PopLev_MN$Longitude) 

# Calculate spatial weights matrix:
knn_MN <- knearneigh(Coordinates_MN, k = 4) # Using K-nearest neighbors, k = 4; this was used as a very liberal account of dispersal in L. saxatilis (i.e., snails can only disperse to the closest four localities)

# Convert to neighbor list:
nb_MN <- knn2nb(knn_MN)

# Calculate test:
moran.mc(PopLev_MN$EcolPC1, nb2listw(nb_MN), nsim = 999)
# OUTPUT: Moran's I = 0.66803, P = 0.001

##### 2.2.2.- Arousa ----

# Select the ria:
PopLev_A <- subset(PopLev, Ria == "Arousa") 

# Extract the coordinates:
Coordinates_A <- cbind(PopLev_A$Latitude, PopLev_A$Longitude) 

# Calculate spatial weights matrix:
knn_A <- knearneigh(Coordinates_A, k = 4) # Using K-nearest neighbors, k = 4

# Convert to neighbor list:
nb_A <- knn2nb(knn_A)

# Calculate test:
moran.mc(PopLev_A$EcolPC1, nb2listw(nb_A), nsim = 999)
# OUTPUT: Moran's I = 0.39799, P = 0.003

##### 2.2.3.- Pontevedra ----

# Select the ria:
PopLev_P <- subset(PopLev, Ria == "Pontevedra") 

# Extract the coordinates:
Coordinates_P <- cbind(PopLev_P$Latitude, PopLev_P$Longitude) 

# Calculate spatial weights matrix:
knn_P <- knearneigh(Coordinates_P, k = 4) # Using K-nearest neighbors, k = 4

# Convert to neighbor list:
nb_P <- knn2nb(knn_P)

# Calculate test:
moran.mc(PopLev_P$EcolPC1, nb2listw(nb_P), nsim = 999)
# OUTPUT: Moran's I = 0.22587, P = 0.013

##### 2.2.4.- Vigo ----

# Select the ria:
PopLev_V <- subset(PopLev, Ria == "Vigo") 

# Extract the coordinates:
Coordinates_V <- cbind(PopLev_V$Latitude, PopLev_V$Longitude) 

# Calculate spatial weights matrix:
knn_V <- knearneigh(Coordinates_V, k = 4) # Using K-nearest neighbors, k = 4

# Convert to neighbor list:
nb_V <- knn2nb(knn_V)

# Calculate test:
moran.mc(PopLev_V$EcolPC1, nb2listw(nb_V), nsim = 999)
# OUTPUT: Moran's I = 0.78253, P = 0.001

#### 2.3.- Plot with magnitude of variation ~ ecosystem clustering ----

# Data frame of values:
Ria <- c("Muros-Noia", "Arousa", "Pontevedra", "Vigo") # Set the rias and coasts
MoranI <- c(0.66803, 0.39799, 0.22587, 0.78253) # Set their corresponding Moran's I
Fulva_cov <- c(0.2362244, 0.09004023, 0.2278874, 0.4232778) # Extracted from the McFadden's pseudo-r^2 of the beta-binomial models
Lineata_cov <- c(0.2996492, 0.1598871, 0.1042499, 0.2683315) # Same as previous line
StrCov <- data.frame(Ria, MoranI, Fulva_cov, Lineata_cov) # Create a data frame
StrCov 

# Turn the data frame into a long format:
StrCov_long <- StrCov %>% 
  pivot_longer(
    cols = c(Fulva_cov, Lineata_cov),
    names_to = "ColorMorph",
    values_to = "Covariation") %>% 
  mutate(
    ColorMorph = if_else(ColorMorph == "Fulva_cov", "Fulva", "Lineata")
  )

# Order the Rias so that they are assigned the same color as in other figures:
StrCov_long <- StrCov_long %>% 
  mutate(
    Ria = factor(Ria, levels = c("Muros-Noia", "Arousa","Pontevedra", "Vigo")))

# Plot with aesthetics:
Fig_SM_StrCov <- ggplot(StrCov_long, aes(x = MoranI, y = Covariation)) +
  geom_point(aes(fill = Ria, shape = ColorMorph), 
             size = 5, alpha = 1,
             #show.legend = FALSE,
  ) +
  geom_smooth(aes(linetype = ColorMorph), method = "lm", se = FALSE, linewidth = 1, color = "black",
              #show.legend = FALSE,
  ) +
  labs(title = "",
       x = "Clustering of ecosystems \n (Moran's I)",
       y = "Strength of trait covariation") +
  theme_test() +
  scale_shape_manual(values = c(21, 23)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_cd(18) +
  scale_y_continuous(limits = c(0, 0.5),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5)) +  # Remove the annoying leading 0s
  scale_x_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1" = 1))  # Remove the annoying leading 0s

Fig_SM_StrCov # View the plot

# Save the plot:
ggsave("CD_MS_SM_Fig_S.svg",
       plot = Fig_SM_StrCov,
       width = 6,
       height = 5)
