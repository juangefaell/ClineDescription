############ ClineDescription - Ecosystem Equivalence ###############

# BASIC INFO ----

# Researcher: [AuthorName]
# Section aim: Are the ecosystems of the Rías Baixas ecologically equivalent?
# Last update: 2026-05-19

# SETUP SECTION ----

## 1.- Packages ----
library(tidyverse) # Basic plots and statistical analyses
library(FactoMineR) # Dimension Reduction Analyses (PCA and MCA)
library(factoextra) # Dimension Reduction Analyses (PCA and MCA)
library(vegan) # PERMANOVA, Bray-Curtis
library(car) # ANOVA and to check assumptions (Non-Constant Variance test)
library(mgcv) # GAM
library(DHARMa) # Model diagnostics
library(cowplot) # Composite figures
library(grid) # Composite figures 

## 2.- Dataset loading ----
setwd("") # Set the directory of the dataset file
read.csv("Data_PopulationLevel.csv")

## 3.- Data handling and overview ----

# PopulationLevel dataset: 
PopLev <- read.csv("CD_DA_RD_PopulationLevel.csv") %>%
  mutate(
    Ria = factor(Ria, levels = c("MurosNoia", "Arousa","Pontevedra", "Vigo")), # Set an order from north to south
    Coast = factor(Coast, levels = c("South", "North"))) # Convert variables into adequate formats
str(PopLev) # Check the structure of the data
view(PopLev) # View the whole dataset in a separate tab

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

# Theme for figures, v. 2:
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

### 1.- Principal Component Analysis (PCA) of ecological variables ----

# Subsection aim: Can we extract meaningful ecological information from the community composition?

# ///

# Format the Muros-Noia ria well:
PopLev <- PopLev %>% mutate(Ria = fct_recode(Ria, "Muros-Noia" = "MurosNoia"))

# Formatting variables for PCA, and giving them separate names:
PopLev <- PopLev %>%
  mutate(
    Ria = factor(Ria, levels = c("Muros-Noia", "Arousa","Pontevedra", "Vigo")),
    'A. nudosum' = as.numeric(Ascophyllum_Coverage),
    'Chthamalus sp.' = as.numeric(Barnacles_Coverage),
    'C. officinalis' = as.numeric(Corallina_Presence),
    'L. obtusata' = as.numeric(FabalisObtusata_Presence),
    'F. vesiculosus' = as.numeric(Fucus_Coverage),
    'G. umbilicalis' = as.numeric(Gibbula_Presence),
    'L. pygmaea' = as.numeric(Lichina_Coverage),
    'L. littorea' = as.numeric(Littorea_Presence),
    'M. neritoides' = as.numeric(Melarhaphe_Presence),
    'M. galloprovincialis' = as.numeric(Mytilus_Coverage),
    'N. lapillus' = as.numeric(Nucella_Presence),
    'P. canaliculata' = as.numeric(Pelvetia_Coverage),
    'P. lineatus' = as.numeric(Phorcus_Presence)
  ) 

# Selecting the variables:
EcolPCA_Vars <- PopLev[, c("A. nudosum", "Chthamalus sp.", 
                            "C. officinalis", "L. obtusata", 
                            "F. vesiculosus", "G. umbilicalis", 
                            "L. pygmaea", "L. littorea",
                            "M. neritoides", "M. galloprovincialis",
                            "N. lapillus", "P. canaliculata","P. lineatus")] 

# Check the format of the variables (IMPORTANT!): 
str(EcolPCA_Vars) 

# Perform PCA:
EcolPCA_Result <- prcomp(EcolPCA_Vars, center = TRUE, scale = TRUE) # Perform PCA, with variables scaled (although here is not strictly necessary due to ordinal data)
summary(EcolPCA_Result) # Display the % of variance explain by each PC
print(EcolPCA_Result) # Display the SD and component loadings of the individual variables
# OUTCOME: PC1 explains 49.19% of variance; PC2: 13.29%; PC3: 8.61%

# Visualization of PCA results (guide: https://www.sthda.com/english/wiki/fviz-pca-quick-principal-component-analysis-data-visualization-r-software-and-data-mining):
fviz_eig(EcolPCA_Result, addlabels = TRUE) # Scree plot (variance explained by each PC)
fviz_pca_biplot(EcolPCA_Result, label = "var", addEllipses = TRUE, ellipse.level = 0.95, repel = TRUE) # Biplot of variables to visualize the directions of variables (individuals can be added)

# Visualization of localities (i.e., individuals):
PCA_Ind <- fviz_pca_ind(EcolPCA_Result, 
                        label = "none", # Add or remove the number of locality ("ind" vs. "none")
                        pointshape = 21, pointsize = 2.5,
                        addEllipses = TRUE, ellipse.level = 0.95, 
                        repel = TRUE,
                        fill = PopLev$Ria,  # Compare distributions of localities by Ria
                        col.ind = "black",
                        axes = c(1, 2)) # Select PCs to display

# Plot of PCA values per ria (Figure 2B):
Fig_2B <- PCA_Ind + 
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(title = "",
       x = "PC1 (49.2%)",
       y = "PC2 (13.3%)") +
  theme_cd(18) +
  scale_y_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_2B # View the plot

# Extract PC scores (PC1 and PC2):
PC_Scores <- as.data.frame(EcolPCA_Result$x)

# Set Ria as grouping factor:
PC_Scores$Ria <- PopLev$Ria

# PERMANOVA (i.e., assessing whether centroids differ in the multivariate space):
PERMANOVA_PCs <- adonis2(PC_Scores[, c("PC1", "PC2")] ~ Ria, data = PC_Scores, method = "euclidean", permutations = 999)

PERMANOVA_PCs # View the result

# Checking dispersion:
anova(betadisper(dist(PC_Scores[, c("PC1", "PC2")]), PC_Scores$Ria, type = "centroid"))
# OUTPUT: No dispersion, so the PERMANOVA results are robust

# Extract the component loadings:
EcolPCA_Loadings_PC1 <- EcolPCA_Result$rotation[, 1] # Extract and select the component loadings ("rotation") for PC1 (corresponding to column 1)
EcolPCA_Loadings_df <- data.frame(Variable = names(EcolPCA_Loadings_PC1), # Create a dataframe of the loadings with appropriate names to create a figure in ggplot2
                                  Loading = EcolPCA_Loadings_PC1) 

EcolPCA_Loadings_df # Check that it worked

# Plot of PCA loadings (Figure 2A): 
Fig_2A <- ggplot(EcolPCA_Loadings_df, aes(x = reorder(Variable, Loading), y = Loading)) + # Create and customize the plot with variables ordered according to contribution
  geom_bar(stat = "identity", aes(fill = Loading), alpha = 1, color = "black") +
  scale_fill_gradient(low = "purple", high = "yellow") +
  coord_flip(ylim = c(-0.5, 0.5)) + # Set the graph in horizontal
  labs(title = "", x = "", y = "Component loading") +
  theme_cd(18) +
  theme(panel.grid.major.y = element_line(color = "grey80", linewidth = 0.8, linetype = "dashed"),
        panel.grid.major.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color = "black", margin = margin(t = 6)),
        axis.text.y = element_text(face = "italic", margin = margin(r = 6)),
        legend.position = "none") +
  scale_y_continuous(breaks = c("-0.5" = -0.5, "-0.25" = -0.25, "0" = 0, "0.25" = 0.25, "0.5" = 0.5)) # Remove the annoying "0.0" axis text 

Fig_2A # View the plot

# Save the "EcolPC1" and "EcolPC2" scores into the dataset file for further analyses:
  # EcolPC1 <- EcolPCA_Result$x[, 1]  # Obtain the PC1 scores (column number 1)
  # EcolPC2 <- EcolPCA_Result$x[, 2]  # Obtain the PC2 scores (column number 2)
  # PopLev <- read.csv("CD_DA_RD_PopulationLevel.csv") # Re-establish original format before saving (IMPORTANT!)
  # PopLev <- cbind(PopLev, EcolPC1, EcolPC2) # Add EcolPC1 and EcolPC2 scores to the dataset
  # str(PopLev) # Check that the format of the variables is the one given by default in R, not the ones used in PCA or MCA; check also that EcolPC1 and EcolPC2 are actually there
  # write.csv(PopLev, "CD_DA_RD_PopulationLevel.csv", row.names = FALSE) # Save a new version of the dataset with EcolPC1 and EcolPC2

### 2.- Ecological PC1 and wave exposure ----

# Subsection aim: Verify that Ecological PC1 is a biological indicator of wave exposure

# ///

# Run a linear model:
EcolPC1_LM_WaveAction <- lm(EcolPC1 ~ Wave_Mean, data = PopLev, na.action = na.omit) # Run a model to check assumptions

summary(EcolPC1_LM_WaveAction) # Summary of the model

# Check the assumptions:
SimulatedResiduals_EcolPC1_LM_WaveAction <- simulateResiduals(EcolPC1_LM_WaveAction)
plot(SimulatedResiduals_EcolPC1_LM_WaveAction)
# OUTPUT: The fit is not perfect...

# Generalized Additive Model (GAM):
EcolPC1_GAM_WaveAction <- gam(EcolPC1 ~ s(Wave_Mean), data = PopLev, method = "REML")

summary(EcolPC1_GAM_WaveAction) # Summary of the model

# Check diagnostics of GAM
gam.check(EcolPC1_GAM_WaveAction)
# OUTPUT: k-index is around 1, so good fit

# Compare the linear and GAM models:
AIC(EcolPC1_LM_WaveAction, EcolPC1_GAM_WaveAction)
# OUTPUT: The GAM model performs better

# Plot of Wave action ~ Ecological PC1 (Figure 2C):
Fig_2D <- ggplot(PopLev, aes(x = Wave_Mean, y = EcolPC1)) + 
  geom_point(aes(fill = Ria), size = 3, alpha = 1, stroke = 1, show.legend = FALSE, color = "black", shape = 21) +
  geom_smooth(method = "gam", se = TRUE, color = "black", aes(group = 1), linewidth = 1) + # Note the GAM fitting here
  labs(title = "",
       x = "Wave action",
       y = "Ecological PC1") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  theme_cd(18) +
  scale_x_reverse(labels = function(x) ifelse(x == 0, "0", x))  # Reverse the x-scale so that it matches the geographical distribution of localities
#coord_fixed(ratio = 0.05) +

Fig_2D # View the plot

### 3.- Ecological PC1 and position of localities ----

# Subsection aim: Do the four Rias vary equally across space?

# ///

# Plot of Ecological PC1 ~ position of localities (Figure 2D):
Fig_2C <- ggplot(data = subset(PopLev, !is.na(DistanceRiver_Normalized)), # Important to omit NAs as the locs from Illa de Arousa don't have distance from river (they are islands)
                 aes(x = DistanceRiver_Normalized, y = EcolPC1)) + 
  geom_point(aes(fill = Coast), size = 3, alpha = 1, stroke = 1, show.legend = FALSE, color = "black", shape = 21) +
  geom_smooth(method = "lm", se = TRUE, color = "black", aes(group = 1), linewidth = 1) +
  labs(title = "",
       x = "Distance from river mouth \n (normalized)",
       y = "Ecological PC1") +
  facet_wrap(~ Ria, ncol = 1) +
  scale_x_reverse(limits = c(1, 0), 
                  breaks = c(0, 0.25, 0.5, 0.75, 1),
                  labels = function(x) ifelse(x == 0, "0", x)) + # Reverse the x-scale so that it matches the geographical distribution of localities and remove the annoying leading 0s
  theme_cd2(18) +
  theme(panel.border = element_rect(linewidth = 1.5))

Fig_2C # View the plot

# ANCOVA: 
Ancova_EcolRias <- lm(EcolPC1 ~ DistanceRiver_Normalized * Ria, data = PopLev)
summary(Ancova_EcolRias)

# Linear model (Muros-Noia):
# Careful with this: PopLev <- PopLev %>% mutate(Ria = fct_recode(Ria, "Muros-Noia" = "MurosNoia"))

EcolPC1_lm_MN <- PopLev %>%
  filter(Ria == "Muros-Noia") %>%
  lm(EcolPC1 ~ DistanceRiver_Normalized, data = .)
summary(EcolPC1_lm_MN)

EcolPC1_lm_MN_res <- resid(EcolPC1_lm_MN) # Extract the residuals to check assumptions
plot(density(EcolPC1_lm_MN_res), main = "Density Plot of EcolPC1 residuals (Muros-Noia)")
shapiro.test(EcolPC1_lm_MN_res)
ncvTest(EcolPC1_lm_MN)
?shapiro.test
# OUTPUT: The residuals are normal and homoscedastic

# Linear model (Arousa):
EcolPC1_lm_A <- PopLev %>% 
  filter(Ria == "Arousa") %>%
  lm(EcolPC1 ~ DistanceRiver_Normalized, data = .)
summary(EcolPC1_lm_A)

EcolPC1_lm_A_res <- resid(EcolPC1_lm_A) # Extract the residuals to check assumptions
plot(density(EcolPC1_lm_A_res), main = "Density Plot of EcolPC1 residuals (Arousa)")
shapiro.test(EcolPC1_lm_A_res)
ncvTest(EcolPC1_lm_A)
# OUTPUT: The residuals are normal and homoscedastic

## Linear model (Pontevedra):
EcolPC1_lm_P <- PopLev %>%
  filter(Ria == "Pontevedra") %>%
  lm(EcolPC1 ~ DistanceRiver_Normalized, data = .)
summary(EcolPC1_lm_P)

EcolPC1_lm_P_res <- resid(EcolPC1_lm_P) # Extract the residuals to check assumptions
plot(density(EcolPC1_lm_P_res), main = "Density Plot of EcolPC1 residuals (Pontevedra)")
shapiro.test(EcolPC1_lm_P_res)
ncvTest(EcolPC1_lm_P)
# OUTPUT: The residuals are normal but not homoscedastic

# Linear model (Vigo):
EcolPC1_lm_V <- PopLev %>% 
  filter(Ria == "Vigo") %>%
  lm(EcolPC1 ~ DistanceRiver_Normalized, data = .)
summary(EcolPC1_lm_V) 

EcolPC1_lm_V_res <- resid(EcolPC1_lm_V) # Extract the residuals to check assumptions
plot(density(EcolPC1_lm_V_res), main = "Density Plot of EcolPC1 residuals (Vigo)")
shapiro.test(EcolPC1_lm_V_res)
ncvTest(EcolPC1_lm_V)
# OUTPUT: The residuals follow a normal distribution and are homoscedastic

### 4.- Ecological differences between ecosystem types ----

# Subsection aim: How different are wave-sheltered and wave-exposed ecosystems in terms of species composition?

# ///

#### 4.1.- Bray-Curtis Index calculations between "WaveSheltered" and "WaveExposed" per ria ----

# NOTE: For the creation of "WaveSheltered" and "WaveExposed" categories, see the Hierarchical clustering in "Supplementary analyses"

# Muros-Noia (Careful with the format: Go to "2.- Dataset loading" and reload the dataset to avoid compatibility problems due to "Muros-Noia" with hyphen):
WaveSheltered_MurosNoia <- subset(PopLev, Cluster == "WaveSheltered" & Ria == "MurosNoia", # Subset and select the variables for WaveExposed
                                select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                           FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                           Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                           Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))

WaveExposed_MurosNoia <- subset(PopLev, Cluster == "WaveExposed" & Ria == "MurosNoia", # Subset and select the variables for WaveExposed
                                select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                           FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                           Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                           Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))

WaveSheltered_MurosNoia_sum <- colSums(WaveSheltered_MurosNoia, na.rm = TRUE) # Calculate the sums of presence values for each species
WaveExposed_MurosNoia_sum <- colSums(WaveExposed_MurosNoia, na.rm = TRUE)

a_MurosNoia <- sum((WaveSheltered_MurosNoia_sum > 0) & (WaveExposed_MurosNoia_sum > 0)) # Calculate the number of species present in both ecosystems (a)
b_MurosNoia <- sum((WaveSheltered_MurosNoia_sum > 0) & (WaveExposed_MurosNoia_sum == 0)) # Calculate the number of species present in WaveSheltered but not in WaveExposed (b)
c_MurosNoia <- sum((WaveExposed_MurosNoia_sum > 0) & (WaveSheltered_MurosNoia_sum == 0)) # Calculate the number of species present in WaveExposed but not in WaveSheltered (c)

BrayCurtis_MurosNoia <- (b_MurosNoia + c_MurosNoia) / (2 * a_MurosNoia + b_MurosNoia + c_MurosNoia) # Calculate the Bray-Curtis similarity index

BrayCurtis_MurosNoia # View the result

# Arousa: 
WaveSheltered_Arousa <- subset(PopLev, Cluster == "WaveSheltered" & Ria == "Arousa", # Subset and select the variables for WaveSheltered
                               select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                          FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                          Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                          Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))

WaveExposed_Arousa <- subset(PopLev, Cluster == "WaveExposed" & Ria == "Arousa", # Subset and select the variables for WaveExposed
                             select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                        FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                        Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                        Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))

WaveSheltered_Arousa_sum <- colSums(WaveSheltered_Arousa, na.rm = TRUE) # Calculate the sums of presence values for each species
WaveExposed_Arousa_sum <- colSums(WaveExposed_Arousa, na.rm = TRUE)

a_Arousa <- sum((WaveSheltered_Arousa_sum > 0) & (WaveExposed_Arousa_sum > 0)) # Calculate the number of species present in both ecosystems (a)
b_Arousa <- sum((WaveSheltered_Arousa_sum > 0) & (WaveExposed_Arousa_sum == 0)) # Calculate the number of species present in WaveSheltered but not in WaveExposed (b)
c_Arousa <- sum((WaveExposed_Arousa_sum > 0) & (WaveSheltered_Arousa_sum == 0)) # Calculate the number of species present in WaveExposed but not in WaveSheltered (c)

BrayCurtis_Arousa <- (b_Arousa + c_Arousa) / (2 * a_Arousa + b_Arousa + c_Arousa) # Calculate the Bray-Curtis similarity index

BrayCurtis_Arousa # View the result

# Pontevedra:
WaveSheltered_Pontevedra <- subset(PopLev, Cluster == "WaveSheltered" & Ria == "Pontevedra", # Subset and select the variables for WaveSheltered
                                   select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                              FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                              Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                              Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))
WaveExposed_Pontevedra <- subset(PopLev, Cluster == "WaveExposed" & Ria == "Pontevedra", # Subset and select the variables for WaveExposed
                                 select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                            FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                            Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                            Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))

WaveSheltered_Pontevedra_sum <- colSums(WaveSheltered_Pontevedra, na.rm = TRUE) # Calculate the sums of presence values for each species
WaveExposed_Pontevedra_sum <- colSums(WaveExposed_Pontevedra, na.rm = TRUE)

a_Pontevedra <- sum((WaveSheltered_Pontevedra_sum > 0) & (WaveExposed_Pontevedra_sum > 0)) # Calculate the number of species present in both ecosystems (a)
b_Pontevedra <- sum((WaveSheltered_Pontevedra_sum > 0) & (WaveExposed_Pontevedra_sum == 0)) # Calculate the number of species present in WaveSheltered but not in WaveExposed (b)
c_Pontevedra <- sum((WaveExposed_Pontevedra_sum > 0) & (WaveSheltered_Pontevedra_sum == 0)) # Calculate the number of species present in WaveExposed but not in WaveSheltered (c)

BrayCurtis_Pontevedra <- (b_Pontevedra + c_Pontevedra) / (2 * a_Pontevedra + b_Pontevedra + c_Pontevedra) # Calculate the Bray-Curtis similarity index

BrayCurtis_Pontevedra # View the result

# Vigo:
WaveSheltered_Vigo <- subset(PopLev, Cluster == "WaveSheltered" & Ria == "Vigo", # Subset and select the variables for WaveSheltered
                               select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                          FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                          Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                          Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))

WaveExposed_Vigo <- subset(PopLev, Cluster == "WaveExposed" & Ria == "Vigo", # Subset and select the variables for WaveExposed
                             select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                        FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                        Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                        Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))

WaveSheltered_Vigo_sum <- colSums(WaveSheltered_Vigo, na.rm = TRUE) # Calculate the sums of presence values for each species
WaveExposed_Vigo_sum <- colSums(WaveExposed_Vigo, na.rm = TRUE) 

a_Vigo <- sum((WaveSheltered_Vigo_sum > 0) & (WaveExposed_Vigo_sum > 0)) # Calculate the number of species present in both ecosystems (a)
b_Vigo <- sum((WaveSheltered_Vigo_sum > 0) & (WaveExposed_Vigo_sum == 0)) # Calculate the number of species present in WaveSheltered but not in WaveExposed (b)
c_Vigo <- sum((WaveExposed_Vigo_sum > 0) & (WaveSheltered_Vigo_sum == 0)) # Calculate the number of species present in WaveExposed but not in WaveSheltered (c)

BrayCurtis_Vigo <- (b_Vigo + c_Vigo) / (2 * a_Vigo + b_Vigo + c_Vigo) # Calculate the Bray-Curtis similarity index

BrayCurtis_Vigo # View the result
# NOTE: The equivalence of Arousa and Vigo values is not an artifact; checked on 2025-11-21)

#### 4.2.- What are the species that differ most between ecosystems? (Rias pooled) ----

# Select cases from WaveSheltered localities: 
WaveSheltered <- subset(PopLev, Cluster == "WaveSheltered", 
                                  select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                             FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                             Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                             Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))
# Select cases from WaveExposed localities:
WaveExposed <- subset(PopLev, Cluster == "WaveExposed",
                                select = c(Ascophyllum_Presence, Barnacles_Presence, Corallina_Presence,
                                           FabalisObtusata_Presence, Fucus_Presence, Gibbula_Presence,
                                           Lichina_Presence, Littorea_Presence, Melarhaphe_Presence,
                                           Nucella_Presence, Pelvetia_Presence, Phorcus_Presence))

# Calculate the sums of presence values for each species:
WaveSheltered_sum <- colSums(WaveSheltered, na.rm = TRUE) # WaveSheltered
WaveExposed_sum <- colSums(WaveExposed, na.rm = TRUE) # WaveExposed

WaveSheltered_sum # Less prevalent species: Corallina (0), Melarhaphe (0), Lichina (1), Nucella (2)
WaveExposed_sum # Less prevalent species: Ascophyllum (0), Littorina obtusata (1), Fucus (2), Littorina littorea (3)

#### 4.3.- Plot of Bray-Curtis differences between wave-sheltered and wave-exposed per ria ----

# Creating the variables and data frame:
Rias <- c("Muros-Noia", "Arousa", "Pontevedra", "Vigo") # Set the rias
Sorensen <- c(0.6666667, 0.5882353, 0.5714286, 0.5882353) # Set their corresponding Sørensen indexes
BrayCurtis <- c(0.3333333, 0.4117647, 0.4285714, 0.4117647) # Set their corresponding Bray-Curtis indexes
Dissimilarity_df <- data.frame(Rias, Sorensen, BrayCurtis) # Create the data frame
Dissimilarity_df$Rias <- factor(Dissimilarity_df$Rias, levels = Rias) # Convert rias to a factor to set the desired order

# Plot of Bray-Curtis across rias (Figure 2E):
Fig_2E <- ggplot(Dissimilarity_df, aes(x = Rias, y = BrayCurtis, group = 1)) +
  geom_point(size = 3, color = "black") +
  geom_line(linewidth = 1, color = "black") +
  labs(title = "",
       x = "Ria",
       y = "Bray-Curtis \n (shelt. vs. exp.)") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = c("0" = 0, "0.25" = 0.25, "0.5" = 0.5, "0.75" = 0.75, "1.0" = 1)) + # Remove the annoying leading 0s
  theme_cd(18) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank())

Fig_2E # View the plot

### 5.- Composite figure (Figure 2) ----

# Set margins:
Margin <- theme(plot.margin = margin(4, 4, 4, 4, "mm"))
Fig_2A <- Fig_2A + Margin
Fig_2B <- Fig_2B + Margin
Fig_2C <- Fig_2C + Margin
Fig_2D <- Fig_2D + Margin
Fig_2E <- Fig_2E + Margin

# Remove legends (they distort the alignment):
Fig_2B <- Fig_2B + theme(legend.position = "none")
Fig_2C <- Fig_2C + theme(legend.position = "none")
Fig_2D <- Fig_2D + theme(legend.position = "none")

# Add tags: 
Fig_2A_Tag <- tags_cd(Fig_2A, "A", size = 22, x = 0.02, y = 0.98)
Fig_2B_Tag <- tags_cd(Fig_2B, "B", size = 22, x = 0.02, y = 0.98)
Fig_2C_Tag <- tags_cd(Fig_2C, "C", size = 22, x = 0.005, y = 0.98)
Fig_2D_Tag <- tags_cd(Fig_2D, "D", size = 22, x = 0.02, y = 0.98)
Fig_2E_Tag <- tags_cd(Fig_2E, "E", size = 22, x = 0.02, y = 1) # Remember the change here so as to not overlap

# Create the left column:
Fig_2_L  <- plot_grid(Fig_2A_Tag, Fig_2D_Tag,
                      ncol = 1,
                      rel_heights = c(1.0, 1.0),
                      align = "v", axis = "l")

# Create the center column:
Fig_2_C  <- plot_grid(Fig_2B_Tag, Fig_2E_Tag,
                      ncol = 1,
                      rel_heights = c(1.0, 1.0),
                      align = "v", axis = "l")

# Create the right column:
Fig_2_R <- plot_grid(Fig_2C_Tag,
                      ncol = 1,
                      rel_heights = c(1.0),
                      align = "v", axis = "l")

# Merge the two columns:
Figure2 <- plot_grid(Fig_2_L, Fig_2_C, Fig_2_R,
                     ncol = 3,
                     rel_widths = c(1.4, 1.1, 1.0),
                     align = "h", axis = "t")

# Save in .svg to export to InkScape:
ggsave("Figure2.svg",
       Figure2,
       width = 460, height = 235, units = "mm")

## Supplementary analyses ----

### 1.- Hierarchical clustering to create groups of localities ----

# Select the variables (presence/absence variables) and without M. galloprovincialis, which is not captured by EcolPC1: 
PopLev_HC <- PopLev[, c("Ria", "LocName", "Ascophyllum_Presence", "Barnacles_Presence", "Corallina_Presence", 
                        "FabalisObtusata_Presence", "Fucus_Presence", "Gibbula_Presence",
                        "Lichina_Presence", "Littorea_Presence", "Melarhaphe_Presence",
                        "Nucella_Presence", "Pelvetia_Presence", "Phorcus_Presence")] 

# Remove Ria and LocName: 
PopLev_HC_Species <- PopLev_HC[,3:14]

# Compute the distance matrix on which HC will be calculated:
HC_DistMat <- vegdist(PopLev_HC_Species, distance = "bray", binary = TRUE) 

# Perform the hierarchical clustering:
HC_Ecol <- hclust(HC_DistMat, method = "average") # Run the actual test
HC_Ecol$labels <- PopLev_HC$LocName # Assign labels with LocName to track localities

# Plot the hierarchical clustering results:
plot(HC_Ecol,
     main = "Hierarchical clustering of localities\n (Rias pooled)",
     xlab = "",
     sub = "", 
     ylab = "Height",
     cex = 1.2,
     lwd = 1,
     col = "black")

# Visualize the corresponding elbow plot:
Fig_SM_Elbow <- fviz_nbclust(PopLev_HC_Species, 
                             FUN = hcut,
                             method = "wss",
                             hc_method = "average",
                             diss = HC_DistMat) +
  labs(x = "Number of clusters",
       y = "Total within-cluster sum of squares") +
  theme_cd(18)

Fig_SM_Elbow # View the plot
# OUTPUT: It seems that the best number of clusters is 3

# Save the plot:
 ggsave("CD_MS_SM_Fig_S6.svg",
       plot = Fig_SM_Elbow,
       width = 6,
       height = 5)

# Assign localities to clusters:
Clusters <- cutree(HC_Ecol, k = 3)

# Add the corresponding variable to the dataset:
PopLev$Cluster <- Clusters # Create the variable and add it to the dataset

# Turn the variable into a character variable (to rightly afterwards plot it):
PopLev$Cluster <- as.character(PopLev$Cluster) 

# Plot to visualize how Cluster relates to EcolPC1 (with aesthetics):
Fig_SM_Clust <- PopLev %>% 
  ggplot(aes(x = Cluster, y = EcolPC1, fill = Cluster)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.1), alpha = 0.5, color = "black") +
  scale_fill_manual(values = c("#63069D", "#C73890", "#EBD809")) +
  labs(x = "Cluster",
       y = "Ecological PC1") +
  theme_cd(18) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = c("-5" = -5, "-2.5" = -2.5, "0" = 0, "2.5" = 2.5))

Fig_SM_Clust # View the plot

# Test differences:
Anova_Clust <- aov(EcolPC1 ~ Cluster, data = PopLev) # Run an ANOVA to formally assess differences

summary(Anova_Clust) # See the results

hist(Anova_Clust$residuals) # The model seems to perform well
# OUTPUT: The clusters map EcolPC1 fairly well 

# Tukey test for pairwise differences:
Tukey_Clust <- TukeyHSD(Anova_Clust)

# Save the plot:
ggsave("CD_MS_SM_Fig_S7.svg",
       plot = Fig_SM_Clust,
       width = 6,
       height = 5)

# Save the "Cluster" variable in the dataset, with new names according to their mapping in the EcolPC1 variable:
# PopLev$Cluster <- recode(PopLev$Cluster,
#                         `1` = "WaveSheltered",
#                         `2` = "WaveIntermediate",
#                         `3` = "WaveExposed") 

# write.csv(PopLev, "CD_DA_RD_PopulationLevel.csv", row.names = FALSE) # Save the updated dataset

### 2.- Multiple Correspondence Analysis (MCA) ----

# Formatting variables for MCA, and giving them separate names:
PopLev <- read.csv("CD_DA_RD_PopulationLevel.csv") # Ensure the original format before performing MCA (IMPORTANT!)
PopLev <- PopLev %>% # Give a factor format and new names to variables of interest
  mutate(
    Ria = factor(Ria, levels = c("MurosNoia", "Arousa","Pontevedra", "Vigo")),
    'Ascophyllum nudosum' = as.factor(Ascophyllum_Presence), # Note that the "_Abundance" variables have been substituted for "Presence/Absence"
    'Barnacles' = as.factor(Barnacles_Presence), # Turned into "Presence/Absence"
    'Corallina officinalis' = as.factor(Corallina_Presence),
    'Littorina obtusata' = as.factor(FabalisObtusata_Presence),
    'Fucus vesiculosus' = as.factor(Fucus_Presence), # Turned into "Presence/Absence"
    'Gibbula umbilicalis' = as.factor(Gibbula_Presence),
    'Lichina pygmaea' = as.factor(Lichina_Presence),
    'Littorina littorea' = as.factor(Littorea_Presence),
    'Melarhaphe neritoides' = as.factor(Melarhaphe_Presence),
    'Mytilus galloprovincialis' = as.factor(Mytilus_Presence), # Turned into "Presence/Absence"
    'Nucella lapillus' = as.factor(Nucella_Presence),
    'Pelvetia canaliculata' = as.factor(Pelvetia_Presence), # Turned into "Presence/Absence"
    'Phorcus lineatus' = as.factor(Phorcus_Presence)
  ) 

# Selecting the variables:
EcolMCA_Vars <- PopLev[, c("Ascophyllum nudosum", "Barnacles", 
                           "Corallina officinalis", "Littorina obtusata", 
                           "Fucus vesiculosus", "Gibbula umbilicalis", "Lichina pygmaea", 
                           "Littorina littorea", "Melarhaphe neritoides", 
                           "Mytilus galloprovincialis", "Nucella lapillus", "Pelvetia canaliculata", 
                           "Phorcus lineatus")] 

# Check the format of the variables (IMPORTANT!): 
str(EcolMCA_Vars) 

# Perform MCA:
EcolMCA_Result <- MCA(EcolMCA_Vars, graph = FALSE) # Perform MCA
summary(EcolMCA_Result) # Display the % variance explain by each dimension

# Scree plot of MCA: 
EcolMCA_ScreePlot <- fviz_eig(EcolMCA_Result, addlabels = TRUE, ylim = c(0, 50)) # Scree plot (variance explained by each dimension)

EcolMCA_ScreePlot <- EcolMCA_ScreePlot +
  theme_minimal() + 
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(color = "black", hjust = 1, size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 10, b =10, l = 65, r = 18))

# Categories' plot of MCA:
EcolMCA_CategoriesPlot <- fviz_mca_var(EcolMCA_Result, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Visualize the quality of representation (cos2) and contributions of variables to Dims 1 and 2
                                       repel = TRUE, labelsize = 4) 

EcolMCA_CategoriesPlot <- EcolMCA_CategoriesPlot +
  theme_minimal() + 
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(color = "black", hjust = 1, size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 10, b =10, l = 65, r = 18))

# Loadings plot:
EcolMCA_LoadingsPlot <- fviz_contrib(EcolMCA_Result, choice = "var", axes = 1, top = 10) 

EcolMCA_LoadingsPlot <- EcolMCA_Loadings_Plot +
  theme_classic() + 
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 15),
        axis.text.y = element_text(color = "black", size = 15),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 10, b =10, l = 65, r = 18))

# Add the tags:
EcolMCA_ScreePlot <- tags_cd(EcolMCA_ScreePlot, "A", size = 22, x = 0.001, y = 0.99)

EcolMCA_LoadingsPlot <- tags_cd(EcolMCA_LoadingsPlot, "B", size = 22, x = 0.001, y = 0.99)

# Composite plot: 
Fig_SM_MCA <- plot_grid(EcolMCA_ScreePlot, EcolMCA_LoadingsPlot,
          ncol = 1, axis = "v", align = "tb")

# Save the plot:
ggsave("CD_MS_SM_Fig_S9.svg",
       plot = Fig_SM_MCA,
       width = 8,
       height = 9)

### 3.- Differences in rock types across rias ----

# Create a contingency table:
Rock_Type_Ria <- table(PopLev$Rock_Type, PopLev$Ria) 

print(Rock_Type_Ria) # View it

# Perform the Fisher's exact test:
fisher.test(Rock_Type_Ria) %>% 
  print()



