# Author: Shawna L. Rowe
# Last Update:  03-Feb-2025
# Associated manuscript: "Evidence of Rapid Evolution in Herbivory Defense Responses with Conserved Trade-offs in Populations of Medicago polymorpha"

# This code is functional as is. Please see the README for more information on version numbers. 

# ========  Packages  ===========
library(broom) 
library(car) 
library(cocor) 
library(corrplot) 
library(cowplot)
library(DHARMa) 
library(emmeans) 
library(Hmisc) 
library(lme4) 
library(lmerTest) 
library(patchwork)
library(plotrix) 
library(RcmdrMisc)
library(tidyverse) 

# ==== Datasets needed for protein analysis ====
ProteinDS <- read.csv("data/processed/ProteinData22Feb2017_edit_10jan24.csv") # first batch of protein data
ProteinNew <- read.csv("data/processed/ProteinAbsorbanceData.csv", header = TRUE, stringsAsFactors = FALSE) # second bath of protein data
PODDS <- read.csv("data/processed/WSU_PODfiles11Sept2017.csv") # Peroxidase
PPODS <- read.csv("data/processed/WSU_PPOfiles11Sept2017.csv") # Polyphenol oxidase

# Datasets of herbivore preference data ====
# 2022 induced dataset for preferences
ICpref <- read.csv("data/processed/InducedRangeDataset3.csv", header = TRUE, stringsAsFactors = FALSE) 
# 2018 data comparing within range induced and const preferences
ConstVInduc <- read.csv("data/processed/ConsVInduc6Sep2018.csv", header = TRUE, stringsAsFactors = FALSE) 
# 2018 data comparing ranges of const preferences
ConstRange <- read.csv("data/processed/ConstRange6Sep2018.csv", header = TRUE, stringsAsFactors = FALSE) 

# ======== functions ===========
# plot mean function for graphing
plot.mean <- function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  ymin <- m - se
  ymax <- m + se
  return(c(y = m, ymin = ymin, ymax = ymax))
}

# ========cafeteria/ preference analysis ======== 
## ==== 2017 cafeteria assay analysis ====
# Herbivore preference for constitutive plants native v invasive
ConstRange <- ConstRange %>%
  select(1, 2, 4, 7:9, 12:14)

ConstRange <- ConstRange %>%
  rename(LfWt.Initial_Invasive = LfWt_Initial_Invasive, LfWt.Final_Invasive = LfWt_Final_Invasive, LfWt.Initial_Native = LfWt_Initial_Native, LfWt.Final_Native = LfWt_Final_Native)

ConstRange <- ConstRange %>%
  dplyr::mutate(Eaten_Invasive = LfWt.Initial_Native - LfWt.Final_Native, Eaten_Native = LfWt.Initial_Invasive - LfWt.Final_Invasive, Proportion_Invasive = Eaten_Invasive / (Eaten_Invasive + Eaten_Native), Proportion_Native = Eaten_Native / (Eaten_Invasive + Eaten_Native))

ConstRange <- ConstRange %>%
  mutate(WinLoss = ifelse(Proportion_Invasive - Proportion_Native > 0, 1, 0))

ConstRange2 <- ConstRange %>%
  pivot_longer(c(3:8, 10:13),
               names_to = c(".value", "Range"),
               names_sep = "_"
  )

# summarise data. Group_by so we know for each caterpillar and range. Summarise all columns with numeric data and the functions performed are mean and std.error
ConstRange3 <- ConstRange2 %>%
  group_by(Herbivore, Range) %>%
  summarise(across(
    .cols = where(is.numeric),
    .fns = list(Mean = ~ mean(., na.rm = TRUE), SE = ~ std.error(., na.rm = TRUE)),
    .names = "{col}_{fn}"
  ))

write.table(ConstRange2, file = "~/Desktop/ConstRange_long-2024-04-28.csv", sep = ",", row.names = F)

ConstVInduc <- ConstVInduc %>%
  rename(LfWt.Initial_Induced = LfWt_Initial_Induced, LfWt.Final_Induced = LfWt_Final_Induced, LfWt.Initial_Const = LfWt_Initial_Const, LfWt.Final_Const = LfWt_Final_Const) %>%
  mutate(Eaten_Induced = ifelse(LfWt.Initial_Induced - LfWt.Final_Induced < 0, 0, LfWt.Initial_Induced - LfWt.Final_Induced), Eaten_Const = ifelse(LfWt.Initial_Const - LfWt.Final_Const < 0, 0, LfWt.Initial_Const - LfWt.Final_Const), Proportion_Induced = Eaten_Induced / (Eaten_Induced + Eaten_Const), Proportion_Const = Eaten_Const / (Eaten_Induced + Eaten_Const))

ConstVInduc2 <- ConstVInduc %>%
  pivot_longer(c(2:7, 10:13),
               names_to = c(".value", "Defense"),
               names_sep = "_"
  )

# summarise data. Group_by so we know for each caterpillar and range. Summarise all columns with numeric data and the functions performed are mean and std.error
ConstVInduc3 <- ConstVInduc2 %>%
  group_by(Herbivore, Range, Defense) %>%
  summarise(across(
    .cols = where(is.numeric),
    .fns = list(Mean = ~ mean(., na.rm = TRUE), SE = ~ std.error(., na.rm = TRUE)),
    .names = "{col}_{fn}"
  ))

#Figure 3: graph proportion eaten by range with caterpillar facet
(old_IC_proportion_range <- ggplot(ConstVInduc3, aes(x = Range, y = Proportion_Mean, color = Defense)) +
    geom_point() +
    geom_errorbar(aes(ymin = Proportion_Mean - Proportion_SE, ymax = Proportion_Mean + Proportion_SE, width = 0.25)) +
    guides(size = "none") +
    theme_bw() +
    facet_wrap(~Herbivore) +
    labs(
      title = "Feeding Preference for Constitutive vs Induced Plants",
      y = "Proportion Consumed",
      x = "",
      color = "Defense"
    ) +
    scale_color_manual(values = c("Const" = "#14C7BA", "Induced" = "#B323A5"),
                       labels = c("Const" = "Constitutive", "Induced" = "Induced")) + # Specify colorblind-friendly colors
    scale_x_discrete(labels = c("Invasive" = "Familiar", "Native" = "Unfamiliar")) +
    facet_wrap(~ Herbivore, labeller = as_labeller(c("Soybean" = "Soybean Looper", "Velvetbean" = "Velvetbean Caterpillar"))) +
    ylim(0.25, 0.75) + 
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "black") +  # Corrected ylim placement
    theme(
      strip.text.x = element_text(color = "black", size = 20, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.position = c(0.9,0.9),
      legend.background = element_rect(colour = "black", fill = "white", linetype = "solid"),
      #remove legend title
      legend.title = element_blank(),
      legend.text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
)
ggsave("Const2017_cafeteria_aug2024.jpg", old_IC_proportion_range, width = 3300, height = 2100, units = "px", dpi = 300)

# save in E and E format. 600 dpi PDF with an min 1800 px canvas size
ggsave("Figure3_const2017_cafeteria_aug2024.pdf", old_IC_proportion_range, width = 6600, height = 4200, units = "px", dpi = 600)

# ====== Induced2022 cafeteria ======
# create win-loss column in binary
# winloss column created by subtracting the
ICpref$win_loss <- ifelse(ICpref$WinnerTest > 0, 1, 0)

# select columns needed for analysis in ICpref data (all induced data)
ICpref2 <- ICpref %>%
  select(Plate_number, Invasive_accession, Native_accession, Pair_number, Caterpillar, Replicate, prop_wtmg_Invasive, Prop_wtmg_Native, win_loss)

ICpref2a <- ICpref2 %>%
  select(Plate_number, Invasive_accession, Native_accession, Pair_number, Caterpillar, Replicate, win_loss)

# names_pattern = "(.*)_accession" removes the "_accession" bit from any matching pattern
ICpref3a <- ICpref2a %>%
  pivot_longer(cols = c(Invasive_accession, Native_accession), names_to = "Range", values_to = "Accession", names_pattern = "(.*)_accession")

ICpref2b <- ICpref2 %>%
  select(Plate_number, Pair_number, Caterpillar, Replicate, prop_wtmg_Invasive, Prop_wtmg_Native, win_loss)

ICpref3b <- ICpref2b %>%
  pivot_longer(cols = c(prop_wtmg_Invasive, Prop_wtmg_Native), names_to = "Range", names_prefix = "[Pp](.*)_wtmg_", values_to = "ProportionEaten")

ICpref3 <- full_join(ICpref3a, ICpref3b)

# change all values of 1 (100% tissue eaten) to 0.9999999 for analysis
ICpref3 <- ICpref3 %>%
  mutate(ProportionEaten = ifelse(ProportionEaten > 0.9999, 0.9999999, ProportionEaten))

ICpref3 <- ICpref3 %>%
  filter(ProportionEaten > 0)

ICpref3$Replicate <- as.factor(ICpref3$Replicate)

# summarise data. Group_by so we know for each caterpillar and range. Summarise all columns with numeric data and the functions performed are mean and std.error
ICpref4 <- ICpref3 %>%
  group_by(Caterpillar, Range) %>%
  summarise(across(
    .cols = where(is.numeric),
    .fns = list(Mean = ~ mean(., na.rm = TRUE), SE = ~ std.error(., na.rm = TRUE)),
    .names = "{col}_{fn}"
  ))

# graph proportion eaten by range with caterpillar facet -- winning figure
(ICpref_proportion_range <- ggplot(ICpref4, aes(x = Range, y = ProportionEaten_Mean)) +
    labs(title = "Proportion of Leaf Tissue Eaten by Herbivore Range", x = expression(italic("M. polymorpha") ~ "Range"), y = "Proportion Eaten") +
    geom_point() +
    geom_errorbar(aes(ymin = ProportionEaten_Mean - ProportionEaten_SE, ymax = ProportionEaten_Mean + ProportionEaten_SE, width = 0.25)) +
    theme_bw() +
    facet_wrap(~Caterpillar, labeller = as_labeller(custom_labels)) +
    scale_x_discrete(labels = c("Invasive" = "Familiar", "Native" = "Unfamiliar")) +
    ylim(0.75, 0.95) +
    scale_color_grey() +
    theme(
      strip.text.x = element_text(color = "black", size = 20, hjust = 0.5, face = "bold"),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 22, hjust = 0.5),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16), panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
)

##---- Figure 2: Panel figure for Cafeteria data Oct24 ----
#graph proportion eaten by each caterpillar within each range. Facet by range and have the x-axis have the herbivore oct24
(ICpref_proportion_caterpillar <- ggplot(ICpref4, aes(x = Caterpillar, y = ProportionEaten_Mean)) +
   geom_point() +
   geom_errorbar(aes(ymin = ProportionEaten_Mean - ProportionEaten_SE, ymax = ProportionEaten_Mean + ProportionEaten_SE, width = 0.25)) +
   labs(title = "Proportion of Leaf Tissue Eaten by Herbivore", x = "", y = "Proportion Consumed") +
   theme_bw() +
   facet_wrap(~Range, labeller = as_labeller(c("Invasive" = "Familiar", "Native" = "Unfamiliar"))) +
   scale_x_discrete(labels = c("soybean" = "Soybean \nLooper", "velvetbean" = "Velvetbean \nCaterpillar")) +
   ylim(0.75, 0.95) +
   theme_bw() +
   scale_colour_grey() +
   theme(strip.text.x = element_text(color = "black", size = 16, face = "bold")) +
   theme(legend.position = "none") +
   theme(strip.text = element_text(size = 16, face = "bold")) +
   theme(axis.text = element_text(size = 14)) +
   theme(axis.title = element_text(size = 16)) +
   theme(plot.title = element_text(size = 18)) +
   theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#Graph the cafeteria data. Facet by caterpillar type and within the panels have familiar vs unfamiliar for only the constitutive tissue oct24
(cafeteria_plot_const <- ggplot(ConstRange2, aes(x = Range, y = Proportion)) +
    stat_summary(fun.data = plot.mean, geom = "errorbar", width = 0.25) +
    stat_summary(fun = mean, geom = "point", size = 2, color = "black", shape = 16, position = position_dodge(width = 0.35)) + # Adds the mean point
    guides(size = "none") +
    facet_wrap(~Herbivore) +
    labs(
      title = "Feeding Preference for Constitutive Plants",
      y = "Proportion Consumed",
      x = "",
      color = "Defense"
    ) +
    scale_x_discrete(labels = c("Invasive" = "Familiar", "Native" = "Unfamiliar")) +
    facet_wrap(~ Herbivore, labeller = as_labeller(c("Soybean" = "Soybean Looper", "Velvetbean" = "Velvetbean Caterpillar"))) +
    ylim(0.4, 0.6) +
    theme_bw() +
    scale_colour_grey() +
    theme(strip.text.x = element_text(color = "black", size = 16, face = "bold")) +
    theme(legend.position = "none") +
    theme(strip.text = element_text(size = 16, face = "bold")) +
    theme(axis.text = element_text(size = 14)) +
    theme(axis.title = element_text(size = 16)) +
    theme(plot.title = element_text(size = 18)) +
    theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#Graph the cafeteria data. Facet by caterpillar type and within the panels have familiar vs unfamiliar for only the induced tissue oct24
(cafeteria_plot_induc <- ggplot(ICpref3, aes(x = Range, y = ProportionEaten)) +
    stat_summary(fun.data = plot.mean, geom = "errorbar", width = 0.25) +
    stat_summary(fun = mean, geom = "point", size = 2, color = "black", shape = 16, position = position_dodge(width = 0.35)) + # Adds the mean point
    guides(size = "none") +
    theme_bw() +
    facet_wrap(~Caterpillar) +
    labs(
      title = "Feeding Preference for Induced Plants",
      y = "Proportion Consumed",
      x = "",
      color = "Defense"
    ) +
    scale_x_discrete(labels = c("Invasive" = "Familiar", "Native" = "Unfamiliar")) +
    facet_wrap(~ Caterpillar, labeller = as_labeller(c("soybean" = "Soybean Looper", "velvetbean" = "Velvetbean Caterpillar"))) +
    ylim(0.5, 0.7) + # Corrected ylim placement
    theme_bw() +
    scale_colour_grey() +
    theme(strip.text.x = element_text(color = "black", size = 16, face = "bold")) +
    theme(legend.position = "none") +
    theme(strip.text = element_text(size = 16, face = "bold")) +
    theme(axis.text = element_text(size = 14)) +
    theme(axis.title = element_text(size = 16)) +
    theme(plot.title = element_text(size = 18)) +
    theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

#panel of cafeteria plots: cafeteria_plot_const, cafeteria_plot_induc, and ICpref_proportion_caterpillar
(cafeteria_panel <- cafeteria_plot_const / cafeteria_plot_induc / ICpref_proportion_caterpillar + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16)))

ggsave("cafeteria_panel_oct24.jpg", cafeteria_panel, width = 2900, height = 3200, units = "px", dpi = 300)
## save in E and E format. 600 dpi PDF with a min 1800 px canvas size
ggsave("Figure2_cafeteria_panel_oct24.pdf", cafeteria_panel, width = 5200, height = 6000, units = "px", dpi = 600)

old_IC_prop.lmer <- lmer(Proportion ~ Range * Herbivore * Defense + (1 | Herbivore / Pair), data = ConstVInduc2)
old_IC_prop2.lm <- lm(Proportion ~ Range * Herbivore * Defense, data = ConstVInduc2)

# summary of the model
summary(old_IC_prop.lmer)

Anova(old_IC_prop.lmer)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# Response: Proportion
# Chisq Df Pr(>Chisq)
# Range                    0.0000  1  1.0000000
# Herbivore                0.0000  1  1.0000000
# Defense                  0.4780  1  0.4893117
# Range:Herbivore          0.0000  1  1.0000000
# Range:Defense           12.2607  1  0.0004626 *** # look at this interaction term more below
#   Herbivore:Defense        0.1337  1  0.7146619
# Range:Herbivore:Defense  0.8946  1  0.3442230
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(old_IC_prop.lmer, pairwise ~ Defense | Range)
emm2 <- emmeans(old_IC_prop.lmer, specs = pairwise ~ f1 | f2, type = "response")
emm2

tidy(joint_tests(old_IC_prop.lmer, by = "Defense"))
# # A tibble: 6 × 6
# term            Defense num.df  den.df statistic p.value
# <chr>           <chr>    <dbl>   <dbl>     <dbl>   <dbl>
#   1 Range           Const        1     152     6.13   0.0144
# 2 Range           Induced      1     152     6.13   0.0144
# 3 Herbivore       Const        1 3939992     0      0.984
# 4 Herbivore       Induced      1 3939992     0      0.984
# 5 Range:Herbivore Const        1     152     0.447  0.505
# 6 Range:Herbivore Induced      1     152     0.447  0.505

tidy(joint_tests(old_IC_prop.lmer, by = "Range")) #
# # A tibble: 6 × 6
# term              Range    num.df  den.df statistic p.value
# <chr>             <chr>     <dbl>   <dbl>     <dbl>   <dbl>
#   1 Herbivore         Invasive      1 1969996     0     1
# 2 Herbivore         Native        1 1969996     0     1
# 3 Defense           Invasive      1      76     8.79  0.00404 #difference between ranges in defense induction? or constitutive....
# 4 Defense           Native        1      76     3.95  0.0505 # these results are expected
# 5 Herbivore:Defense Invasive      1      76     0.86  0.357 #no diff bw herbs in range
# 6 Herbivore:Defense Native        1      76     0.168 0.683

tidy(joint_tests(old_IC_prop.lmer, by = "Herbivore"))
# # A tibble: 6 × 6
# term          Herbivore  num.df den.df statistic p.value
# <chr>         <chr>       <dbl>  <dbl>     <dbl>   <dbl>
#   1 Range         Soybean         1     76     0     1.00
# 2 Range         Velvetbean      1     76     0     1
# 3 Defense       Soybean         1     76     0.053 0.818
# 4 Defense       Velvetbean      1     76     0.559 0.457
# 5 Range:Defense Soybean         1     76     3.27  0.0747 #marginal in the interaction --
# 6 Range:Defense Velvetbean      1     76     9.89  0.00237 # strong difference in the interaction for velvetbean: to the extent we see significant difference in the interaction, it seems to be more driven by the velvetbean caterpillar than the looper

leveneTest((Proportion) ~ Range, data = subset(ConstVInduc2, Defense == "Const")) #  checking the variances


leveneTest((Proportion) ~ Range, data = subset(ConstVInduc2, Defense == "Induced"))

## ======== 2022 cafeteria analysis ========
# create win-loss column in binary
# winloss column created by subtracting the
ICpref$win_loss <- ifelse(ICpref$WinnerTest > 0, 1, 0)

# select columns needed for analysis in ICpref data (all induced data)
ICpref2 <- ICpref %>%
  select(Plate_number, Invasive_accession, Native_accession, Pair_number, Caterpillar, Replicate, prop_wtmg_Invasive, Prop_wtmg_Native, win_loss)

ICpref2a <- ICpref2 %>%
  select(Plate_number, Invasive_accession, Native_accession, Pair_number, Caterpillar, Replicate, win_loss)

# names_pattern = "(.*)_accession" removes the "_accession" bit from any matching pattern
ICpref3a <- ICpref2a %>%
  pivot_longer(cols = c(Invasive_accession, Native_accession), names_to = "Range", values_to = "Accession", names_pattern = "(.*)_accession")

ICpref2b <- ICpref2 %>%
  select(Plate_number, Pair_number, Caterpillar, Replicate, prop_wtmg_Invasive, Prop_wtmg_Native, win_loss)

ICpref3b <- ICpref2b %>%
  pivot_longer(cols = c(prop_wtmg_Invasive, Prop_wtmg_Native), names_to = "Range", names_prefix = "[Pp](.*)_wtmg_", values_to = "ProportionEaten")

ICpref3 <- full_join(ICpref3a, ICpref3b)

# change all values of 1 (100% tissue eaten) to 0.9999999 for analysis
ICpref3 <- ICpref3 %>%
  mutate(ProportionEaten = ifelse(ProportionEaten > 0.9999, 0.9999999, ProportionEaten))

ICpref3 <- ICpref3 %>%
  filter(ProportionEaten > 0)

ICpref3$Replicate <- as.factor(ICpref3$Replicate)

# summarise data. Group_by so we know for each caterpillar and range. Summarise all columns with numeric data and the functions performed are mean and std.error
ICpref4 <- ICpref3 %>%
  group_by(Caterpillar, Range) %>%
  summarise(across(
    .cols = where(is.numeric),
    .fns = list(Mean = ~ mean(., na.rm = TRUE), SE = ~ std.error(., na.rm = TRUE)),
    .names = "{col}_{fn}"
  ))

# Custom labelling function
custom_labels <- function(x) {
  x <- as.character(x) # Ensure the levels are character strings
  x[x == "soybean"] <- "Soybean Looper" # Change 'soybean' to 'Soybean Looper'
  x[x == "velvetbean"] <- "Velvetbean Caterpillar" # Change 'velvetbean' to 'Velvetbean Caterpillar'
  # Custom labelling function
  x[x == "Invasive"] <- "Familiar" # Change 'Invasive' to 'Familiar'
  x[x == "Native"] <- "Unfamiliar" # Change 'Native' to 'Unfamiliar'
  return(x)
}

ICpref_wl.glmer <- glmer(win_loss ~ Range * Caterpillar + (1 | Pair_number), family = binomial(link = logit), ICpref3)

Anova(ICpref_wl.glmer)
# Analysis of Deviance Table (Type II Wald chisquare tests)
#
# Response: win_loss
# Chisq Df Pr(>Chisq)
# Range             0.0021  1     0.9638
# Caterpillar       0.7317  1     0.3923
# Range:Caterpillar 0.0010  1     0.9742

ICpref_prop.lmer <- lmer(ProportionEaten ~ Range * Caterpillar + (1 | Pair_number), data = ICpref3)

Anova(ICpref_prop.lmer)
# Analysis of Deviance Table (Type II Wald chisquare tests)
#
# Response: ProportionEaten
# Chisq Df Pr(>Chisq)
# Range              3.0604  1    0.08022 .  # range is marginally significant
# Caterpillar       26.3462  1  2.854e-07 *** #caterpillar type is significant for amount eaten
#   Range:Caterpillar  1.3939  1    0.23775
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ======== Protein Data analysis ========
## ======= First ========
ProteinDS$Hour <- factor(ProteinDS$Hour, levels = c("0hr", "4hr", "24hr"))
ProteinDS$Defense <- ifelse(ProteinDS$Time == 0, "Constitutive", "Induced")
ProteinDS$Range <- as.factor(ProteinDS$Range)
ProteinDS$Defense <- as.factor(ProteinDS$Defense)
ProteinDS$Site <- as.factor(ProteinDS$Site)

#comparing 0hr (constitutive) with 24hr (induced)
Pro.mod <- lmer(log(Protein) ~ Range * Defense + (1|Genotype), ProteinDS)
anova(Pro.mod)
# Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Range          1.2321  1.2321     1  93.68  3.5646   0.06212 .  #marginally significant diffs in protein induction by range
# Defense       22.7415 22.7415     1 757.93 65.7914 2.015e-15 ***
# Range:Defense  0.0541  0.0541     1 757.93  0.1566   0.69239    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#do comparison by hour instead of range 
Pro_hour.mod <- lmer(log(Protein) ~ Range * Hour + (1|Genotype), ProteinDS)
anova(Pro_hour.mod)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Range       1.1584  1.1584     1  89.99  3.4261   0.06746 .  
# Hour       24.4193 12.2097     2 755.82 36.1122 1.051e-15 ***
# Range:Hour  0.6432  0.3216     2 755.82  0.9512   0.38673    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

emmeans(Pro.mod, pairwise ~ Defense | Range)
# $emmeans
# Range = Invasive:
#   Defense      emmean     SE  df lower.CL upper.CL
# Constitutive   11.4 0.1400 157     11.1     11.7
# Induced        10.9 0.1265 106     10.7     11.2
# 
# Range = Native:
#   Defense      emmean     SE  df lower.CL upper.CL
# Constitutive   11.1 0.0641 155     11.0     11.3
# Induced        10.7 0.0579 104     10.6     10.8
# 
# Degrees-of-freedom method: kenward-roger 
# Results are given on the log (not the response) scale. 
# Confidence level used: 0.95 
# 
# $contrasts
# Range = Invasive:
#   contrast               estimate     SE  df t.ratio p.value
# Constitutive - Induced    0.485 0.1039 758   4.669  <.0001
# 
# Range = Native:
#   contrast               estimate     SE  df t.ratio p.value
# Constitutive - Induced    0.440 0.0471 759   9.354  <.0001
# 
# Degrees-of-freedom method: kenward-roger 
# Results are given on the log (not the response) scale. 

# Define a dodge position object with a specific width to dodge by
dodge <- position_dodge(width = 0.1)

tidy(joint_tests(Pro.mod, by = "Defense"))
# A tibble: 2 × 6
# term  Defense      num.df den.df statistic p.value
# <chr> <chr>         <dbl>  <dbl>     <dbl>   <dbl>
#   1 Range Constitutive      1   156.      3.26  0.0731
# 2 Range Induced           1   105.      2.80  0.0974

tidy(joint_tests(Pro.mod, by = "Range"))
# # A tibble: 2 × 6
# term    Range    num.df den.df statistic  p.value
# <chr>   <chr>     <dbl>  <dbl>     <dbl>    <dbl>
#   1 Defense Invasive      1   758.      21.8 3.58e- 6
# 2 Defense Native        1   759.      87.5 9.27e-20

leveneTest(Protein ~ Range, data = subset(ProteinDS, Defense == "Constitutive")) #  checking the variances

leveneTest(log(Protein) ~ Range, data = subset(ProteinDS, Defense == "Constitutive")) #  checking the variances
# Levene's Test for Homogeneity of Variance (center = median)
#        Df F value Pr(>F)
# group   1   0.899 0.3439
#       280    

Pro.mod2 <- lmer(log(Protein) ~ Range + Defense + (1|Genotype), ProteinDS) # checking alternative model -- without interaction winner
anova(Pro.mod2)
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF  DenDF  F value  Pr(>F)    
# Range    1.183   1.183     1  89.99   3.4252 0.06749 .  # marginally significant diffs in protein induction by range
# Defense 37.726  37.726     1 759.81 109.2637 < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(Pro.mod, Pro.mod2) # comparing the two models
# refitting model(s) with ML (instead of REML)
# Data: ProteinDS
# Models:
#   Pro.mod2: log(Protein) ~ Range + Defense + (1 | Genotype)
# Pro.mod: log(Protein) ~ Range * Defense + (1 | Genotype)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# Pro.mod2    5 1687.7 1711.5 -838.87   1677.7                     
# Pro.mod     6 1689.6 1718.1 -838.79   1677.6 0.1569  1      0.692


## ======== POD analysis from 2018 data ======
PODDS$Hour <- factor(PODDS$Hour, levels = c("0hr", "4hr", "24hr"))
PODDS$Defense <- ifelse(PODDS$Time == 0, "Constitutive", "Induced")
PODDS$Range <- as.factor(PODDS$Range)
PODDS$Defense <- as.factor(PODDS$Defense)
PODDS$Site <- as.factor(PODDS$Site)
PODDS$TransfPOD <- log(PODDS$AbsFreshWeight + 4)

# with interaction -- winning model
POD_hour.mod <- lmer(TransfPOD ~ Range * Hour + (1 | Genotype), PODDS)
anova(POD_hour.mod)
# Type III Analysis of Variance Table with Satterthwaite's method
#             Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)
# Range        0.057   0.057     1  89.81   0.322 0.57185
# Hour       120.923  60.462     2 755.66 339.589 < 2e-16 ***
# Range:Hour   0.841   0.421     2 755.66   2.363 0.09484 .  # marginal significance in the interaction
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tidy(joint_tests(POD_hour.mod, by = "Hour"))
# # A tibble: 3 × 6
# term  Hour  num.df den.df statistic p.value
# <chr> <chr>  <dbl>  <dbl>     <dbl>   <dbl>
#   1 Range 0hr        1   150.     0       0.987
# 2 Range 4hr        1   150.     2.36    0.126
# 3 Range 24hr       1   150.     0.001   0.979

## ======== PPO analysis from 2018 ======
PPODS$Hour <- factor(PPODS$Hour, levels = c("0hr", "4hr", "24hr"))
PPODS$Defense <- ifelse(PPODS$Time == 0, "Constitutive", "Induced")
PPODS$Range <- as.factor(PPODS$Range)
PPODS$Defense <- as.factor(PPODS$Defense)
PPODS$Site <- as.factor(PPODS$Site)
PPODS$TransfPPO <- log(PPODS$AbsFreshWeight + 1)
PPODS <- na.omit(PPODS)

PPO.mod <- lmer(TransfPPO ~ Range * Defense + (1 | Genotype), PPODS)
anova(PPO.mod)
# Type III Analysis of Variance Table with Satterthwaite's method
#               Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)
# Range         0.4215  0.4215     1  95.58  2.8049   0.09724 .  #marginal sig range
# Defense       3.5785  3.5785     1 757.68 23.8150 1.293e-06 ***
# Range:Defense 0.0609  0.0609     1 757.68  0.4053   0.52456
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tidy(joint_tests(PPO.mod, by = "Defense"))
# # A tibble: 2 × 6
# term  Defense      num.df den.df statistic p.value
# <chr> <chr>         <dbl>  <dbl>     <dbl>   <dbl>
#   1 Range Constitutive      1   198.      2.77  0.0978 #marginal significance maybe driven by constit levels?
# 2 Range Induced           1   114.      1.65  0.202


tidy(joint_tests(PPO.mod, by = "Range"))
# A tibble: 2 x 6
# term    Range    num.df den.df statistic     p.value
# <chr>   <chr>     <dbl>  <dbl>     <dbl>       <dbl>
#   1 Defense Invasive      1   758.      9.17 0.00255
# 2 Defense Native        1   760.     26.5  0.000000340

#---- Figure 4: Protein, POD, and PPO induction ----

#facet for Protein induction
(Pro_plot <- ggplot(ProteinDS, aes(x = Hour, y = Protein, fill = Hour)) + 
   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) + 
   stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "black") +  # Change fill to black
   facet_wrap(~Range, labeller = labeller(Range = c(Invasive = "Familiar", Native = "Unfamiliar"))) + 
   labs(title = "Protein Content Induction", y = "Protein (mg/g FW)") + 
   theme_bw() + 
   scale_colour_grey() + 
   theme(strip.text.x = element_text(color = "black", size = 16, face = "bold")) + 
   theme(legend.position = "none") + 
   theme(strip.text.x = element_text(size = 16, face = "bold")) + 
   theme(axis.text = element_text(size = 14)) + 
   theme(axis.title = element_text(size = 16)) + 
   theme(axis.title.x = element_blank()) +
   theme(plot.title = element_text(size = 18)) + 
   theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

# facet for POD induction
(POD_plot <- (PODfig <- ggplot(PODDS, aes(x = Hour, y = TransfPOD, fill = Hour)) +
                stat_summary(fun.data = plot.mean, geom = "errorbar", width = 0.25) + 
                stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "black") +  # Add this line to plot means
                facet_wrap(~Range, labeller = labeller(Range = c(Invasive = "Familiar", Native = "Unfamiliar"))) +
                labs(title = "Peroxidase Induction", x = NULL, y = "Absorbance / g FW") +
                theme_bw() +
                scale_colour_grey() +
                theme(strip.text.x = element_text(color = "black", size = 16, face = "bold")) +
                theme(legend.position = "none") +
                theme(strip.text.x = element_text(size = 16, face = "bold")) +
                theme(axis.text = element_text(size = 14)) +
                theme(axis.title = element_text(size = 16)) +
                theme(plot.title = element_text(size = 18)) +
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())))

# facet for PPO induction
(PPO_plot <- (PPOfig <- ggplot(PPODS, aes(x = Hour, y = TransfPPO, fill = Hour)) +
                stat_summary(fun.data = plot.mean, geom = "errorbar", width = 0.25) +
                stat_summary(fun = mean, geom = "point", shape = 21, size = 2, color = "black", fill = "black") +  # Add this line to plot means
                facet_wrap(~Range, labeller = labeller(Range = c(Invasive = "Familiar", Native = "Unfamiliar"))) +
                labs(title = "Polyphenol Oxidase Induction", x = NULL, y = "Absorbance / g FW") +
                theme_bw() +
                scale_colour_grey() +
                theme(strip.text.x = element_text(color = "black", size = 16, face = "bold")) +
                theme(legend.position = "none") +
                theme(strip.text = element_text(size = 16, face = "bold")) +
                theme(axis.text = element_text(size = 14)) +
                theme(axis.title = element_text(size = 16)) +
                theme(plot.title = element_text(size = 18)) +
                theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank())))

(figure4 <- Pro_plot / POD_plot / PPO_plot + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 18)))
ggsave("AssayFig14august24_all_bioc.jpg", figure3, width = 2900, height = 3200, units = "px", dpi = 300)

# save in E and E preferred format. 600 dpi PDF with an 1800 px canvas size 
ggsave("Figure4_assayFig14august24_all_bioc.pdf", figure3, width = 5200, height = 6000, units = "px", dpi = 600)

# ========= Correlation analysis ========



## ======== Inducibility & Correlations data prep ========
PODtrade <- read.csv("data/processed/WSU_PODfiles11Sept2017.csv")

# remove unwanted columns:
PODtrade[, c(1, 3, 6:13)] <- list(NULL)

# remove hour 4
PODtrade <- subset(PODtrade, Time != 4)

# Giving appropriate column names
PODtrade$trt <- ifelse(PODtrade$Time == 0, 1, 2)
PODtrade$res <- log(PODtrade$AbsFreshWeight + 4)

# Now ready to set up dataset according to script

# A: Both ranges together
# 1) Genotype is level of permutation

PODtradeAll_Gen <- PODtrade
# change name of replicate column to "rep"
PODtradeAll_Gen <- rename(PODtradeAll_Gen, rep = Replicate)

uniqueGeno <- unique(PODtradeAll_Gen$Sample)
GenDF <- data.frame(uniqueGeno, fam = 1:length(uniqueGeno))

PODtradeAll_Gen$fam <- GenDF[match(PODtradeAll_Gen$Sample, GenDF$uniqueGeno), 2]
unique(PODtradeAll_Gen$fam)
# family 95 (W0477) does not have time 0 data
PODtradeAll_Gen <- subset(PODtradeAll_Gen, fam != 95)
str(PODtradeAll_Gen)
PODtradeAll_Gen[, c(1, 3:7)] <- list(NULL)
write.table(PODtradeAll_Gen, file = "POD_spur_All_Geno.csv", sep = ",", row.names = F)

# # From running the permute code, we have successfully detected a tradeoff! (see code beginning at line _____)
### ==== B: Invasive range =====
# 1) Genotype is level of permutation

PODtrade_I_Gen <- subset(PODtrade, Range == "Invasive")

uniqueGeno <- unique(PODtrade_I_Gen$Sample)
GenDF <- data.frame(uniqueGeno, fam = 1:length(uniqueGeno))

PODtrade_I_Gen$fam <- GenDF[match(PODtrade_I_Gen$Sample, GenDF$uniqueGeno), 2]
table(PODtrade_I_Gen$fam)

PODtrade_I_Gen[, c(1, 3:7)] <- list(NULL)
write.table(PODtrade_I_Gen, file = "POD_spur_Inv_Geno.csv", sep = ",", row.names = F)
# Negative tradeoff is real!


### ==== B: Native range =====
# 1) Genotype is level of permutation

PODtrade_N_Gen <- subset(PODtrade, Range != "Invasive")

uniqueGeno <- unique(PODtrade_N_Gen$Sample)
GenDF <- data.frame(uniqueGeno, fam = 1:length(uniqueGeno))

PODtrade_N_Gen$fam <- GenDF[match(PODtrade_N_Gen$Sample, GenDF$uniqueGeno), 2]
table(PODtrade_N_Gen$fam)
# family 79 (W0477) does not have time 0 data
PODtrade_N_Gen <- subset(PODtrade_N_Gen, fam != 79)

PODtrade_N_Gen[, c(1, 3:7)] <- list(NULL)
write.table(PODtrade_N_Gen, file = "POD_spur_Nat_Geno.csv", sep = ",", row.names = F)

# Negative tradeoff is real!

## ==== PPO data prep ========
# PPO correlations using WSU PPOfiles11Sept2017
PPOtrade <- read.csv("data/processed/WSU_PPOfiles11Sept2017.csv")

# remove unwanted columns:
PPOtrade[, c(1, 3, 6:13)] <- list(NULL)

# remove hour 4
PPOtrade <- subset(PPOtrade, Time != 4)

# Giving appropriate column names
PPOtrade$trt <- ifelse(PPOtrade$Time == 0, 1, 2)
PPOtrade$res <- log(PPOtrade$AbsFreshWeight + 4)

# Now ready to set up dataset according to script

# A: Both ranges together
# 1) Genotype is level of permutation

PPOtradeAll_Gen <- PPOtrade
PPOtradeAll_Gen <- rename(PPOtradeAll_Gen, rep = Replicate)

uniqueGeno <- unique(PPOtradeAll_Gen$Sample)
GenDF <- data.frame(uniqueGeno, fam = 1:length(uniqueGeno))

PPOtradeAll_Gen$fam <- GenDF[match(PPOtradeAll_Gen$Sample, GenDF$uniqueGeno), 2]
table(PPOtradeAll_Gen$fam)
# family 95 (W0477) does not have time 0 data
PPOtradeAll_Gen <- subset(PPOtradeAll_Gen, fam != 95)
str(PPOtradeAll_Gen)
PPOtradeAll_Gen[, c(1, 3:7)] <- list(NULL)
write.table(PPOtradeAll_Gen, file = "PPO_spur_All_Geno.csv", sep = ",", row.names = F)

### ==== B: Invasive range =====
# 1) Genotype is level of permutation

PPOtrade_I_Gen <- subset(PPOtrade, Range == "Invasive")

uniqueGeno <- unique(PPOtrade_I_Gen$Sample)
GenDF <- data.frame(uniqueGeno, fam = 1:length(uniqueGeno))

PPOtrade_I_Gen$fam <- GenDF[match(PPOtrade_I_Gen$Sample, GenDF$uniqueGeno), 2]
table(PPOtrade_I_Gen$fam)

PPOtrade_I_Gen[, c(1, 3:7)] <- list(NULL)
write.table(PPOtrade_I_Gen, file = "PPO_spur_Inv_Geno.csv", sep = ",", row.names = F)
# seems marginal but probably real

### ==== B: Native range =====
# 1) Genotype is level of permutation

PPOtrade_N_Gen <- subset(PPOtrade, Range != "Invasive")

uniqueGeno <- unique(PPOtrade_N_Gen$Sample)
GenDF <- data.frame(uniqueGeno, fam = 1:length(uniqueGeno))

PPOtrade_N_Gen$fam <- GenDF[match(PPOtrade_N_Gen$Sample, GenDF$uniqueGeno), 2]
table(PPOtrade_N_Gen$fam)
# family 79 (W0477) does not have time 0 data
PPOtrade_N_Gen <- subset(PPOtrade_N_Gen, fam != 79)

PPOtrade_N_Gen[, c(1, 3:7)] <- list(NULL)
write.table(PPOtrade_N_Gen, file = "PPO_spur_Nat_Geno.csv", sep = ",", row.names = F)
# negative tradeoff is real!

## ==== Protein data prep ========
# Protein correlations using WSU Profiles11Sept2017
Protrade <- read.csv("data/processed/ProteinData22Feb2017.csv")
# make new column for defense type based on value of hour column using case_when and mutate

Protrade <- Protrade %>%
  mutate(Defense = case_when(
    Hour == "0hr" ~ "Constitutive",
    Hour == "4hr" ~ NA,
    Hour == "24hr" ~ "Induced"
  ))

# columns to keep: Sample, Replicate, Range, Defense, Hour, Protein

# remove unwanted columns: 3,7:9,12,13,14,15

Protrade[, c(1, 3, 4, 7:10, 12:15)] <- list(NULL)

# remove hour 4
Protrade <- subset(Protrade, Hour != "4hr")

# Now ready to set up dataset according to script

# A: Both ranges together
# 1) Genotype is level of permutation

ProtradeAll_Gen <- Protrade
ProtradeAll_Gen <- rename(ProtradeAll_Gen, rep = Replicate)

uniqueGeno <- unique(ProtradeAll_Gen$Sample)
GenDF <- data.frame(uniqueGeno, fam = 1:length(uniqueGeno))

ProtradeAll_Gen$fam <- GenDF[match(ProtradeAll_Gen$Sample, GenDF$uniqueGeno), 2]
table(ProtradeAll_Gen$fam)
# family 95 (W0477) does not have time 0 data
ProtradeAll_Gen <- subset(ProtradeAll_Gen, fam != 95)
head(ProtradeAll_Gen)
ProtradeAll_Gen[, -c(3, 8:11)] <- list(NULL)
write.table(ProtradeAll_Gen, file = "Pro_spur_All_Geno.csv", sep = ",", row.names = F)

# From running the permute code, we have successfully detected a tradeoff!


# ##==== B: Invasive range =====
# 1) Genotype is level of permutation

Protrade_I_Gen <- subset(Protrade, Range == "Invasive")
head(Protrade_I_Gen)

uniqueGeno <- unique(Protrade_I_Gen$Sample)
GenDF <- data.frame(uniqueGeno, fam = 1:length(uniqueGeno))

Protrade_I_Gen$fam <- GenDF[match(Protrade_I_Gen$Sample, GenDF$uniqueGeno), 2]
table(Protrade_I_Gen$fam)

Protrade_I_Gen[, -c(3, 8:10)] <- list(NULL)
write.table(Protrade_I_Gen, file = "Pro_spur_Inv_Geno.csv", sep = ",", row.names = F)

# It's real!
# ##==== B: Native range =====
# 1) Genotype is level of permutation

Protrade_N_Gen <- subset(Protrade, Range != "Invasive")

uniqueGeno <- unique(Protrade_N_Gen$Sample)
GenDF <- data.frame(uniqueGeno, fam = 1:length(uniqueGeno))

Protrade_N_Gen$fam <- GenDF[match(Protrade_N_Gen$Sample, GenDF$uniqueGeno), 2]
table(Protrade_N_Gen$fam)
# family 79 (W0477) does not have time 0 data
Protrade_N_Gen <- subset(Protrade_N_Gen, fam != 79)
head(Protrade_N_Gen)
Protrade_N_Gen[, -c(3, 8:10)] <- list(NULL)
write.table(Protrade_N_Gen, file = "Pro_spur_Nat_Geno.csv", sep = ",", row.names = F)



## ======== Inducibility & Correlations analysis ========
# Randomization test for non-random association between measures of constitutive and induced resistance

# Uses two measures of association:
# 1) correlation between induced resistance (damaged minus control family means) and constitutive
#    resistance (control family mean)
# 2) slope of a regression of damaged vs. control family means

# See Morris, W.F., et al. Oikos (2006) 112: 102-110.

###########################  DATA FILE REQUIREMENTS  ###################################
# 														   #
# The data file must be a csv (comma-separated values) file with 4 columns,		   #
# all of which must contain numeric data:								   #
# 1) family, genotype, or species index, going from 1 to the number of groups;	   #
# 		numbers should not be skipped								   #
# 2) treatment: 1 = control 2 = damaged treatment						   #
# 3) replicate plant number										   #
# 4) measure of resistance (assumed to be a non-negative number)				   #
# 														   #
# Column headings in datafile must be "fam", "trt", "rep", "res" (without the quotes)  #
# 														   #
########################################################################################
# This script has a toggle (see "gonormal" below) for two options in cases
# in which negative resistances are likely to arise in the randomization
# procedure.

# Written by W.F. Morris
# Last revised: 6 Dec 2009

rm(list = ls(all = TRUE)) # clear memory
graphics.off() # close open figures

# *************************** USER DEFINED PARAMETERS *******************************

# number of bootstrap replicates to perform
maxreps <- 10000

# number of families, genotypes, or species in the data set
numfams <- 94

# Chose one of the following options for what to do if, for a family, Imean, the randomly
# generated mean of induced resistance, is negative and less in absolute value than Cmin,
# the lowest randomly generated constitutive resistance (which would frequently lead to
# negative resistances in the randomization):
#   gonormal = 1 -> generate normally distributed induced resistances but
#       discard values less than -Cmin
#   gonormal = 0 -> try generating new family mean constitutive and induced
#       resistances until Imean > -Cmin; this will be slower.
gonormal <- 1

# name of comma-separated values file containing the data
fname <- "POD_spur_All_Geno.csv"

# ***********************************************************************************


q <- numfams - 1
tcrit <- qt(0.975, numfams) # critical value of t statistic with numfams degrees of freedom
chi2critlo <- qchisq(.975, numfams - 1) # critical lower value of chi-square statistic with numfams-1 degrees of freedom
chi2critup <- qchisq(.025, numfams - 1) # critical upper value of chi-square statistic with numfams-1 degrees of freedom

data <- read.csv(fname)
attach(data) # column headings in datafile must be fam, trt, rep, res

Nc <- Nd <- Cest <- Dest <- Cc <- Cd <- matrix(0, numfams, 1)

for (f in 1:numfams) {
  I <- which(fam == f & trt == 1) # find control plants for fam f
  Nc[f] <- length(I)
  x <- res[I]
  Cest[f] <- mean(x) # find the fam mean resistance
  Cc[f] <- sd(x) / Cest[f] # and coefficient of variation

  I <- which(fam == f & trt == 2) # repeat for damage trt.
  Nd[f] <- length(I)
  x <- res[I]
  Dest[f] <- mean(x)
  Cd[f] <- sd(x) / Dest[f]
}

Iest <- Dest - Cest # difference measure of induced resistance
CVs <- (Nc * Cc + Nd * Cd) / (Nc + Nd) # weighted mean CV across trts., by fam
CV <- mean(CVs) # mean weighted CV across families
sigma2 <- log(1 + CV^2)
sigma <- sigma2^.5 # SD of exponent in lognormal resistance

# Plot the data to make sure it looks correct: comparing first graph with second to make sure the second looks like the first. The procedure can't take negative results so a monte carlo resampling procedure is used to generate a distribution of possible results.

split.screen(matrix(c(.25, .75, .5, 1, .25, .75, 0, .5), 2, 4, byrow = TRUE))

screen(1)
xx <- seq(0, 2, .01)
plot(Cest, Iest,
  type = "p",
  xlab = "Constitutive Resistance",
  ylab = "Induced Resistance",
  xlim = c(min(Cest) - .01, max(Cest) + .01),
  ylim = c(min(Iest) - .01, max(Iest) + .01)
)
lines(c(min(Cest) - .01, max(Cest) + .01), c(0, 0), lty = "dashed")
lines(c(0, 0), c(min(Iest) - .01, max(Iest) + .01), lty = "dashed")

screen(2)
xx <- seq(min(min(Dest), min(Cest)), max(max(Dest), max(Cest)), .1)
plot(Cest, Dest,
  type = "p",
  xlab = "Control Resistance",
  ylab = "Damaged Resistance",
  xlim = c(min(Cest) - .01, max(Cest) + .01),
  ylim = c(min(Dest) - .01, max(Dest) + .01)
)
lines(xx, xx)

# compute confidence interval for population mean and among-fam variance in constitutive resistance
Cmeanhat <- mean(Cest)
Cvarhat <- var(Cest)
SE.Cmean <- (Cvarhat / numfams)^.5
CI.Cmean <- c(Cmeanhat - tcrit * SE.Cmean, Cmeanhat + tcrit * SE.Cmean)
CI.Cvar <- c(Cvarhat * q / chi2critlo, Cvarhat * q / chi2critup)

# compute confidence interval for population mean and among-fam variance in induced resistance
Imeanhat <- mean(Iest)
Ivarhat <- var(Iest)
SE.Imean <- (Ivarhat / numfams)^.5
CI.Imean <- c(Imeanhat - tcrit * SE.Imean, Imeanhat + tcrit * SE.Imean)
CI.Ivar <- c(Ivarhat * q / chi2critlo, Ivarhat * q / chi2critup)

# observed tradeoff measures
corrobs <- cor(Iest, Cest) # correlation measure of tradoff
bobs <- cov(Dest, Cest) / var(Cest) # slope measure of tradeoff

print("Observed correlation between induced (damage minus control)and constitutive (control) resistance", quote = FALSE)
corrobs
# [,1]
# [1,] -0.6195606
print("Observed slope of regression of damaged vs. control fam means", quote = FALSE)
bobs
# [,1]
# [1,] 0.06584537

correst <- best <- matrix(0, 1, maxreps)

# generate pseudodata with no correlation between C and I

rep <- 0

while (rep < maxreps) {
  done <- 0

  # generate grand fam mean and var for C and I using their confidence intervals

  Cmean <- Inf
  while (Cmean < CI.Cmean[1] || Cmean > CI.Cmean[2]) {
    Cmean <- Cmeanhat + SE.Cmean * rnorm(1)
  }

  Cvar <- Inf
  while (Cvar < CI.Cvar[1] || Cvar > CI.Cvar[2]) {
    Cvar <- Cvarhat * rchisq(1, df = q) / q
  }

  Imean <- Inf
  while (Imean < CI.Imean[1] || Imean > CI.Imean[2]) {
    Imean <- Imeanhat + SE.Imean * rnorm(1)
  }

  Ivar <- Inf
  while (Ivar < CI.Ivar[1] || Ivar > CI.Ivar[2]) {
    Ivar <- Ivarhat * rchisq(1, df = q) / q
  }
  Isd <- Ivar^.5

  # generate (lognormal) resistance measures for each fam in control and damage treatments
  sig2 <- log(1 + Cvar / Cmean^2)
  C <- exp(rnorm(numfams, log(Cmean) - 0.5 * sig2, sig2^.5)) # generate lognormal C's
  Cmin <- min(C)
  if (Imean > -Cmin) { # abs(Imean)>Cmin, so mean of the following lognormal is positive
    done <- 1
    rep <- rep + 1
    sig2 <- log(1 + Ivar / (Imean + Cmin)^2)
    I <- exp(rnorm(numfams, log(Imean + Cmin) - 0.5 * sig2, sig2^.5)) - Cmin # generate I's as shifted lognormal, with Imin=-Cmin
  } else {
    if (gonormal) { # Imean<=-Cmin, so mean of preceding lognormal would be neg; instead draw normal I's and discard if < -Cmin
      done <- 1
      rep <- rep + 1
      nf <- 0
      while (nf < numfams) {
        Itemp <- rnorm(1, Imean, Isd)
        if (Itemp > -Cmin) {
          nf <- nf + 1
          I[nf] <- Itemp
        }
      }
    }
  }
  # NOTE: The preceding "else" statement will cause the mean and
  # variance of I to differ from Imean and Ivar

  if (done) {
    D <- C + I
    Cest <- Dest <- matrix(0, numfams, 1)
    for (f in 1:numfams) { # generate lognormal resistance measures for each plant, and take means
      Cest[f] <- mean(exp(rnorm(Nc[f], log(C[f]) - .5 * sigma2, sigma))) # use of sigma2 and sigma assumes fixed CV across families and trts.
      Dest[f] <- mean(exp(rnorm(Nd[f], log(D[f]) - .5 * sigma2, sigma)))
    }
    Iest <- Dest - Cest
    correst[rep] <- cor(Iest, Cest)
    best[rep] <- cov(Dest, Cest) / var(Cest) # slope of regression of Dest on Cest
  }
} # while(rep<maxreps)

corrlow <- quantile(correst, .05)
slopelow <- quantile(best, .05)

print("Lower 5th percentile of the bootstrap distribution of the correlation coefficient", quote = FALSE)
corrlow
# 5%
# -0.1631732
print("Lower 5th percentile of the bootstrap distribution of the regression slope", quote = FALSE)
slopelow
# 5%
# 0.8350824

# for plotting, set lower limit at largest value divisible by .2 below estimates
ccmin <- round(min(correst), 1)
if (ccmin > min(correst)) ccmin <- ccmin - .1
if (!(ccmin %% .2 < 1e-10)) ccmin <- ccmin - .1
# for plotting, set upper limit at smallest value divisible by .2 above estimates
ccmax <- round(max(correst), 1)
if (ccmax < max(correst)) ccmax <- ccmax + .1
if (!(ccmax %% .2 < 1e-10)) ccmax <- ccmax + .1
dcc <- (ccmax - ccmin) / 50
edgecc <- seq(ccmin, ccmax, dcc)
cctks <- seq(ccmin, ccmax, .2)

yy <- hist(correst, breaks = edgecc, plot = FALSE)
yc <- yy$counts
lowcc <- min(which(yc > 0))
highcc <- max(which(yc > 0))
yccmax <- 1.1 * max(yc)
Fcc <- cumsum(yc) / maxreps
cci <- min(which(Fcc >= 0.05)) # index of critical correlation coeff.
cco <- min(which(edgecc >= as.numeric(corrobs))) # index of observed correlation coeff.
pccobs <- Fcc[cco]
print("Approximate probability of observed or smaller correlation", quote = FALSE)
pccobs
# [1] 0.4473
bmin <- floor(min(best))
bmax <- ceiling(max(best))
db <- (bmax - bmin) / 50
edgeb <- seq(bmin, bmax, db)
btks <- seq(bmin, bmax, 1)

yy <- hist(best, breaks = edgeb, plot = FALSE)
lowb <- min(which(yy$counts > 0))
highb <- max(which(yy$counts > 0))
ybmax <- 1.1 * max(yy$counts)
Fb <- cumsum(yy$counts) / maxreps
bi <- min(which(Fb >= 0.05)) # index of critical regression slope
bo <- min(which(edgeb >= as.numeric(bobs))) # index of observed regression slope
pbobs <- Fb[bo]
print("Approximate probability of observed or smaller regression slope", quote = FALSE)
pbobs
# [1] 0.9994

# Plot the bootstrap distribution of the correlation coefficient;
# plot the lower 5th percentile as a vertical dashed line
# plot the observed value as a vertical solid line

# windows()
graphics.off()
hist(correst,
  breaks = edgecc, ylab = "Number of replicates", xlab = "Correlation coefficient",
  main = "Solid:Observed; Dashed:CumProb=0.05",
  xlim = c(ccmin, ccmax), ylim = c(0, yccmax), lwd = 1.5, ps = 16, cex.main = 1, xaxt = "n"
)
axis(1, cctks)
lines(edgecc[cci] * c(1, 1), c(0, yccmax), lwd = 3, lty = "dashed", col = "blue")
lines(edgecc[cco] * c(1, 1), c(0, yccmax), lwd = 2, col = "red") # result = normal (expected)


# Plot the bootstrap distributions of the correlation coefficient and regression slope,
# as well as the cumulative distributions of the two measures;
# plot the lower 5th percentile as a vertical dashed line;
# plot the observed value as a vertical solid line

# windows()
graphics.off()
split.screen(c(2, 2))

screen(1)
hist(correst,
  breaks = edgecc, ylab = "Number of replicates", xlab = "Correlation coefficient",
  main = "Solid: Observed; Dashed: CumProb=0.05",
  xlim = c(ccmin, ccmax), ylim = c(0, yccmax), lwd = 1.5, ps = 16, cex.main = .9, xaxt = "n"
)
axis(1, cctks)
lines(edgecc[cci] * c(1, 1), c(0, yccmax), lwd = 3, lty = "dashed")
lines(edgecc[cco] * c(1, 1), c(0, yccmax), lwd = 3)

screen(2)
hist(best,
  breaks = edgeb, ylab = "Number of replicates", xlab = "Regression slope",
  main = "Solid: Observed; Dashed: CumProb=0.05",
  xlim = c(bmin, bmax), ylim = c(0, ybmax), lwd = 1.5, ps = 16, cex.main = .9,
  xaxp = c(bmin, bmax, bmax - bmin)
)
lines(edgeb[bi] * c(1, 1), c(0, ybmax), lwd = 3, lty = "dashed")
lines(edgeb[bo] * c(1, 1), c(0, ybmax), lwd = 3)

screen(3)
plot(edgecc[-length(edgecc)], Fcc,
  type = "l", ylab = "Cumulative probability", xlab = "Correlation coefficient",
  main = "Horiz. line shows Pr(observed correlation)",
  xlim = c(ccmin, ccmax), ylim = c(0, 1), lwd = 1.5, ps = 16, cex.main = .9, xaxt = "n"
)
axis(1, cctks)
lines(edgecc[cco] * c(1, 1), c(0, Fcc[cco]), lwd = 2)
lines(c(ccmin, edgecc[cco]), Fcc[cco] * c(1, 1), lwd = 2)

screen(4)
plot(edgeb[-length(edgeb)], Fb,
  type = "l", ylab = "Cumulative probability", xlab = "Regression slope",
  main = "Horiz. line shows Pr(observed slope)",
  xlim = c(bmin, bmax), ylim = c(0, 1), lwd = 1.5, ps = 16, cex.main = .9, xaxt = "n"
)
axis(1, btks)
lines(edgeb[bo] * c(1, 1), c(0, Fb[bo]), lwd = 2)
lines(c(bmin, edgeb[bo]), Fb[bo] * c(1, 1), lwd = 2)



# #======== Inducibility figures ========
# load .csv files as dataframes
ProCorr <- read.csv("data/processed/Protein_Correlations9Dec2018.csv")
PodCorr <- read.csv("data/processed/POD_Correlations9Dec2018.csv")
PpoCorr <- read.csv("data/processed/PPO_Correlations9Dec2018.csv")

# Example of how to calculate inducibility (using Protein data)
# ProCorr <- within(ProCorr, DegrInduc <- (ProteinHr24 - ProteinHr0)/ProteinHr0)
# ProCorr <- within(ProCorr, {
#   Tr_Hour0 <- exp(ProteinHr0)
#   Tr_Hour24 <- exp(ProteinHr24)
#   Tr_DegInduc <- (Tr_Hour24 - Tr_Hour0)/Tr_Hour0
#   NewTr_DegInduc <- Tr_Hour24 - Tr_Hour0
#   NewDegrInduc <- ProteinHr24 - ProteinHr0 # Degree induction method 2
# })

# Protein correlations
PCnat <- subset(ProCorr, Range == "Native")
PCinv <- subset(ProCorr, Range != "Native")

cor.test(PCnat$Tr_Hour0, PCnat$NewTr_DegInduc, method = "kendall") # Tr = Transformed , NewTr = New Transformed and degree inducibility
# Kendall's rank correlation tau #comparing it within a range (below is overall, not range specific)
#
# data:  PCnat$Tr_Hour0 and PCnat$NewTr_DegInduc
# z = -12.179, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0 #correlation between time 0 and the inducibility protein (using degree of inducibility is because you need to take into account the starting point to make claims about the final point of whatever you're testing ) This is native
# sample estimates:
#        tau
# -0.5649965

cor.test(PCinv$Tr_Hour0, PCinv$NewTr_DegInduc, method = "kendall")
# Kendall's rank correlation tau
#
# data:  PCinv$Tr_Hour0 and PCinv$NewTr_DegInduc
# T = 103, p-value = 8.738e-08
# alternative hypothesis: true tau is not equal to 0 #same as above but invasive. Tau is the correlation coefficient and is negative
# sample estimates:
#        tau
# -0.6098485

cor.test(ProCorr$Tr_Hour0, ProCorr$NewTr_DegInduc, method = "kendall")
# Kendall's rank correlation tau
#
# data:  ProCorr$Tr_Hour0 and ProCorr$NewTr_DegInduc
# z = -13.01, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0 #overall
# sample estimates:
#        tau
# -0.5604476

Pro.list <- list(PCinv, PCnat)
cocor(~ ProteinHr0 + NewDegrInduc | ProteinHr0 + NewDegrInduc, Pro.list)
# Comparison between r1.jk (ProteinHr0, NewDegrInduc) = -0.5873 and r2.hm (ProteinHr0, NewDegrInduc) = -0.7691
# Difference: r1.jk - r2.hm = 0.1818
# Data: Pro.list: j = ProteinHr0, k = NewDegrInduc, h = ProteinHr0, m = NewDegrInduc
# Group sizes: n1 = 33, n2 = 210
# Null hypothesis: r1.jk is equal to r2.hm
# Alternative hypothesis: r1.jk is not equal to r2.hm (two-sided)
# Alpha: 0.05
#
# fisher1925: Fisher's z (1925)
#   z = 1.7639, p-value = 0.0777 #worth reporting
#   Null hypothesis retained #cocor looks at the correlations and looks to see if there is a difference between correlations. It's a test of differences between severity of correlations.
#
# zou2007: Zou's (2007) confidence interval
# 95% confidence interval for r1.jk - r2.hm: -0.0153 0.4679
# Null hypothesis retained (Interval includes 0)

(ProTradeFig <- ggplot(ProCorr, aes(x = ProteinHr0, y = NewDegrInduc)) +
  geom_point(shape = 21, size = 2, fill = "grey") +
  facet_wrap(~Range) +
  theme_pubr() +
  theme(strip.background = element_blank(), axis.text = element_text(size = 18), strip.text = element_text(size = 18), axis.title = element_text(size = 20), legend.text = element_text(size = 17), legend.title = element_text(size = 18)) +
  labs(y = "Inducibility", x = "Constitutive level (Protein)") +
  geom_smooth(method = "lm", se = F))

# updated corrplot
(ProTradeFig <- ggplot(ProCorr, aes(x = ProteinHr0, y = NewDegrInduc)) +
  geom_point(aes(fill = Range), shape = 21, size = 3) + # Use fill for internal color based on 'Range'
  scale_fill_grey() +
  geom_smooth(method = "lm", se = FALSE, color = "black") + # Ensure trendline is black
    facet_wrap(~Range, labeller = labeller(Range = c(Invasive = "Familiar", Native = "Unfamiliar"))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    legend.position = "none",
    legend.background = element_rect(colour = "black", fill = "white", linetype = "solid"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    strip.text.x = element_text(color = "black", size = 20, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  labs(y = "Inducibility", x = "Constitutive level (Protein)")
)
ProTradeFig

# POD correlations
# removed constitutive values less than 0 as unlikely
PodCorrOut <- subset(PodCorr, PODHr0 > 0)

head(PodCorrOut)


PodCnat <- subset(PodCorrOut, Range == "Native")
PodCinv <- subset(PodCorrOut, Range != "Native")

cor.test(PodCnat$Tr_Hour0, PodCnat$Tr_DegInduc_2, method = "kendall")
# Kendall's rank correlation tau
#
# data:  PodCnat$Tr_Hour0 and PodCnat$Tr_DegInduc_2
# z = -3.7968, p-value = 0.0001466
# alternative hypothesis: true tau is not equal to 0 # negative correlation between starting and degree of inducibility for POD within native
# sample estimates:
#       tau
# -0.178762

cor.test(PodCinv$Tr_Hour0, PodCinv$Tr_DegInduc_2, method = "kendall")
# Kendall's rank correlation tau
#
# data:  PodCinv$Tr_Hour0 and PodCinv$Tr_DegInduc_2
# T = 105, p-value = 1.353e-07
# alternative hypothesis: true tau is not equal to 0 #also within invasive
# sample estimates:
#        tau
# -0.6022727

cor.test(PodCorrOut$Tr_Hour0, PodCorrOut$Tr_DegInduc_2, method = "kendall")
# Kendall's rank correlation tau
#
# data:  PodCorrOut$Tr_Hour0 and PodCorrOut$Tr_DegInduc_2
# z = -5.0188, p-value = 5.198e-07
# alternative hypothesis: true tau is not equal to 0 #overall
# sample estimates:
#       tau
# -0.218964

Pod.list <- list(PodCinv, PodCnat)
cocor(~ PODHr0 + DegrInduc_2 | PODHr0 + DegrInduc_2, Pod.list)
# Results of a comparison of two correlations based on independent groups
#
# Comparison between r1.jk (PODHr0, DegrInduc_2) = -0.7486 and r2.hm (PODHr0, DegrInduc_2) = -0.6159
# Difference: r1.jk - r2.hm = -0.1327
# Data: Pod.list: j = PODHr0, k = DegrInduc_2, h = PODHr0, m = DegrInduc_2
# Group sizes: n1 = 33, n2 = 204
# Null hypothesis: r1.jk is equal to r2.hm
# Alternative hypothesis: r1.jk is not equal to r2.hm (two-sided)
# Alpha: 0.05
#
# fisher1925: Fisher's z (1925)
#   z = -1.2847, p-value = 0.1989 # not seemingly different in correlation comparison
#   Null hypothesis retained
#
# zou2007: Zou's (2007) confidence interval
# 95% confidence interval for r1.jk - r2.hm: -0.2847 0.0851
# Null hypothesis retained (Interval includes 0)

(PodTradeFig <- ggplot(PodCorrOut, aes(x = PODHr0, y = DegrInduc_2)) +
  geom_point(shape = 21, size = 2, fill = "grey") +
  facet_wrap(~Range) +
  theme_pubr() +
  theme(strip.background = element_blank(), axis.text = element_text(size = 18), strip.text = element_blank(), axis.title = element_text(size = 20), legend.text = element_text(size = 17), legend.title = element_text(size = 18)) +
  labs(y = "Inducibility", x = "Constitutive level (POD)") +
  geom_smooth(method = "lm", se = F))

# updated POD corrplot
(PodTradeFig <- ggplot(PodCorrOut, aes(x = PODHr0, y = DegrInduc_2)) +
  geom_point(aes(fill = Range), shape = 21, size = 3) + # Use fill for internal color based on 'Range'
    scale_fill_grey() +
  geom_smooth(method = "lm", se = FALSE, color = "black") + # Ensure trendline is black
    facet_wrap(~Range, labeller = labeller(Range = c(Invasive = "Familiar", Native = "Unfamiliar"))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    strip.text.x = element_text(color = "black", size = 20, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  labs(y = "Inducibility", x = "Constitutive level (POD)"))

PodTradeFig

# PPO correlations
PpoCorrOut <- subset(PpoCorr, PPOHr0 > 0)

PpoCnat <- subset(PpoCorrOut, Range == "Native")
PpoCinv <- subset(PpoCorrOut, Range != "Native")

cor.test(PpoCnat$Tr_Hour0, PpoCnat$Tr_DegInduc_2, method = "kendall")
# Kendall's rank correlation tau
#
# data:  PpoCnat$Tr_Hour0 and PpoCnat$Tr_DegInduc_2
# z = -10.952, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0 #correlation between time 0 and the inducibility protein (using degree of inducibility is because you need to take into account the starting point to make claims about the final point of whatever you're testing ) This is native
# sample estimates:
#        tau
# -0.5131691

cor.test(PpoCinv$Tr_Hour0, PpoCinv$Tr_DegInduc_2, method = "kendall")
# Kendall's rank correlation tau
#
# data:  PpoCinv$Tr_Hour0 and PpoCinv$Tr_DegInduc_2
# z = -3.9055, p-value = 9.402e-05
# alternative hypothesis: true tau is not equal to 0 #same as above but invasive
# sample estimates:
#        tau
# -0.4781792

cor.test(PpoCorrOut$Tr_Hour0, PpoCorrOut$Tr_DegInduc_2, method = "kendall")
# Kendall's rank correlation tau
#
# data:  PpoCorrOut$Tr_Hour0 and PpoCorrOut$Tr_DegInduc_2
# z = -11.444, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0  #overall
# sample estimates:
#        tau
# -0.4972129

Ppo.list <- list(PpoCnat, PpoCinv)
cocor(~ PPOHr0 + DegrInduc_2 | PPOHr0 + DegrInduc_2, Ppo.list)
# Results of a comparison of two correlations based on independent groups
#
# Comparison between r1.jk (PPOHr0, DegrInduc_2) = -0.768 and r2.hm (PPOHr0, DegrInduc_2) = -0.6781
# Difference: r1.jk - r2.hm = -0.0899
# Data: Ppo.list: j = PPOHr0, k = DegrInduc_2, h = PPOHr0, m = DegrInduc_2
# Group sizes: n1 = 206, n2 = 33
# Null hypothesis: r1.jk is equal to r2.hm
# Alternative hypothesis: r1.jk is not equal to r2.hm (two-sided)
# Alpha: 0.05
#
# fisher1925: Fisher's z (1925)
#   z = -0.9702, p-value = 0.3319 #same as pod
#   Null hypothesis retained
#
# zou2007: Zou's (2007) confidence interval
# 95% confidence interval for r1.jk - r2.hm: -0.3369 0.0731
# Null hypothesis retained (Interval includes 0)

(PpoTradeFig <- ggplot(PpoCorrOut, aes(x = PPOHr0, y = DegrInduc_2)) +
  geom_point(shape = 21, size = 2, fill = "grey") +
  facet_wrap(~Range) +
  theme_pubr() +
  theme(strip.background = element_blank(), axis.text = element_text(size = 18), strip.text = element_blank(), axis.title = element_text(size = 20), legend.text = element_text(size = 17), legend.title = element_text(size = 18)) +
  labs(y = "Inducibility", x = "Constitutive level (PPO)") +
  geom_smooth(method = "lm", se = F))

# updated corr plot
(PpoTradeFig <- ggplot(PpoCorrOut, aes(x = PPOHr0, y = DegrInduc_2)) +
  geom_point(aes(fill = Range), shape = 21, size = 3) + # Use fill for internal color based on 'Range'
    scale_fill_grey() +
  geom_smooth(method = "lm", se = FALSE, color = "black") + # Ensure trendline is black
    facet_wrap(~Range, labeller = labeller(Range = c(Invasive = "Familiar", Native = "Unfamiliar"))) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 20),
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    strip.text.x = element_text(color = "black", size = 20, hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.title = element_text(size = 22, hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  ) +
  labs(y = "Inducibility", x = "Constitutive level (PPO)")
)
PpoTradeFig

# Plot all three together
cowplot::plot_grid(ProTradeFig, PodTradeFig, PpoTradeFig,
  labels = "AUTO", align = "v", nrow = 3, label_size = 20,
  rel_heights = c(1.2, 1.2, 1.2)
) # Spacing controlled by smaller values for empty plots
ggsave("Constit_Tradeoff16aug2024.pdf")
# save in E and E format. Save as a pdf, 600 dpi, min of 1800 px
ggsave("Figure5_onstit_Tradeoff16aug2024.jpg", width = 1800, units = px, height = 1800, dpi = 600)

# updated combo corr plot
corr_combo <- plot_grid(ProTradeFig, PodTradeFig, PpoTradeFig, labels = c("A", "B", "C"), align = "v", ncol = 1)
corr_combo

corr_combo2 <- (ProTradeFig / PodTradeFig / PpoTradeFig) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 18))

corr_combo2
ggsave("Constit_Tradeoff_panel_16aug2024.jpg", corr_combo2, width = 3800, height = 3100, units = "px", dpi = 300)
ggsave("Figure5_constit_Tradeoff16aug2024.pdf",corr_combo2, width = 5200, height = 5400, units = "px",  dpi = 600)

#---- generate correlations for combined dataframe export ----
# Adding a new column to each dataframe to indicate the variable type (Protein, POD, or PPO)
ProCorr$Type <- "Protein"
PodCorrOut$Type <- "POD"
PpoCorrOut$Type <- "PPO"

# Selecting and renaming the columns to have a consistent format across dataframes
ProCorr_selected <- ProCorr %>%
  select(Sample, Replicate, Range, Constitutive = ProteinHr0, Inducibility = NewDegrInduc, Type)

PodCorr_selected <- PodCorrOut %>%
  select(Sample, Replicate, Range, Constitutive = PODHr0, Inducibility = DegrInduc_2, Type)

PpoCorr_selected <- PpoCorrOut %>%
  select(Sample, Replicate, Range, Constitutive = PPOHr0, Inducibility = DegrInduc_2, Type)

# Combining all three dataframes into one
combined_data <- bind_rows(ProCorr_selected, PodCorr_selected, PpoCorr_selected)

# Viewing the combined dataframe
head(combined_data)

# Optionally, save this combined dataframe to a CSV file
write.csv(combined_data, "Constitutive_Inducibility_CombinedData.csv", row.names = FALSE)
