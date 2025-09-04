# AC_age_gr_dietsize.R
# This script reads in the Arctic cod size, Gill raker and diet item sizes 
# as measured during dissection under the stereomicroscope in the lab.
# Arctic cod size is converted to age using von Bertalanffy with 
# the relationship and parameters derived by Forster et al. 2021 (as per 
# Malizia et al. 2023). Gillraker data is used to calculate GR density, 
# and diet item sizes are used to estimate mean diet item size for each 
# Individual.

# Written by dr June 2024

## clearing all instances (data and variables)
rm(list = ls())

## loading needed libraries 
{
  library(dplyr)
  library(emmeans)
  library(ggplot2)
  library(lmerTest)
  library(visreg)
  library(stats)
  library(car)
  library(rstatix)
  library(ggpubr)
}

# Setting the working directory to where data and scripts are located
setwd("Path/to your file/")

# Reading in the diet data from the file. Looking for "AC_size_gr_dietsize.csv" file.
acdiet <- read.csv("your file here.csv", header=TRUE, stringsAsFactors = T, na.strings = "NA")

# Ordering the diet data by Fish_ID
acdiet <- acdiet[order(acdiet$Fish_ID),]

# Looking at the data
head(acdiet)

# Dataset is odd with individuals repeated to get size of items in stomachs. 
# Need to collapse data back into individuals as rows and variables as columns.

# Removing "Gut_Content" from named columns 
acdiet2 <- acdiet[, !names(acdiet) %in% "Gut_Content"]

# use dplyr to pipe data into group_by cmd, and summarise data
# by individual Fish_ID - but not site
newdata <- acdiet2 %>%
  group_by(across(where(is.factor))) %>%
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = T)))

# look at newdat to verify it looks good
newdata

# Write data out to summarise for later use in other analyses
#write.table(newdata,file = "compiled_dat_csv", 
#            append = F, quote = F)

# Test relationship between length and GR density
# both variables need log transformation to fit linear assumptions 
gr_l1<-lm(log(Raker_Density_T) ~ log(Body_Length), data = na.omit(newdata))

# summarise relationship
summary(gr_l1)

# Plot out diagnostics demonstrating fit
plot(gr_l1)

# plot out relationship for SI showing the data
ggplot(newdata, aes(x = log(Body_Length), y = log(Raker_Density_T)))+
  geom_point(size = 6, col = "skyblue3", alpha = 0.6, stroke = 1) +
  labs(x="log(Body length(mm))", y="log(Gill raker density)") +
  geom_smooth(method = "lm", se = T, color = "firebrick4", linetype = "dashed") +
  theme(aspect.ratio = 0.8) +
  theme_classic() +
  scale_y_continuous(limits = c(-0.1,1.5), breaks = seq(-0.1, 1.5, by = 0.2)) +
  scale_x_continuous(limits = c(4, 6), breaks = seq(4, 6, by = 0.2)) +
  theme(axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 24),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"))


################################## AGE CLASS #################################
# Adding age class of ind to newdata using total lengths in von Bertalanffy 
# formula with parameters derived from Forster et al. 2020

# Arctic cod VBGF parameters
Linf_ac <- 290.7462
K_ac <- 0.2175
t0_ac <- -0.8603

# Reversed function for age estimates for given lengths in mm.
VBage_est <- ((log(1 - (as.numeric(newdata$Body_Length) / Linf_ac))) / -K_ac) + t0_ac

# Round the VBage_est to the nearest integer to create year classes.
VBage <- ceiling(VBage_est)

# Add VBage to newdata
newdata$VBage <- as.factor(VBage)

# Check number of individuals in each age class
table(newdata$VBage)

# Filter age classes with too few individuals
bsac <- newdata[which(newdata$VBage != "6" & newdata$VBage != "5" & newdata$VBage != "14"),]

# Filter individuals from BB, not BS
bsac <- bsac[which(bsac$site != "BB"),]

# Check number of individuals in each age class and remove the 
# empty levels of the factor
table(bsac$VBage)
bsac$VBage <- droplevels(bsac$VBage) 

## Run basic linear model to assess relationship between length and VBage:
agelen <- lm(Body_Length ~ VBage, data = bsac)

# Use results of lm above to generate the means, SEs, and 95% CIs 
# for means of each age class
agelenmat <- as.data.frame(emmeans(agelen, "VBage",type = "response"))

# print the age-len df to see and report values in a table
agelenmat

# Do the same for the gillraker density - which we use later
grage <- lm(Raker_Density_T ~ VBage, data = bsac)
gragemat <- as.data.frame(emmeans(grage, "VBage",type = "response"))
gragemat

# Plot length ~ age relationship from the VBGF for Arctic cod using 
# Beaufort Sea fish only 
ggplot(agelenmat, aes(x = emmean, y = as.factor(VBage))) +
  geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL), width = 0.3, linewidth = 1) +
  geom_point(shape = 19, size = 6) +
  geom_jitter(bsac, mapping = aes(x = Body_Length, y = VBage, color = VBage), 
              position = position_jitter(0.3), size = 6, stroke = 1, alpha = 0.6) +
  labs(x = "Length (mm)", y = "VBage estimate", color = "VBage") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")) +
  scale_x_continuous(limits = c(50, 200), breaks = seq(50, 200, by = 25)) +
  theme(aspect.ratio = 0.80) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 24),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"),
        legend.position = "none")


################# Testing Gill raker density differences by Ageclass #########   

# Levene's test assesses Gill raker density departure from normality 
levgrd <- leveneTest(bsac$Raker_Density_T ~ as.factor(bsac$VBage))
levgrd

# Results > 0.05, gill raker density not considered deviating from 
# normal distribution

# Ran linear model between age vs. Raker_Density_T above
# called it grage. Can re-summarise here
summary(grage)

# Shapiro-Wilk's test for normality of residuals from linear model. 
shapiro_test(residuals(grage))

# Results sign (< 0.05), so likely deviating from normal distribution

# Assess deviations in the qqplot.
ggqqplot(residuals(grage))

# Boundary exceeded in many cases, so data deviate from normal

# Therefore, use glm with Gamma family distribution to account for 
# skew in data.
gr_age <- glm(Raker_Density_T ~ VBage, data = bsac, family = Gamma(link = "identity"))

# Can try others to make sure
gr_age2 <- glm(Raker_Density_T ~ VBage, data = bsac, family = gaussian(link = "identity"))

# Compare the 3 models with AIC
AIC(grage,gr_age, gr_age2)

# Summarise relationship and see the estimates
summary(gr_age)

# Use anova to test hypothesis. Here, error dist means little as all 
# give similar results (can switch "Chisq" to see). Go with Chisq with 
# deviance estimating goodness of fit G.
anova(gr_age, test = "Chisq")

# Sign difference over all ages - now use emmeans to do pairwise assessments
# using estimated marginal means, which are the means for each level of the
# factor, adjusted for other predictors in the model.
emmeans_gr <- emmeans(gr_age, ~ VBage)

# Use the pairs to do pairwise assessments
phpair_gr <- pairs(emmeans_gr)

# display results adjusted for multiple comparisons
phpair_gr

grdd_res <- as.data.frame(emmeans_gr)

## Generate stripcharts with 95%CI boxplots of the gillraker counts by age class
ggplot(grdd_res, aes(x = VBage, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_point(shape = 19, size = 6) +
  geom_jitter(bsac, mapping = aes(x = VBage, y = Raker_Density_T, 
                                    color = VBage), position = position_jitter(0.3), size = 6, stroke = 1, alpha = 0.6) +
  labs(x = "Age class", y = "Gill raker density (#/mm)", color = "maj_class") +
  scale_color_manual(name = "Age class", values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")) +
  scale_y_continuous(limits = c(1, 4), breaks = seq(1, 4, by = 0.5)) +
  theme(aspect.ratio = 0.70) +
  annotate("text", x = 1, y = 4, label = paste("X"), color = "black", size = 10, family = "Courier") +
  annotate("text", x = 2, y = 4, label = paste("Y"), color = "black", size = 10, family = "Courier") +
  annotate("text", x = c(3,4), y = c(4,4), label = paste("Z"), color = "black", size = 10, family = "Courier") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 24),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"),
        legend.position = "none")


############################ Diet item size ##################################

# Levene's test assesses diet item size (dis) departure from normality 
levdis <- leveneTest(bsac$Content_Size ~ as.factor(bsac$VBage))
levdis

# Run a linear model to test for significant differences between size of 
#diet items among VB age classes.
disdif <- lm(Content_Size ~ VBage, data = bsac)

# Summarising the linear model to get the coefficients
summary(disdif)

# Running a Shapiro-Wilk's test to assess the normality of residuals from 
# the linear model. 
# Results that are not sign. (< 0.05), are not considered to 
# be deviating from a normal distribution.
shapiro_test(residuals(disdif))

# Results are significant, so our data is deviated from normal. 
# Meaning parametric tests, unless transformed, would not be appropriate.

# Quickly assess major deviations in the qqplot.
ggqqplot(residuals(disdif))

# Trying the Gamma glm, but first have to verify that the data do not have negative 
# numbers or 0s, as these are not translatable in Gamma distributions.
bsac <- bsac %>%
  filter(Content_Size > 0)

# Trying the glm with Gamma error distribution
disdif1 <- glm(Content_Size ~ VBage, na.action = na.exclude, 
              family = Gamma(link = "identity"), data = bsac)

# Can also try other Gamma link functions
disdif2 <- glm(Content_Size ~ VBage, na.action = na.exclude, 
              family = Gamma(link = "inverse"), data = bsac)

# Can also try other error distributions link functions
disdif3 <- glm(Content_Size ~ VBage, na.action = na.exclude, 
               family = gaussian(link = "inverse"), data = bsac)


# Compare the 3 models with AIC (there is a warning here)
# but both disdif1 and 2 are identical
AIC(disdif,disdif1,disdif2, disdif3)

# Summarise relationships and see the estimates
summary(disdif1)
summary(disdif2)

# Use anova to test hypothesis. Here, error dist. different from normal, so 
# use "LRT" to test significance of deviance estimating goodness of fit G.
# Since the disdif 1 has fewer scoring iterations and lower dispersion
# go with that one.
anova(disdif1, test = "LRT")

# Sign difference overall ages - now use emmeans to do pairwise assessments
# using estimated marginal means, which are the means for each level of the
# factor, adjusted for other predictors in the model.
emmeans_dis <- emmeans(disdif1, ~ VBage)

# Use the pairs to do pairwise assessments
phpair_dis <- pairs(emmeans_dis)

# display results adjusted for multiple comparisons
phpair_dis

## Generate a dataframe of the results from an emmeans call that estimates the 
## coefficients and their 95% CIs. 
disdif_res <- as.data.frame(emmeans(disdif1, "VBage", type = "response"))

## The tabulated transformed coefficients and 95% CI can be visualised: 
disdif_res

## Use ggplot2 to generate stripcharts with 95%CI boxplots of the Mean Diet Size counts by age class
ggplot(disdif_res, aes(x = VBage, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_point(shape = 19, size = 6) +
  geom_jitter(bsac, mapping = aes(x = VBage, y = Content_Size, 
                                  color = VBage), position = position_jitter(0.3), size = 6, stroke = 1, alpha = 0.6) +
  labs(x = "Age class", y = "Mean diet item size (in mm)", color = "maj_class") +
  scale_color_manual(name = "Age class", values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 2.5)) +
  theme(aspect.ratio = 0.70) +
  annotate("text", x = 1, y = 15, label = paste("X"), color = "black", size = 12, family = "Courier") +
  annotate("text", x = c(2,3,4), y = c(15,15,15), label = paste("Y"), color = "black", size = 12, family = "Courier") +
    theme_classic() +
  theme(axis.text.x = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"),
        legend.position = "none")

###### END OF THESE ANALYSES

# Write data out to summarise for later use in other analyses
write.table(newdata,file = "ACcompileddat.csv", 
            append = F, quote = F)

