# AC_SIA.R 
# Script that reads in stable isotope data for determining from individual Arctic
# cod collected from the Beaufort Sea and assess them for differences among the 
# estimated age classes (determined from their length and calculated using the
# von Bertalanffy growth equation). 

# The data are plotted, and lm used to test differences among age classes overall 
# then tested in pairwise assessments. Plotted data have 95% ellipses plotted, 
# and these are used to calculate the ellipse eccentricity. Eccentricity gives an
# estimate of diet specificity or constraints relative to the isotopes used 
# to describe it.

# written by dr June 2024

# Clearing all instances
rm(list=ls())

#Loading needed libraries
{
  library(tidyverse)
  library(ggplot2)
  library(plyr)
  library(mvtnorm) # References rmvnorm()
  library(emmeans)
  library(visreg)
  library(car)
  library(dplyr)
  library(SIBER)
  library(rjags)
  library(lsr)
  library(conover.test)
  library(FSA)
  library(parallel)
  library(gt)
  library(lmerTest)
  library(lme4)
}

# Setting the core allocation to run parallel computations
options(mc.cores = parallel::detectCores())
cores <- detectCores(logical = T)

# Looking for the metadata which contains SI data "AC_metadata_062024.csv"
setwd("/path/to your data/here/")

# load the data. Looking for "AC_metadata_062024.csv" here.
sidat <- read.csv("your file here.csv", header = T, stringsAsFactors = T)
attach(sidat)

#sidat$VBage <- as.factor(VBage)

# filter sidat to exclude individuals with NA, those from BB, and 
# those younger than age class 4
fil_sidat <- sidat %>%
  filter(!is.na(calDN1514)) %>% # Remove rows with NA in the calDN1514
  filter(site == "BS") %>%
  filter(VBage <= 4)

# review the data
head(fil_sidat)

# set working colours similar to those used in other scripts.
mycol <- c("#E69F00", "#56B4E9", "#009E73","#F0E442")

################################ C and N #############################
# Get the ranges of the C and N isotope ratio values over all individual
# for plotting
range(fil_sidat$calDC1312)
range(fil_sidat$calDN1514)

# plot the data with the 95% ellipses with ggplot2
ggplot(fil_sidat, aes(x = calDC1312, y = calDN1514, color = as.factor(VBage), fill = as.factor(VBage))) + 
  geom_point(size = 6.5, stroke = 1, alpha = 0.7) +
  stat_ellipse(aes(fill = as.factor(VBage)), geom = "polygon", color = NA, alpha = 0.3, level = 0.95) +
  stat_ellipse(aes(color = as.factor(VBage)), geom = "path", linewidth = 1, level = 0.95) +
  scale_fill_manual(values = mycol) +
  scale_color_manual(values = mycol) +
  labs(x = expression({delta}^13*C~'ppm'), y = expression({delta}^15*N~'ppm')) +
  #geom_hline(yintercept = 15.3, linewidth = 1, color = "firebrick4", linetype = "dotted") +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  scale_x_continuous(limits = c(-25, -20), breaks = seq(-25, -20, by = 1)) +
  scale_y_continuous(limits = c(13, 18), breaks = seq(13, 18, by = 1)) +
  annotate("text", x = c(-23.25,-23,-22.6,-22.2), y = c(18,18,18,18), label = paste(c("X","X","XY","Y")), color = c("#E69F00", "#56B4E9", "#009E73","yellow3"), size = 12, family = "Courier") +
  annotate("text", x = c(-20,-20,-20,-20), y = c(14.3,15.3,15.8,16.2), label = paste(c("X","Y","Z","Z")), color = c("#E69F00", "#56B4E9", "#009E73","yellow3"), size = 12, family = "Courier") +
  theme(axis.text.x = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth = 1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3, "cm")) +
  theme(legend.text = element_text(size = 18)) +
  theme(legend.title = element_blank())

# Calculating ellipse areas and eccentricities CN

# setting an empty matrix to store the ellipse area and eccentricity 
# doing it for C and N first
CN_elin <- matrix(data = NA, 4, 5)

# Using split function to split data by VBage
agedata <- split(fil_sidat, fil_sidat$VBage)

# Looping through all age classes to calculate ellipse area
# eccentricity. Ec ~ 1 means narrow, Ex ~ 0 means broad

for (a in 1:4) { # looping through 4 age classes
  # set temporary df to play with data
  tmp <- as.data.frame(agedata[a], col.names = NULL) 
  
# calculate the means of the SI vectors for the age class
  mu <- cbind(mean(tmp$calDC1312), mean(tmp$calDN1514))
  
# calculate the covariance matrix for the two isotope ratios
  covm <- cov(cbind(tmp$calDC1312,tmp$calDN1514))
  
# get the eigen values
  evals<-eigen(covm)$values
  
# and the eigen vectors
  evecs<-eigen(covm)$vectors
  
# get the values for the 95% CI of the data
  c<-sqrt(qchisq(0.95, 2))
  
# use the evals to calculate area
  elar<-c*(sqrt(evals[1])*sqrt(evals[2]))*pi
  
# and the ecc.
  elx<-sqrt(evals[1]^2-evals[2]^2)/evals[1]
  
# make a new matrix to store the bootstraps
  elit <- matrix(data = NA, 1000, 2)
  
# set internal loop to calculate bootstrap standard errors (bse)
  for (b in 1:1000) { # doing 1000 iterations
    
# randomly sample stable isotope ratios with replacement 
    xboot1 <- sample(tmp$calDC1312, replace = T)
    xboot2 <- sample(tmp$calDN1514, replace = T)

# combine and calculate cov of sample
    tcov <- cov(cbind(xboot1, xboot2))
    
# estimate ellipse axes of bootstrap sample
    teval <- eigen(tcov)$values
    
# calculate sample area and ecc.
    belar <- c*(sqrt(teval[1])*sqrt(teval[2]))*pi
    belx <- sqrt(teval[1]^2-teval[2]^2)/teval[1]
    
# store bootstrap result in elit matrix 
    elit[b,1] <- belar
    elit[b,2] <- belx
  }  

# calulate the sd (which is the se of bootstrap) to get bse
bse <- cbind(sd(elit[,1]),sd(elit[,2]))
# bind the bse values by column along with age and calculated elar and elx
# to make a df
CN_elin[a,]<-cbind(a,(elar),bse[1],elx,bse[2])

}

# Make the CN_elin matrix into df 
CN_elin<-as.data.frame(CN_elin)

# rename the columns to better names
CN_elin <- CN_elin %>%
  dplyr::rename(agecl = V1, el_area = V2, a_bse = V3, el_ecc = V4, ecc_bse = V5)

# display CN_elin to show results
CN_elin

# remove the looped items to be able to reuse without issues
rm(tmp,mu,covm,evals, evecs, elar, elx, elit, 
           xboot1, xboot2, tcov, teval, belar, belx)

################################ S and N #############################
# Get the ranges of the S and N isotope ratio values over all individual
# for plotting
range(fil_sidat$calDS3432)
range(fil_sidat$calDN1514)

# plot the data with the 95% ellipses with ggplot2
ggplot(fil_sidat, aes(x = calDS3432, y = calDN1514, colour=as.factor(VBage))) + 
  geom_point(size = 6.5, stroke = 1, alpha = 0.7) + 
  stat_ellipse(aes(fill = as.factor(VBage)), geom = "polygon", color = NA, alpha = 0.3, level = 0.95) +
  stat_ellipse(aes(color = as.factor(VBage)), geom = "path", linewidth = 1, level = 0.95) +
  scale_fill_manual(values = mycol) +
  labs(x = expression({delta}^34*S~'ppm'), y = expression({delta}^15*N~'ppm')) +
  scale_color_manual(values = mycol) +
  #geom_hline(yintercept = 15.3, linewidth = 1, color="firebrick4", linetype="dotted") +
  scale_x_continuous(limits = c(16.5, 21), breaks = seq(16.5, 21, by = 1)) +
  scale_y_continuous(limits = c(13, 18), breaks = seq(13, 18, by = 1)) +
  annotate("text", x = c(20.8,20.8,20.8,20.8), y = c(14.3,15.3,15.8,16.2), label = paste(c("X","Y","Z","Z")), color = c("#E69F00", "#56B4E9", "#009E73","yellow3"), size = 12, family = "Courier") +
  theme(aspect.ratio = 0.70) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm")) +
  theme(legend.text = element_text(size = 18))+
  theme(legend.title = element_blank())

# Calculating ellipse areas and eccentricities SN

# setting an empty matrix to store the ellipse area and eccentricity 
# doing it for S and N first
SN_elin <- matrix(data = NA, 4, 5)

# Looping through all age classes to calculate ellipse area
# eccentricity. Ec ~ 1 means narrow, Ex ~ 0 means broad

for (a in 1:4) { # looping through 4 age classes
  
# set temporary df to play with data
  tmp <- as.data.frame(agedata[a], col.names = NULL) 
  
# calculate the means of the SI vectors for the age class
  mu <- cbind(mean(tmp$calDS3432), mean(tmp$calDN1514))
  
# calculate the covariance matrix for the two isotope ratios
  covm <- cov(cbind(tmp$calDS3432,tmp$calDN1514))
  
# get the eigen values
  evals<-eigen(covm)$values
  
# and the eigen vectors
  evecs<-eigen(covm)$vectors
  
# get the values for the 95% CI of the data
  c<-sqrt(qchisq(0.95, 2))
  
# use the evals to calculate area
  elar<-c*(sqrt(evals[1])*sqrt(evals[2]))*pi
  
# and the ecc.
  elx<-sqrt(evals[1]^2-evals[2]^2)/evals[1]
  
# make a new matrix to store the bootstraps
  elit <- matrix(data = NA, 1000, 2)
  
# set internal loop to calculate bootstrap standard errors (bse)
  
for (b in 1:1000) { # doing 1000 iterations
    
# randomly sample stable isotope ratios with replacement 
    xboot1 <- sample(tmp$calDS3432, replace = T)
    xboot2 <- sample(tmp$calDN1514, replace = T)
    
# combine and calculate cov of sample
    tcov <- cov(cbind(xboot1, xboot2))
    
# estimate ellipse axes of bootstrap sample
    teval <- eigen(tcov)$values
    
# calculate sample area and ecc.
    belar <- c*(sqrt(teval[1])*sqrt(teval[2]))*pi
    belx <- sqrt(teval[1]^2-teval[2]^2)/teval[1]
    
# store bootstrap result in elit matrix 
    elit[b,1] <- belar
    elit[b,2] <- belx
  }  
  
# calulate the sd (which is the se of bootstrap) to get bse
  bse <- cbind(sd(elit[,1]),sd(elit[,2]))
  
# bind the bse values by column along with age and calculated elar and elx
# to make a df
  SN_elin[a,]<-cbind(a,(elar),bse[1],elx,bse[2])
}

# Make the SN_elin matrix into df 
SN_elin<-as.data.frame(SN_elin)

# rename the columns to better names
SN_elin <- SN_elin %>%
  dplyr::rename(agecl = V1, el_area = V2, a_bse = V3, el_ecc = V4, ecc_bse = V5)

# display SN_elin to show results
SN_elin

# remove the looped items to be able to reuse without issues
rm(tmp,mu,covm,evals, evecs, elar, elx, elit, 
   xboot1, xboot2, tcov, teval, belar, belx)

################################ C and S #############################
# Get the ranges of the C and S isotope ratio values overall individual
# for plotting
range(fil_sidat$calDC1312)
range(fil_sidat$calDS3432)

# plot the data with the 95% ellipses with ggplot2
ggplot(fil_sidat, aes(x = calDC1312, y = calDS3432, colour=as.factor(VBage))) + 
  geom_point(size = 6.5, stroke = 1, alpha = 0.7) +
  stat_ellipse(aes(fill = as.factor(VBage)), geom = "polygon", color = NA, alpha = 0.3, level = 0.95) +
  stat_ellipse(aes(color = as.factor(VBage)), geom = "path", linewidth = 1, level = 0.95) +
  scale_fill_manual(values = mycol) +
  labs(x = expression({delta}^13*C~'ppm'), y = expression({delta}^34*S~'ppm')) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")) +
  scale_x_continuous(limits = c(-25, -20), breaks = seq(-25, -20, by = 1)) +
  scale_y_continuous(limits = c(16.5, 20.7), breaks = seq(16.5, 20.5, by = 1)) +
  annotate("text", x = c(-23.25,-23,-22.6,-22.1), y = c(20.6,20.6,20.6,20.6), label = paste(c("X","X","XY","Y")), color = c("#E69F00", "#56B4E9", "#009E73","yellow3"), size = 12, family = "Courier") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 26),
        axis.title.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"))+
  theme(legend.text = element_text(size = 18))+
  theme(legend.title = element_blank())

# Calculating ellipse areas and eccentricities CS

# setting an empty matrix to store the ellipse area and eccentricity 
# doing it for C and S first
CS_elin <- matrix(data = NA, 4, 5)

# Looping through all age classes to calculate ellipse area
# eccentricity. Ec ~ 1 means narrow, Ex ~ 0 means broad

for (a in 1:4) { # looping through 4 age classes
  
# set temporary df to play with data
  tmp <- as.data.frame(agedata[a], col.names = NULL) 
  
# calculate the means of the SI vectors for the age class
  mu <- cbind(mean(tmp$calDC1312), mean(tmp$calDS3432))
  
# calculate the covariance matrix for the two isotope ratios
  covm <- cov(cbind(tmp$calDC1312, tmp$calDS3432))
  
# get the eigen values
  evals<-eigen(covm)$values
  
# and the eigen vectors
  evecs<-eigen(covm)$vectors
  
# get the values for the 95% CI of the data
  c<-sqrt(qchisq(0.95, 2))
  
# use the evals to calculate area
  elar<-c*(sqrt(evals[1])*sqrt(evals[2]))*pi
  
# and the ecc.
  elx<-sqrt(evals[1]^2-evals[2]^2)/evals[1]
  
# make a new matrix to store the bootstraps
  elit <- matrix(data = NA, 1000, 2)
  
# set internal loop to calculate bootstrap standard errors (bse)
  for (b in 1:1000) { # doing 1000 iterations
    
# randomly sample stable isotope ratios with replacement 
    xboot1 <- sample(tmp$calDC1312, replace = T)
    xboot2 <- sample(tmp$calDS3432, replace = T)
    
# combine and calculate cov of sample
    tcov <- cov(cbind(xboot1, xboot2))
    
# estimate ellipse axes of bootstrap sample
    teval <- eigen(tcov)$values
    
# calculate sample area and ecc.
    belar <- c*(sqrt(teval[1])*sqrt(teval[2]))*pi
    belx <- sqrt(teval[1]^2-teval[2]^2)/teval[1]
    
# store bootstrap result in elit matrix 
    elit[b,1] <- belar
    elit[b,2] <- belx
  }  
  # calulate the sd (which is the se of bootstrap) to get bse
  bse <- cbind(sd(elit[,1]),sd(elit[,2]))
  
# bind the bse values by column along with age and calculated elar and elx
# to make a df
  CS_elin[a,]<-cbind(a,(elar),bse[1],elx,bse[2])
}

# Make the CS_elin matrix into df 
CS_elin<-as.data.frame(CS_elin)

# rename the columns to better names
CS_elin <- CS_elin %>%
  dplyr::rename(agecl = V1, el_area = V2, a_bse = V3, el_ecc = V4, ecc_bse = V5)

# display SN_elin to show results
CS_elin

# remove the looped items to be able to reuse without issues
rm(tmp,mu,covm,evals, evecs, elar, elx, elit, 
   xboot1, xboot2, tcov, teval, belar, belx)

############################# Testing N signal ############################
# Convert VBage into a factor variable
fil_sidat$VBage <- as.factor(fil_sidat$VBage)

# order the loc west to east
fil_sidat$loc <- factor(fil_sidat$loc, levels = c("MSS","MSH","FKB","BKI",
                                                  "DLB", "MTI", "DUS", "ULU"))

# tabulate the N data to see how it's distributed by location and age.
tabdat <- as.data.frame.matrix(table(fil_sidat$loc, fil_sidat$VBage))

# Convert row names (locations) to a stub column
tabdat <- rownames_to_column(tabdat, var = "Location")

gt(tabdat, rownames_to_stub = F, auto_align = T) %>%
tab_spanner(label = "Age", columns = -Location) %>%
  tab_header(title = "Table X. Distribution of Arctic cod ages by Location") %>%
  cols_align(align = "center",columns = everything())

# Check if data is normally distributed
# Use a histogram to see if there is a normal bell-shaped curve with nitrogen data
hist(fil_sidat$calDN1514, main = "Histogram of Nitrogen isotopes")

# Run a linear model and test residuals for normality using 
# Shapiro-Wilk's test
lmN <- lm(calDN1514 ~ VBage, data = fil_sidat)

# Shapiro tests normality and distribution of lm residuals
shapiro.test(lmN$residuals)

# The data do not deviate from normality and so can use normal 
# based stats to assess differences.

# Want to assess whether there is a difference among age classes, 
# but also need to consider location to address reviewer concern. 
# We run a mixed-effect model that looks at locations as a random factor.
lmN2 <- lmer(calDN1514 ~ VBage + (1|loc), data = fil_sidat)

# Summarise the model to estimate parameters (do not trust t or p-values)
summary(lmN2)

# Looks like location loc (has little overall importance - 
# low amounts of var = 0.025 relative to Residuals 
# var = 0.301.

# Run an anova to test hypothesis of no differences in age
anova(lmN2, type = 1)

# results show important differences among VBage.

# Use emmeans to get the coefficients for each group and 
# to run the pairwise analyses (t-tests). Used the lmer.df 
# = satterthwaite 
emmeans(lmN2, pairwise ~ VBage, lmer.df = "satterthwaite")

# Can use visreg to plot the distributions of just the N ratios
# with 95% CI around means 
visreg(lmN2, line.par = list(col = "firebrick4", lty = 2),
       points = list(size = 5, col = "skyblue3", alpha = .6, stroke = 1),
       scale = "linear", rug = F, gg = T) +
  #aes(colour = factor(VBage)) +
  labs(x = "Age", y = expression({delta}^15*N~'ppm')) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442"))+
  scale_y_continuous(limits = c(13, 18), breaks = seq(13, 18, by = 0.5)) +
  theme(aspect.ratio = 0.8) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm")) +
  theme(legend.text = element_text(size = 17))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom", legend.justification = "center")

# To test the overall impact of N differences by location can 
# isolate set of individuals of same age class from as many locations
# as possible and test them for structure.

# From tabulated data - it seems we could isolate age 1 fish, at
# MSS, FKB, DLB, MTI, and DUS. 
Nta1 <- fil_sidat %>%
  filter(fil_sidat$VBage == 1) %>%
  filter(!loc %in% c("MSH", "BKI", "ULU"))

# eliminates the unwanted sampling location as levels
droplevels(Nta1$loc)

# Run basic lm to test residuals
lmNta1 <- lm(calDN1514 ~ loc, data = Nta1)

# Shapiro tests normality and distribution of lm residuals
shapiro.test(lmNta1$residuals)

# Test is significant - meaning the residuals are not 
# normally distributed. So can't easily run lm here 
#, and an alternative might be better

# Run Kruskal-Wallis test 
ksNta1 <- kruskal.test(calDN1514 ~ loc, data = Nta1)
ksNta1

# Overall KW test shows a high Chisqr value.
# Run pairwise test to see which groups are different

# Can use Dunn's test is the most reliable in terms of proper correction of p-values.
dunnTest(calDN1514 ~ loc, data = Nta1, method = "by")

# Calculate the mean and 95% CI for Nta1$calDN1514 per location
mN <- by(Nta1$calDN1514, Nta1$loc, na.omit(mean))
sdN <- by(Nta1$calDN1514, Nta1$loc, na.omit(sd))
nN <- by(Nta1$calDN1514, Nta1$loc, na.omit(length))
cil <- mN - 1.96 * (sdN / sqrt(nN))
cih <- mN + 1.96 * (sdN / sqrt(nN))

# Put the summary stats into a single dataframe.
sumstat <- cbind.data.frame(loc = as.factor(levels(droplevels(Nta1$loc))), 
                            mN = na.omit(as.numeric(mN)), cil = na.omit(as.numeric(cil)), 
                            cih = na.omit(as.numeric(cih)))

# Then plot using ggplot with jitter and geom_pointrange
ggplot(Nta1, aes(x = loc, y = calDN1514, colour = loc)) +
  geom_jitter(size = 5, width = 0.15, stroke = 1, alpha = 0.7) +
  geom_pointrange(data = sumstat,
                  aes(x = loc, y = mN, ymin = cil, ymax = cih),
                  position = position_dodge(width = 0.3),
                  size = 1.4, colour = "black") +
  # 95% CI is not appropriate here because the data are not normal
  labs(x = "Location", y = expression({delta}^15*N~'ppm')) + 
  scale_y_continuous(limits = c(13, 15.5), breaks = seq(13, 15.5, by = 0.5)) +
  scale_color_manual(values = c("darkgray", "firebrick4", "skyblue4", "orange", "springgreen4")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 1.3),
    axis.ticks = element_line(linewidth = 1.3),
    axis.ticks.length = unit(.3, "cm"),
    legend.text = element_text(size = 17),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )

############################# Testing C signal ############################
# Should do the same for calDC1312 to see if data are normal 
hist(fil_sidat$calDC1312, main = "Histogram of Carbon isotopes")

# Run a linear model and test residuals for normality using 
# Shapiro-Wilk's test
lmC <- lm(calDC1312 ~ VBage, data = fil_sidat)

# Shapiro tests normality and distribution of lm residuals
shapiro.test(residuals(lmC))

# Want to assess whether there is a difference among age classes, 
# but as above, need to consider location to address reviewer concern. 
# We run a mixed-effect model that looks at locations as random factor.
lmC2 <- lmer(calDC1312 ~ VBage + (1|loc), data = fil_sidat)

# Summarise the linear model to estimate parameters, but do not trust 
# t or p-values
summary(lmC2)

# As above it looks like location loc (has little overall importance - 
# low amounts of var = 0.055 relative to Residuals 
# var = 0.281.

# Run an anova to test hypothesis of no differences in age
anova(lmC2, type = 1)

# A significant difference in carbon signal is seen as well.
# emmeans to run the pairwise differences
emmeans(lmC2, pairwise ~ VBage, lmer.df = "satterthwaite")
 
# Only groups 1 vs. 4 and 2 vloc# Only groups 1 vs. 4 and 2 vs. 4 are different. The rest 
# are not different

# Adjust the p-values for multiple comparisons using false 
# discovery rate
pairwise.t.test(fil_sidat$calDC1312, fil_sidat$VBage, p.adj = "BY")

visreg(lmC2, line.par = list(col = "firebrick4", lty = 2),
       points = list(size = 5, col = "skyblue3", alpha = .6, stroke = 1),
       scale = "linear", rug = F, gg = T) +
  labs(x = "Age", y = expression({delta}^13*C~'ppm')) +
  scale_y_continuous(limits = c(-25, -20), breaks = seq(-25, -20, by = 1)) +
  theme(aspect.ratio = 0.8) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 30),
        axis.text.x = element_text(size = 28),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 28))

# Can just use the data isolated previously for the N data 
# but analyse the calDC1312 

# Run basic lm to test residuals
lmCta1 <- lm(calDC1312 ~ loc, data = Nta1)

# Shapiro tests normality and distribution of lm residuals
shapiro.test(lmCta1$residuals)

# In this case, the data do not deviate from normality 
# so can go ahead and interpret lm. 1st, summarise
summary(lmCta1)


# Then test the Ho of no difference among locations using anova
anova(lmCta1)

# Use emmeans to get the coefficients for each group and 
# to run the pairwise analyses (t-tests). 
emmeans(lmCta1, pairwise ~ loc)

# Calculate the mean and 95% CI for Nta1$calDC1312 per location
mC <- by(Nta1$calDC1312, Nta1$loc, na.omit(mean))
sdC <- by(Nta1$calDC1312, Nta1$loc, na.omit(sd))
nC <- by(Nta1$calDC1312, Nta1$loc, na.omit(length))
cilc <- mC - 1.96 * (sdC / sqrt(nC))
cihc <- mC + 1.96 * (sdC / sqrt(nC))

# Put the summary stats into a single dataframe.
sumstatC <- cbind.data.frame(loc = as.factor(levels(droplevels(Nta1$loc))), 
                            mC = na.omit(as.numeric(mC)), cilc = na.omit(as.numeric(cilc)), 
                            cihc = na.omit(as.numeric(cihc)))

# Then plot using ggplot with jitter and geom_pointrange
ggplot(Nta1, aes(x = loc, y = calDC1312, colour = loc)) +
  geom_jitter(size = 5, width = 0.15, stroke = 1, alpha = 0.7) +
  geom_pointrange(data = sumstatC,
                  aes(x = loc, y = mC, ymin = cilc, ymax = cihc),
                  position = position_dodge(width = 0.3),
                  size = 1.4, colour = "black") +
  labs(x = "Location", y = expression({delta}^13*C~'ppm')) + 
  scale_y_continuous(limits = c(-25, -20), breaks = seq(-25, 20, by = 0.5)) +
  scale_color_manual(values = c("darkgray", "firebrick4", "skyblue4", "orange", "springgreen4")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 1.3),
    axis.ticks = element_line(linewidth = 1.3),
    axis.ticks.length = unit(.3, "cm"),
    legend.text = element_text(size = 17),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )

############################# Testing S signal ############################
# Use a histogram to see if there is normal data 
hist(fil_sidat$calDS3432, main = "Histogram of Sulfur isotopes")

# Run a linear model and test residuals for normality using 
# Shapiro-Wilk's test for both loc = location, and VBage = age.
lmS <- lm(calDS3432 ~ VBage, data = fil_sidat)

# Shapiro tests normality and distribution of lm residuals
shapiro.test(residuals(lmS))

# Data do not conform well with normality. Could use non-parametric 
# test but also could use the glm with family Gamma, and identity link,
# which works well for positive skew in data. As for diet size and GRs.

# However, the data also need to account for location as we have done 
# already in the other SI data. Here we can use a glm with a mixed 
# effect (random effect of location).
lmS2 <- glmer(calDS3432 ~ VBage + (1|loc), family = Gamma(link = "identity"), data = fil_sidat)

# Summarise the model to see the effects. 
summary(lmS2)

# Found Anova from car package that can be used with glmer 
# So, test Ho of no difference among the VBages. 
Anova(lmS2, type = "III", singular.ok = T, test.statistic = "Chisq")
anova(lmS2, test.statistic = "Chisq")

# Accordingly, the amount of variance explained by both the 
# fixed and random factors is small. And while the intercept is not 0
# There is no difference among VBages, even when location is considered. 

# Can use emmeans to get the untransformed coefficients and their 95%
# the df is Inf, because they are approximated as very large from a 
# Laplace transformation, and this is typical for a glmer.
emmeans(lmS2, "VBage", type = "response")
emmeans(lmS2, pairwise ~ VBage, lmer.df = "satterthwaite")

visreg(lmS2, line.par = list(col = "firebrick4", lty = 2),
       points = list(size = 5, col = "skyblue3", alpha = .6, stroke = 1),
       scale = "linear", rug = F, gg = T) +
  labs(x = "Age", y = expression({delta}^34*S~'ppm')) +
  scale_y_continuous(limits = c(18.4, 19.4), breaks = seq(18.4, 19.4, by = 0.2)) +
  theme(aspect.ratio = 0.8) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 30),
        axis.text.x = element_text(size = 28),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 28))

# Could also assess differences by loc in calDS3432, by using the 
# reduced dataset of just age 1 fish.
lmSl <- lm(calDS3432 ~ loc, data = Nta1)

# Shapiro tests normality and distribution of lm residuals
shapiro.test(residuals(lmSl))

# No deviation from normality. So, can summarise and test.
summary(lmSl)

# Coefficients are similar. test Ho with anova.
anova(lmSl)

# Confirms no differences in S signal in the data whether 
# treated for VBage or location. 
sumstatS <- as.data.frame(emmeans(lmSl, "loc", type = "response"))
sumstatS

# Then plot using ggplot with jitter and geom_pointrange
ggplot(Nta1, aes(x = loc, y = calDS3432, colour = loc)) +
  geom_jitter(size = 5, width = 0.15, stroke = 1, alpha = 0.7) +
  geom_pointrange(data = sumstatS,
                  aes(x = loc, y = emmean, ymin = lower.CL, ymax = upper.CL),
                  position = position_dodge(width = 0.3),
                  size = 1.4, colour = "black") +
  labs(x = "Location", y = expression({delta}^34*S~'ppm')) + 
  scale_y_continuous(limits = c(17, 21), breaks = seq(17, 21, by = 0.5)) +
  scale_color_manual(values = c("darkgray", "firebrick4", "skyblue4", "orange", "springgreen4")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(linewidth = 1.3),
    axis.ticks = element_line(linewidth = 1.3),
    axis.ticks.length = unit(.3, "cm"),
    legend.text = element_text(size = 17),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.justification = "center"
  )

############################# SIBER ################################
# Below are analyses that can assess the amount of overlap of the ellipses 
# using the SIBER package in R. Some of these cmds take a while to run.

# first, create a siber object from the data 
siNC <- cbind.data.frame(iso1 = fil_sidat$calDC1312, iso2 = fil_sidat$calDN1514, 
                         group = fil_sidat$VBage, community = rep(1,nrow(fil_sidat)))


siNS <- cbind.data.frame(iso1 = fil_sidat$calDS3432, iso2 = fil_sidat$calDN1514, 
                         group = fil_sidat$VBage, community = rep(1,nrow(fil_sidat)))


siSC <- cbind.data.frame(iso1 = fil_sidat$calDC1312, iso2 = fil_sidat$calDS3432, 
                         group = fil_sidat$VBage, community = rep(1,nrow(fil_sidat)))



# Once the data are in the format create SIBER object
sidfNC <- createSiberObject(siNC)
sidfNS <- createSiberObject(siNS)
sidfSC <- createSiberObject(siSC)

# this remakes 95% ellipses for the data as above
group.ellipses.args  <- list(n = 1000, p.interval = 0.95, lty = 1, lwd = 2)

# Here, SIBER plots the data as well- helps to check if we get same patterns 
# as biplots (which are prettier).
par(mfrow=c(1,1))
plotSiberObject(sidfNC,
                ax.pad = 1, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L", las = 1,
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'ppm'),
                ylab = expression({delta}^15*N~'ppm'),
                cex = 0.5)

plotSiberObject(sidfNS,
                ax.pad = 1, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L", las = 1,
                iso.order = c(1,2),
                xlab = expression({delta}^34*S~'ppm'),
                ylab = expression({delta}^15*N~'ppm'),
                cex = 0.5)

plotSiberObject(sidfSC,
                ax.pad = 1, 
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = F, group.hull.args,
                bty = "L", las = 1,
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'ppm'),
                ylab = expression({delta}^34*S~'ppm'),
                cex = 0.5)

# Calculate the summary statistics for each group 
# Total Area - TA
# Standard Ellipse Area - SEA
# Standard Ellipse Area for small sample sizes - SEAc

# once plotted can look at size of each ellipse and determine 
# amount of overlap
MLNC <- groupMetricsML(sidfNC)
MLNS <- groupMetricsML(sidfNS)
MLSC <- groupMetricsML(sidfSC)

# print group ML data to screen
MLNC
MLNS
MLSC

# Identify the ellipses
a1 <-"1.1"
a2 <-"1.2"
a3 <-"1.3"
a4 <-"1.4"

####################### These next parts can take a while to run
# Parameter options for running gibbs sampler (jags) 
parms <- list()
parms$n.iter <- 2 * 10^6
parms$n.burin <- 1 * 10^4
parms$n.thin <- 100
parms$n.chains <- 5

# defining the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posteriorNC <- siberMVN(sidfNC, parms, priors)
ellipses.posteriorNS <- siberMVN(sidfNS, parms, priors)
ellipses.posteriorSC <- siberMVN(sidfSC, parms, priors)

# siberEllipses calculates the posterior probabilities of overlap
SEA.B_NC <- siberEllipses(ellipses.posteriorNC)
SEA.B_NS <- siberEllipses(ellipses.posteriorNS)
SEA.B_SC <- siberEllipses(ellipses.posteriorSC)

# Possibly add CI to ellipses/means?
cr.p <- c(0.95, 0.99)
SEA.B.NC.modes <- lapply(as.data.frame(SEA.B_NC),
  function(x,...){tmp <- hdrcde::hdr(x)$hdr},
  prob = cr.p)

# Setting empty matrix to store data  
BONC <- matrix(data = NA, 16, 3)
BONS <- matrix(data = NA, 16, 3)
BOSC <- matrix(data = NA, 16, 3)

# Initialising counter
d <- 0

# loop through the possible permutations to calculate ellipse overlap
for (i in 1:4) { 
  for (j in 1:4) {
    d <- d + 1 #increasing counter
    # Calculate Bayesian area and overlap of 
    tmp1 <- bayesianOverlap(i, j, ellipses.posteriorNC,draws = 10, p.interval = 0.95, n = 360)
    tmp2 <- bayesianOverlap(i, j, ellipses.posteriorNS,draws = 10, p.interval = 0.95, n = 360)
    tmp3 <- bayesianOverlap(i, j, ellipses.posteriorSC,draws = 10, p.interval = 0.95, n = 360)
    
    tmp1 <- as.data.frame(tmp1)
    tmp2 <- as.data.frame(tmp2)
    tmp3 <- as.data.frame(tmp3)
    
    mtmp1 <- apply(tmp1,2,mean)
    mtmp2 <- apply(tmp2,2,mean)
    mtmp3 <- apply(tmp3,2,mean)
    
    BONC[d,] <- mtmp1
    BONS[d,] <- mtmp2
    BOSC[d,] <- mtmp3
  }
}
# convert bayesian matrix to df
BONC <- as.data.frame(BONC)
BONS <- as.data.frame(BONS)
BOSC <- as.data.frame(BOSC)

# rename the columns to better names
BONC <- BONC %>%
  dplyr::rename(area1 = V1, area2 = V2, b_over = V3)

# rename the columns to better names
BONS <- BONS %>%
  dplyr::rename(area1 = V1, area2 = V2, b_over = V3)

# rename the columns to better names
BOSC <- BOSC %>%
  dplyr::rename(area1 = V1, area2 = V2, b_over = V3)

# Set sequences iding the comparisons
el_a <- rep(1:4, each = 4)
el_b <- rep(1:4, times = 4)

# Recombine the df with indices
BONC <- cbind(el_a, el_b, BONC)
BONS <- cbind(el_a, el_b, BONS)
BOSC <- cbind(el_a, el_b, BOSC)


BONC$bo_1 <- (BONC$b_over/BONC$area1)*100
BONC

BONS$bo_1 <- (BONS$b_over/BONS$area1)*100
BONS

BOSC$bo_1 <- (BOSC$b_over/BOSC$area1)*100
BOSC
# To access the overlap - print out the BOXX dfs 

#### END

