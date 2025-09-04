# AC_diet_item.R
# This script reads in Arctic cod diet data compiled during diet analyses.
# Data are listed as the size of 18 different possible diet categories recovered
# from the stomachs of each individual. Each item was measured to get diet sizes
# but were also counted. So entries are one or more single rows per individual.
# The data thus have to be compiled by individual and then assessed using the 
# points method as used in Genner et al. 2003 and in Roy et al. 2007.

# Points data are converted to a dissimilarity matrix, and then evaluated using
# Adonis2/PERMANOVA to assess overall differences as provided through the vegan
# and the pairwiseAdonis packages to test for pairwise comparisons among age classes.

# written by DR June 2024

# Clearing instances and resetting environment
rm(list=ls())

## Loading needed libraries
{
  library(devtools)
  library(ggplot2)
  library(tidyverse)
  library(parallel)
  library(dplyr)
  library(vegan)
  library(tidyr)
  library(pairwiseAdonis)
  library(ggbiplot)
  library(lessR)
}
# install the pairwiseAdonis package from GithHub (may only be needed once).
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
  library(pairwiseAdonis)


# Setting the core allocation to run parallel computations
options(mc.cores = parallel::detectCores())
cores <- detectCores(logical = T)

# Setting working directory
setwd("/your/path/to the file/here/")

# Reading the raw diet data - looking for the file "AC_dietitem_2024.csv"
# which is available on the 
Diet <- read.csv("You file here.csv", header = T, stringsAsFactors = T)

# Viewing and attaching the data
head(Diet)
attach(Diet)

# Making sure to set VBage as a factor
Diet$VBage <- as.factor(Diet$VBage)

# Use dplyr to pipe data into group_by cmd, and summarise data
# by individual Fish_ID - but leave site untouched - effectively 
# aggregates data by Fish_ID.
diet <- Diet %>%
  group_by(across(where(is.factor))) %>%
  summarize(across(where(is.numeric), ~ sum(.x, na.rm = T)))

# Making sure data have no NAs
any_missing <- apply(diet[, -c(1:3)], 2, function(x) any(is.na(x)))
any_missing

# Looking at the number of unique items in the different individuals
udi <- apply(diet[, -c(1:3)], 1, function(row) length(unique(row[row > 0])))

# Setting a sequence of numbers for the fish
ind <- seq(1,nrow(diet), by = 1)

# Combine with the udi, fish_ID, and VBage data from the diet data
ibi <- cbind.data.frame(indiv = ind, udi=udi, fishid=diet$Fish_ID, age=diet$VBage, site=diet$site)

############ calculate the mean and 95%CI of diet items in overall data excluding BB
sumibi<-ibi[ibi$site != "BB",]
mudi <- mean(sumibi$udi)
sdudi <- sd(sumibi$udi)
n <- nrow(sumibi)
SEM <- sdudi/sqrt(n)

# Determine the critical value for 99 degrees of freedom
# Use qt function to get the critical value for a 95% confidence interval
# The 95% confidence interval corresponds to 97.5th percentile for the two-tailed test
critval <- qt(0.975, df = n - 1)

# Calculate the margin of error
ME <- critval * SEM

# Calculate the confidence interval
CI_lower <- mudi - ME
CI_upper <- mudi + ME

# Print the results
CI <- c(CI_lower, CI_upper)
CI
##############

# Tabulate the data to get a feel for the numbers.
agedit <- table(ibi$age, ibi$udi)
agedit

# Make proportions table
ageditp <- prop.table(agedit,1)
ageditp
addmargins(ageditp)

# Plot the data to see the frequency of number of unique items in each fish
ggplot(ibi, aes(x=udi)) + 
  geom_bar(stat="count", fill = "skyblue3", colour = "black", linewidth = 1) +
  labs(x = "# of unique diet items", y = "Frequency") +
  scale_x_continuous(limits = c(-1, 8), breaks = seq(-1, 8, by = 1)) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 10)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 24),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"))


# Plots the data a bit differently to see what the mean and spread of  
# unique items per age class. This makes it look like we have a spread 
# of different diet items per age class, but this is a function of the jitter
# can be a bit misleading.
ggplot(ibi, aes(x=as.factor(age), y=udi)) +
  geom_jitter(color = "darkorchid4", size = 5, width = 0.2, 
              height = 0.1, alpha = 0.8, stroke = 1) +
  stat_summary(fun.data = mean_sdl, colour = "grey49", 
               na.rm = T, size = 1.5, linewidth = 1.5, alpha = 0.8, stroke = 1) +
  labs(x = "Age", y = "# of unique diet items") + 
  scale_y_continuous(limits = c(-1, 8), breaks = seq(-1, 8, by = 1)) +
  #scale_x_continuous(limits = c(1, 6), breaks = seq(1, 6, by = 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 24),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"))

## Setting age as a factor for remaining analyses
diet$VBage <- as.factor(diet$VBage)

# Use prop.table to calculate the proportions of diet items/ind 
# essentially devising volumetric-based data and making dp 
# (diet proportion dataframe)
dp <- prop.table(as.matrix(diet[,-c(1:3)]), 1)
dp <- replace(dp, is.nan(dp), 0.000)

## set the max number of point categories needed
ptsn<-max(udi)+1
ptsn

# set an empty vector
pts <- vector()

# Setting the first two values to 0 and 1, respectively
pts[1]<-0
pts[2]<-1

# Short loop to make the points for up to 8 categories (# unique items)
for (a in 3:ptsn) {
  pts[a]<-pts[a-1]*2
}

# Look at the points available per individual
pts <- rev(pts)
pts

# Function replacing non-zero elements of dp with ranks
replace_with_rank <- function(row) {
  nozeroel <- row[row != 0]
  ranks <- rank(-nozeroel)
  vals <- ranks
  row[row != 0] <- vals
  return(row)
}

# Apply the function "replace_with_rank" to the transposed
# matrix of dp (diet proportions) by row
dp2 <- t(apply(dp, 1, replace_with_rank))

# rename dp2 to dp3 for use in points allocation
dp3 <- dp2

# Assign the points to the ranks of dp2 starting with 1 = 64, 2 = 32,... 
for (b in 1:length(pts)) {
 dp3[dp2 == b] <- pts[b] 
}

# Look at and compare dp2 and dp3 to make sure
dp3

# replace the pts in dp3 to the calculated point proportion/ind
dp4 <- proportions(dp3, 1)

# replace any individuals that are nan to 0.000
dp4 <- replace(dp4, is.nan(dp4), 0.000)

# Combine dp4 matrix with Fish_ID, VBage, and site
dp5<-cbind.data.frame(fish_id=diet$Fish_ID, vbage=diet$VBage, site=diet$site, dp4)

# Remove age 5 and age 6 fish from anosim (no power on these)
dp6 <- dp5[dp5$vbage != "5",]
dp6 <- dp6[dp6$vbage != "6",]

# Check number of individuals in each age class and remove the 
# empty levels of the factor
table(dp6$vbage)
dp6$vbage <- droplevels(dp6$vbage) 

# Verify dropped levels
table(dp6$vbage)

# Remove individuals from "BB" site from analyses
dpbs <- dp6[dp6$site != "BB",]
dpbb <- dp6[dp6$site != "BS",]

# Focusing on BS
table(dpbs$vbage)

# Short function to get the mean of the columns of data.frame
bycolmean<-function(x) {
  tdf<-apply(x, 2, mean)
  return(tdf)
}

# Use the by cmd to split the data by age and then use
# the bycolmean to get means for each diet item (could do proportions too)
dietmeans<-by(dpbs[,-c(1:3)], dpbs$vbage, bycolmean)
dietmeans

# Return the data back into a single matrix/array
agediet <- do.call(cbind, dietmeans)
agediet

# Get the proportions of agediet in bs individuals
pagediet<-prop.table(agediet,2)
addmargins(pagediet)
# Convert to long format

## Set matrix as data frame to manipulate 
pagediet<-as.data.frame(pagediet)

## Set the new dataframe columns to specific age names
colnames(pagediet) <- c("1", "2", "3", "4")

## Set a new column to also have the names of all the diet items
pagediet$Row <- rownames(pagediet)

## Use tidyverse to restructure the agediet df to a long format 
plagediet <- pagediet %>%
  gather(key = "age", value = "diet", -Row)

# display the data.frame to make sure it appears as expected
plagediet

# Define your custom colors
cuscols <- c("darkslategray","dodgerblue","darkslateblue",
             "blue","powderblue","firebrick4",
             "gold4","darkorange", "yellow",
             "darkolivegreen1","skyblue3","aquamarine",
             "plum","bisque2","springgreen4",
             "salmon","darkorchid3","deeppink2","red")

# Produce a stacked barplot using ggplot but with colours set above 
# also make axes darker and enlarge writing. 

########## This is a figure that should appear in the analyses and shows the 
# mean proportions of diet items in the stomachs of differently aged AC from the 
# bs
ggplot(plagediet, aes(fill=Row, y=diet, x=age)) + 
  geom_bar(position="stack", stat="identity")+
  labs(x = "Age", y = "% volume of items in stomachs")+
  theme(aspect.ratio = 0.8)+
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, by = 0.2)) +
  scale_fill_manual(values = cuscols) +  # Use scale_fill_manual for custom colors
  annotate("text", x = c(1,3), y = c(1.05,1.05), label = paste("X"), color = "black", size = 12, family = "Courier") +
  annotate("text", x = 2, y = 1.05, label = paste("Y"), color = "black", size = 12, family = "Courier") +
  annotate("text", x = 4, y = 1.05, label = paste("XY"), color = "black", size = 12, family = "Courier") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 26),
        axis.text.x = element_text(size = 26),
        axis.title.y = element_text(size = 26),
        axis.text.y = element_text(size = 26),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom")

#################### Diet item analyses #########################################
# Section is greatly aided by https://rpubs.com/hafezahmad/948799

# Remove the rows (fish) that have empty stomachs
dpbb_ne <- dpbs[-which(rowSums(dpbs[,-c(1:3)]) == 0),]

# To use PERMANOVA/adonis, need to convert data to dissimilarity matrix
# use vegdist (vegan) to convert data to Bray-Curtis dissimilarity
acddis<-vegdist(dpbb_ne[,-c(1:3)], method='bray', na.rm = T)

# Run the pairwise.adonis2 to test pairwise differences in diets
# Stronger alternative to ANOSIM, which is not strong
ageclass_dd <- pairwise.adonis2(acddis ~ vbage, data = dpbb_ne, permutations = 10000,
                     sqrt.dist = T, add = T, parallel = cores)

# observe the pairwise results 
ageclass_dd

#### END
