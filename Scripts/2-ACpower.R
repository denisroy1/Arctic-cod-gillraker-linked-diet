# ACpower.R
# Short script used to calculate the power of our analyses and 
# test the adequacy of our sample sizes for them.
# In this case, we are interested in looking at the sample sizes 
# required to make inferences of gill raker differences among 4 age
# groups of Arctic cod (AC) from the Beaufort Sea. 

# written by DR - Summer 2023

## clearing all instances (data and variables)
rm(list = ls())

## loading the libraries needed 
{
library(pwr)
library(ggplot2)
}
## power test for one way anova, 
# k = # of groups, 
# sig. level = signficance level
# and f = effect sizes

# The smaller the effect size, the better, because it means that we would have 
# the power to detect small effect sizes if they exist. 

# Set a vector from 0-1 by 0.02 to estimate the change in effect size  
efsi <- seq(0.02,1,0.02)

# Set a vector from 1-50 to estimate the change in sample size
sasi <- seq(1,50,1)

# Set empty vectors to store the effect sizes and sample sizes
bsasi<-vector()
befsi<-vector()

# Short loop to calculate out the changes in sample size (bsasi) and 
# effect size (befsi) based on changing sasi and efsi. 
for (i in 1:50) {
  if (i < 3) { bsasi[i] == 0 && befsi[i] == 0.001 } else { 
  tsasi<-pwr.anova.test(k = 4, f = efsi[i], sig.level = 0.05, p = 0.9)
  bsasi[i]<-tsasi$n
  tefsi<-pwr.anova.test(k = 4, n = sasi[i], sig.level = 0.05, p = 0.8)
  befsi[i]<-tefsi$f
  }
}

# make a dataframe of the results for plotting
pdf<-cbind.data.frame(mss=sasi, mes=efsi, ess=bsasi, ees=befsi)

# remove the first few values as ss < 3 are not useful
pdf<-pdf[-c(1:3),]

# find the closest value to 0.35, which is a pretty good effect size
bees<- which.min(abs(pdf$ees - 0.35))

# find what these values are in sample size and estimated effect size
xseg<-pdf$mss[bees]
yseg<-pdf$ees[bees]

# plot the relationship with mss and ees to see where the best sample size 
# would be.
ggplot(pdf, aes(x = mss, y = ees)) +
  geom_line(linetype = "dashed") +
  geom_point(shape = 19, size = 4, colour = "skyblue3", alpha = 0.5) +
  labs(x = "Sample size", y = "Estimated effect size") +
  scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  geom_segment(aes(x = xseg, y = 0, xend = xseg, yend = yseg), linetype = "dotted", colour = "firebrick4") +
  geom_segment(aes(x = xseg, y = yseg, xend = 0, yend = yseg), linetype = "dotted", colour = "firebrick4") +
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

# Figure shows that aiming for a sample size around 20-25 would give a good chance to
# observe a difference with ~ 80% probability at the 0.05 level of significance.

