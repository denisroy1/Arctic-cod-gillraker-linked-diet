#### Mapping Script ####
# Beaufort Sea sampling sites for Launay et al. 2025
# This short script uses ggOceanMaps written by Mikko to visualise 
# our sampling sites used for the Gill raker/diet data analyses.

citation("ggOceanMaps")

# Written by dr, 2025-08-27

#### 1. Working Environment ####

# Clearing all instances
rm(list = ls())

# Loading graphics libraries to use
{
  library(ggOceanMaps)
  library(ggspatial) # for data plotting
}

# Set working directory to get the associated files. Mostly interested
# in getting the sampling site positions, and where to plot them on the
# figure(s).
setwd("Enter/your path/to the file/here")

#### 2. Data import ####

# Read in sample site coordinates from a file.File also contains
# lats longs of where the names should appear in the produced map.
sscoor <- read.csv("your map coords here.csv", header = T, stringsAsFactors = T)

# Assign bslabs a dataframe housing the information for the 
# Beaufort Sea labels. Adjust X positions, adjust Y positions, and text you want on the map.
bslabs <- data.frame(Long = c(-122.5, -141.6), Lat  = c(75, 71), name = c("Beaufort\nSea", "Beaufort\nSea"))

# Assign land as a dataframe housing the information for the 
# major land features.
land <- data.frame(Long = c(-121.9, -136.5,-112, -121), Lat = c(73, 68.3, 71.2, 68.7), name = c("Banks\nIsland", "Northwest Territories", 
           "Victoria\nIsland", "Nunavut"))

# These cmds change the names of the sites with long names to appear 
# on more than a single line.
sscoor$site <- gsub("MacKenzie Shelf Slope", "MacKenzie\nShelf\nSlope", sscoor$site)
sscoor$site <- gsub("MacKenzie Shelf", "MacKenzie\nShelf", sscoor$site)

#### 3. Mapping ####

# Make the basemap as per directions in ggOceanMaps outlining the Beaufort Sea 
# region. Add to it the sample site locations and labels, and the major land
# features. Also include a north-pointing arrow, a scale, and a legend for 
# the depth profile
basemap(limits = c(-145, -110, 68, 76), shapefiles = "Arctic", rotate = TRUE, 
        bathymetry = TRUE, bathy.style = "rcg")+
  geom_spatial_point(data = sscoor, aes(x = Long, y = Lat), 
                     colour = "firebrick4", size = 4, stroke = 1, alpha = 0.8)+
  geom_spatial_text(data = sscoor, aes(x = Long_lab, y = Lat_lab, label = site),
                    size = 5.85, fontface = "bold", colour = "firebrick4", lineheight = 0.9) +
  geom_spatial_text(data = bslabs, aes(x = Long, y = Lat, label = name),
                    size = 7, fontface = "italic", colour = "black") +
  geom_spatial_text(data = land, aes(x = Long, y = Lat, label = name),
                    size = 7.8, fontface = "bold.italic", colour = "black") +
 annotation_scale(location = "br", text_cex = 1.5) + 
 annotation_north_arrow(location = "tr", which_north = "true")+
  theme(axis.text.x = element_text(size = 17),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 17),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm")) +
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))

#### 4. Inset #### 
# Making the inset map in a similar way as above. Here, the area is larger 

# Initiate a section of the map that will be the focus of the main map.
inset <- data.frame(lon = c(-145, -145, -110, -110), lat = c(67, 76, 76, 67))

# As above, set land2 as a data frame outlining the coords of the placement 
# of the major land features.
land2 <- data.frame(Long = c(-112, -40), Lat = c(63, 73.5), name = c("Canada", "Greenland"))

# Here, replot a larger-scale map which will be used as the inset to our 
# main map.
basemap(limits = c(-145, -40, 55, 85), shapefiles = "Arctic", rotate = TRUE, 
        bathymetry = TRUE, bathy.style = "rcg")+
  annotation_scale(location = "bl", text_cex = 3) + 
  annotation_north_arrow(location = "tr", which_north = "true")+
  ggspatial::geom_spatial_polygon(data = inset, aes(x = lon, y = lat), 
                                  fill = NA, color = "red", size = 3)+
  geom_spatial_text(data = land2, aes(x = Long, y = Lat, label = name),
                    size = 20, fontface = "bold.italic", colour = "black") +
  theme(axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title = element_text(face = "bold"),
        axis.line = element_line(linewidth=1.3),
        axis.ticks = element_line(linewidth = 1.3),
        axis.ticks.length = unit(.3,"cm")) +
  theme(legend.text = element_text(size = 12))

# The two maps made should be copied to the clipboard and expanded to 
# a width of ~ 14-1600 and a height of 800 (maintaining ratios).

# The large-scale map can be cropped and placed in the top right corner
# of the main map as the inset.

#### 5. END ####
