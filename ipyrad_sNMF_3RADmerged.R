################################################################################
## Script to run Population assignment in R using DAPC, IBD, and sNMF ##
################################################################################
## https://wyoibc.github.io/popgen_workshop/week4/DAPC_snmf/dapc_ibd_snmf.html 
## NOTE: Make sure you have all relevant outfiles after running Ipyrad 
################################################################################
################################################################################
## Install and Load relevant packages ##

#install.packages(c("plotrix", "rworldmap", "BiocManager", "vcfR", "fossil"))
## load up relevant packages for running DAPC, IBD, and sNMF 
library(adegenet)
library(LEA) # for running sNMF
library(plotrix)
library(mapdata)
library(rworldmap)
library(vcfR)
library(fossil)
library(MASS)
library(plyr)
library(tidyverse)
#install.packages("tidyverse")
library("Hmisc") # for imputing data for IBD
library(cowplot)
library(magrittr)
library(tidyr)
################################################
## For editing genid files
#install.packages("dartR")
#BiocManager::install("SNPRelate")
library(SNPRelate)
library("dartR")
################################################

## Packages for mapping
library(ggspatial)
library("rnaturalearth")
library("rnaturalearthdata")
library("sf")
library(tools)
library(ggplot2)  # ggplot() fortify()
#install.packages("ggplot2")
library(dplyr)  # %>% select() filter() bind_rows()
library(rgdal)  # readOGR() spTransform()
library(raster)  # intersect()
library(ggsn)  # north2() scalebar()
library(paletteer)
################################################################################
## Then we will specify a number of file paths and read in a few files so that 
## we don’t have to repeatedly hardcode file paths farther down in the script. 
## This makes it easier to reuse the script on different datasets or the same 
## data with different filtering schemes without having to search through the 
## script for every time an absolute file path is specified.
################################################################################

## Set up an object containing the path to the data ##

###########################################################################################################################
############################################ Ipyrad pipeline on 3RAD Merged data ########################################## 
###########################################################################################################################

###########################################################################################################################
#################################################### DATA DIRECTORIES #####################################################
###########################################################################################################################
## 3RADmerged: All samples from merged 3RAD_R1 and 3RAD_final - 296 individuals ##

  data_dir <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/allsamples"                        # Mac
  data_dir <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/allsamples"                       # PC
  
###########################################################################################################################
## noreponly: All 3RADmerged samples with no repeats - 273 individuals ##
  
  data_dir <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly"                         # Mac
  data_dir <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly"                        # PC
# VCF Filtered 
  vcffiles <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly/vcftools"
  
###########################################################################################################################
## norep_nolowloci: All 3RADmerged samples excluding repeated individuals and indiv. with low loci - 256 individuals ##
  
  data_dir <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/norepnolowloci"                    # Mac
  data_dir <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/norepnolowloci"                   # PC
  
###########################################################################################################################
## withgpsnorep: All 3RADmerged samples with GPS data excluding repeated individuals - 256 individuals ##
  
  data_dir <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/withgps"                           # Mac
  data_dir <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/withgps"                          # PC
  
###########################################################################################################################
## CMonly: Morphological Morelets from 3RADmerged samples excluding repeated individuals - 86 individuals ##
  
  data_dir <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/CMonly"                            # Mac
  data_dir <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/CMonly"                           # PC

###########################################################################################################################
## CMonly_struct: Morelets (<90%) determined from STRUCTURE using norepindiv dataset from Floyd - 103 individuals ##
  
  data_dir <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/CMonly_struct"                     # Mac
  data_dir <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/CMonly_struct"                    # PC
  
###########################################################################################################################
## allacutus: Morphological acutus from Cayes and Mainland from no repeated samples with GPS data - 122 individuals ##
  
  data_dir <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/allacutus"                         # Mac
  data_dir <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/allacutus"                        # PC

###########################################################################################################################
## acutus_struct: all acutus (<90%) determined from STRUCTURE using norepindiv dataset from Floyd - 95 individuals ##
  
  data_dir <- "/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/acutus_struct"                     # Mac
  data_dir <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/acutus_struct"                    # PC
  
###########################################################################################################################
###########################################################################################################################
## For filtered through snpFiltR data 
  filteredVCF <- paste0(data_dir, "/filteredVCF")
  filteredVCF
  data_dir <- filteredVCF
  data_dir
#####################################################################################
## Set as our working directory ##
setwd(data_dir)
getwd()
list.files()
#####################################################################################
## Set up object for Output Directory ##
#####################################################################################

out_dir <- paste0(data_dir, "/Pop_structr_out")

if(!dir.exists(out_dir)){ # check if the directory exists
  dir.create(out_dir)   # and create it if it does not
}

######################################################################################
## Set up an object that contains the base file name of files in the output directory. 
## Data files are all this basename with varying extensions
## we won't call this 'basename' because that is a function in R

# 3RADmerged: All samples from merged 3RAD_R1 and 3RAD_final - 296 individuals ##
basefile <- "3RADmerged"
  
# 3RADmerged_noreponly: All 3RADmerged samples with no repeats - 273 individuals ##
basefile <- "3RADmerged_noreponly"
  
# norep_nolowloci: All 3RADmerged samples excluding repeated individuals and indiv. with low loci - 256 individuals ##
basefile <- "norep_nolowloci"
  
# 3RADmerged_withgpsnorep: All 3RADmerged samples with GPS data excluding repeated individuals - 256 individuals ##
basefile <- "3RADmerged_withgpsnorep"
  
# 3RADmerged_CMonly: Morphological Morelets from 3RADmerged samples excluding repeated individuals - 86 individuals ##
basefile <- "3RADmerged_CMonly"

# CMonly_structpop: Morelets (<90%) determined from STRUCTURE using norepindiv samples excluding repeated individuals - 103 individuals ##
basefile <- "CMonly_structpop"

# allacutus_structpop_norepindiv: all acutus (<90%) determined from STRUCTURE using norepindiv samples excluding repeated individuals - 95 individuals
basefile <- "allacutus_structpop_norepindiv"   

######################################################################################
## Set up paths to input files using the base file name specified above ##
path_ugeno<-paste0(data_dir,"/", basefile,".ugeno")
path_ustr<-paste0(data_dir,"/", basefile,".ustr")
path_vcf<-paste0(data_dir,"/", basefile,".vcf")
path_coords<-paste0(data_dir,"/", basefile,"_coords.csv")
data_dir

## Path for filtered data 
path_vcf <-paste0(filteredVCF,"/", basefile,".filtered.85snp.vcf.gz") # 85% SNP
filtered_thinned_vcf <-paste0(filteredVCF,"/", basefile,".LDthinned.vcf")
path_vcf <- filtered_thinned_vcf
filtered_thinned_geno <-paste0(filteredVCF,"/", basefile,".LDthinned.geno")
path_ugeno <- filtered_thinned_geno
path_ugeno
######################################################################################
## Read in coordinates and metadata ##
coords<-read.csv(path_coords, header=TRUE, row.names=NULL)
coords

## for tidyverse/dplyr
coords.tbl <- read_csv(path_coords)
coords.tbl

## Get in MetaData ## 
md_csv<-paste0(data_dir,"/", basefile,"_metadata.csv")
md <- read_csv(md_csv)
names(md)
glimpse(md)
################################################################################
### CREATING METADATA CSV - if haven't already ###
#metadata <- read_csv("C:/Users/helen/Dropbox/Dissertation/DATA/3RAD_final/FINALMASTERSamples_Metadata.csv")
metadata <- read_csv("/Users/hwsung/Dropbox/Dissertation/DATA/3RAD_final/FINALMASTERSamples_Metadata.csv")

# Using metadata from norepindiv
metadata <- read_csv("/Users/hwsung/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly/3RADmerged_noreponly_metadata.csv")
metadata$Sample

# read in list of Sample names
sampleslist <- read_tsv("/Users/hwsung/Dropbox/Dissertation/DATA/Floyd_subset/unfiltered_subset/norepindiv_1000SNPS_subset_CApopstruct.txt", col_names = F)
names(sampleslist) <- "Sample"
sampleslist

# match metadata for only working samples 
md <- sampleslist %>%
  left_join(metadata, by = "Sample")
md 

# Fix up column headers
names(md)
names(md)[18] <- "CrocID"
names(md)[6] <- "Monitoring.Unit"

md_csv<-paste0(data_dir,"/", basefile,"_metadata.csv")
write.csv(md, file = md_csv, row.names = FALSE)

# make coords dataframe for IBD and mapping
coords.tbl <- pops %>% dplyr::select(CrocID, Sample, Longitude, Latitude)
coords.tbl

coords <- as.data.frame(coords.tbl)
coords
################################################################################
# make dataframe for 'population' later
glimpse(md)
md
pops <- md %>% dplyr::select(CrocID, Sample, Morph_Species, Longitude, Latitude, Monitoring.Unit, Route, Site, Subdivision)
pops
names(pops)

################################################################################
## Set up some colors for plotting farther down
colors_2<- c("green", "red") # colors for plotting 2 populations
colors_3 <- c("green", "red", "blue") # colors for plotting 3 populations
colors_4 <- c("green", "red", "blue", "darkred") # colors for plotting 4 populations

# Pick colors from paletteer package
colors_2 <- paletteer_d("RColorBrewer::Set1", 2)
colors_3 <- paletteer_d("RColorBrewer::Set1", 3)
colors_4 <- paletteer_d("RColorBrewer::Set1", 4)

################################################################################
################ Discriminant analysis of principal components #################
################################################################################
## Discriminant analysis of principal components, or DAPC, is a method of 
## classifying individuals into clusters that does not include any explicit 
## population genetic model. We'll use the unlinked Structure-formatted file for 
## this, but because Structure-formatted files can come in a variety of different 
## configurations, we need to tell the function how many individuals and loci 
## are present in the file. Read in the 'ugeno' file and use the dimensions of 
## that file to get the number of individuals and snps.
################################################################################

## read in the geno file to get the number of individuals and snps for this assembly
geno_txt<-readLines(path_ugeno)

## The number of lines is the number of loci
nums_snps<-length(geno_txt) 
nums_snps # 88559 loci --> 3RADmerged_noreponly
          # 87739 loci --> 3RADmerged_withgpsnorep


## the number of columns is the number of individuals
num_ind<-length(strsplit(geno_txt[[1]], "")[[1]]) # here we split apart the first line into individual characters and count the length
num_ind 

## Read in the Structure file ##
  ## A quirk of read.structure function is that it requires the structure file to have the
  ## file extension ".stru" - do some copying to make a new file with this extension

## Use a regular expression substitution to generate the new file name
path_stru<-gsub(".ustr", ".stru", path_ustr)
file.copy(path_ustr, path_stru) # make a copy of the file with the new name

## Now we can finally read in this file
DAPC_ustr<-read.structure(path_stru, n.ind=num_ind, n.loc=nums_snps, onerowperind = FALSE, col.lab=1, col.pop=0, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
DAPC_ustr ## look at how the data is structured for Adegenet
DAPC_ustr@pop <- as.factor(md$Morph_Species)

## Get the individual names in the order that they show up in the various files 
## this is important farther down for getting coordinates into the right order for plotting
ind_names<-rownames(DAPC_ustr@tab)
ind_names

## Run DAPC ##
  ## start by determining how many clusters we want to use, if best fit >= max.n.clust, 
  ## expand max number of clusters to ensure that we aren't artificially limiting 
  ## number of clusters that our data can fall into.

grp <- find.clusters(DAPC_ustr, max.n.clust=8) # test up to 8 clusters, 
  ## spits out plot of variance by each PC in the PCA transformation, 
  ## we want to retain all PCs up to point at which the cumulative variance plateaus, 
  ## if no clear plateau then retail all PCs 
      # input number of pc found = # of individuals 
  ## input number of clusters we want to retain, with the fit of each number of 
  ## clusters estimated by BIC. lowest BIC = best fit or look for inflection point 
  ### w/ sharpest decrease in BIC (the elbow in the BIC curve as function of K)
      # input = 2
grp3 <- find.clusters(DAPC_ustr, max.n.clust=8, n.pca = num_ind, n.clust = 3) # group by 3 instead of 2
grp4 <- find.clusters(DAPC_ustr, max.n.clust=8, n.pca = num_ind, n.clust = 4) # group by 4 instead of 2

## Output from find.clusters is a list
names(grp)
grp$Kstat  # What is the lowest kstat?
grp$stat
head(grp$grp, 10)
grp$size
table(pop(DAPC_ustr), grp$grp)
# Rows correspond to our morph group, columns to inferred groups 
table.value(table(pop(DAPC_ustr), grp$grp), col.lab=paste("inf_groups", 1:length(levels(grp$grp))),
            row.lab=levels(pop(DAPC_ustr)))

## grp object now contains the groupings of these individuals into the 2 clusters we decided on, and we can use the function dapc to describe these groups.
dapc1 <- dapc(DAPC_ustr, grp$grp) # run DAPC
  # input number of pc found = 275; 230; 256
  # We will then be asked how many discriminant functions to retain. 
      # With only 2 groups, only 1 is possible --> input number = 1
dapc2 <- dapc(DAPC_ustr, grp3$grp) # run DAPC with 2 discriminant functions from 3 groups
dapc3 <- dapc(DAPC_ustr, grp4$grp) # run DAPC with 3 discriminant functions from 4 groups

## membership probabilities based on the retained discriminant functions 
# stored in dapc objects in the slot "posterior" 
dapc1
class(dapc1$posterior)
dim(dapc1$posterior)
round(head(dapc1$posterior),2) # Each row corresponds to an individual, each column to a group.
summary(dapc1)

class(dapc2$posterior)
dim(dapc2$posterior)
round(head(dapc2$posterior),3) # Each row corresponds to an individual, each column to a group.
summary(dapc2)
round(dapc2$posterior,3)

dim(dapc3$posterior)
round(head(dapc3$posterior),3) # Each row corresponds to an individual, each column to a group.
summary(dapc3)
round(dapc3$posterior,3)
##########################################################
################# Plotting out everything ################
##########################################################

setwd(out_dir) # set to out directory 

## Plot out our confidence in assigning individuals to each of the groups w/dapc1 ##

## plot the DAPC the ugly way
scatter(dapc1, col=colors_2,  bg="white",
        legend=FALSE, posi.da = "bottomright",
        solid=.5)
scatter(dapc2, posi.da="bottomright",  bg="white",
        pch=17:22, cstar=0, col=colors_3, scree.pca=TRUE,
        posi.pca="bottomleft")
scatter(dapc3, posi.da = "bottomleft")
scatter(dapc3,1,1, col=colors_4, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)

## Plot from $assign.per.pop indicating proportions of successful reassignment (based on the discriminant functions) of individuals to their original clusters. 
  # Large values indicate clear-cut clusters, while low values suggest admixed groups.
  # Heat colors represent membership probabilities (red=1, white=0); blue crosses represent the prior cluster provided to DAPC.
  # DAPC classification is consistent with the original clusters (blue crosses are on red rectangles),
  # if blue cross on yellow, then there is discrepancy where its classified in red group but DAPC would assign to yellow group
  # useful when prior biological groups are used, as one may infer admixed or misclassified individuals.
assignplot(dapc1)
assignplot(dapc2)
assignplot(dapc3)

## To find most "admixed" individuals, consider admixed individual having no more than 90% prob of membership in single cluster
temp <- which(apply(dapc1$posterior,1, function(e) all(e<0.9)))
temp

## Plot in STRUCTURE-like way with with bar chart
compoplot(dapc1, posi="bottomright",
          txt.leg=paste("Cluster", 1:length(grp$size)), lab=rownames(dapc1),
          n.col=1, xlab="individuals", show.lab = TRUE, col = colors_2)

compoplot(dapc2, posi="bottomright",
          txt.leg=paste("Cluster", 1:length(grp3$size)), lab=rownames(dapc2),
          n.col=1, xlab="individuals", show.lab = TRUE, col = colors_3)

compoplot(dapc3, posi="bottomright",
          txt.leg=paste("Cluster", 1:length(grp4$size)), lab=rownames(dapc3),
          n.col=1, xlab="individuals", show.lab = TRUE, col = colors_4)

## Plot using spatial pattern on map ##
## first process the coordinates (remember that we read these into R earlier) 
## to make sure that the coordinates are in the right order and that every 
## individual we want to plot has locality data.

## We'll first do a little processing of the coordinates (remember that we read these into R earlier) 
##    to make sure that the coordinates are in the right order and that every individual we want to plot 
##    has locality data.

## make sure there aren't any individuals that don't have coordinates
ind_names[which(!ind_names %in% coords[,"Sample"])]

## match up the coordinates to the order of the individuals from genetic data
match_coords<-match(ind_names, coords[,"Sample"])
coords<-coords[match_coords,]

## Plot out the weird map I don't like
map("worldHires", "Belize", xlim=c(-89.5,-87.5), ylim=c(15.8,18.5),col="gray90", fill=TRUE)

## Add pie charts showing the probability of membership in each cluster for each sample
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc1$posterior[x,1], dapc1$posterior[x,2]), radius=0.1, col=colors_2)}

###############################################
## Better Map from my own script ##
###############################################

world <- getMap(resolution = "low")
world_belize <- world[world@data$ADMIN == "Belize", ]
world_belize
#bz <- CRS("+init=EPSG:2028")  #projected to Belize UTM
Belize1<-getData("GADM", country="BZ", level=1)
Belize1
class(Belize1)

plot(Belize1) # better map of Belize with district lines 

# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc1$posterior[x,1], dapc1$posterior[x,2]), radius=0.05, col=colors_2)}
# DAPC2
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc2$posterior[x,1], dapc2$posterior[x,2], dapc2$posterior[x,3]), radius=0.05, col=colors_3)}
# DAPC3
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc3$posterior[x,1], dapc3$posterior[x,2], dapc3$posterior[x,3], dapc3$posterior[x,4]), radius=0.05, col=colors_4)}

## MAKE PDF image of maps
DAPC_Map <-paste0("DAPC_Map_", basefile,".pdf")
DAPC_Map
pdf(file=DAPC_Map, width=10, height=10) # open the pdf plotting device
plot(Belize1) # Plot out the map
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x], 
                                        c(dapc1$posterior[x,1], dapc1$posterior[x,2]), radius=0.05, col=colors_2)}
dev.off() # close the pdf plotting device

DAPC_Map <-paste0("DAPC2_Map_", basefile,".pdf")
DAPC_Map
pdf(file=DAPC_Map, width=10, height=10) # open the pdf plotting device
plot(Belize1) # Plot out the map
# DAPC2
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc2$posterior[x,1], dapc2$posterior[x,2], dapc2$posterior[x,3]), radius=0.05, col=colors_3)}
dev.off() # close the pdf plotting device

DAPC_Map <-paste0("DAPC3_Map_", basefile,".pdf")
DAPC_Map
pdf(file=DAPC_Map, width=10, height=10) # open the pdf plotting device
plot(Belize1) # Plot out the map
# DAPC3
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(dapc3$posterior[x,1], dapc3$posterior[x,2], dapc3$posterior[x,3], dapc3$posterior[x,4]), radius=0.05, col=colors_4)}
dev.off() # close the pdf plotting device

## because this is a classification problem and isn’t explicitly modeling admixture, 
## we only see the confidence with which each sample is assigned to each cluster. 
## I.e., this does not indicate that there is no evidence of admixture among populations. 

## Map with ggplot ##

################################################################################
## DAPC Plot with population assignments ##

#https://luisdva.github.io/rstats/dapc-plot/
# Load libraries. Install first if needeed
library(readxl)    # CRAN v1.3.1
library(janitor)   # CRAN v2.1.0
library(dplyr)     # CRAN v1.0.7
library(tidyr)     # CRAN v1.1.3
library(ggplot2)   # CRAN v3.3.5
library(forcats)   # CRAN v0.5.1
library(stringr)   # CRAN v1.4.0
library(ggh4x)     # [github::teunbrand/ggh4x] v0.2.0.9000
library(paletteer) # CRAN v1.4.0
library(extrafont) # CRAN v0.17


# create an object with membership probabilities based on the retained discriminant functions
postprobs_dapc1 <- as.data.frame(round(dapc1$posterior, 3))  # Each row corresponds to an individual, each column to a group.
postprobs_dapc1

postprobs_dapc2 <- as.data.frame(round(dapc2$posterior, 3))  # Each row corresponds to an individual, each column to a group.
postprobs_dapc2

postprobs_dapc3 <- as.data.frame(round(dapc3$posterior, 3))  # Each row corresponds to an individual, each column to a group.
postprobs_dapc3

pops<- arrange(pops, Sample) # rearrange to match gendata
pops

# put probabilities in a tibble with IDS and labels for sites
clusters <- tibble::rownames_to_column(postprobs_dapc2, var = "Sample") %>%
  left_join(pops, by = "Sample")
clusters 

head(clusters)
head(coords)
head(postprobs_dapc2)

# melt into long format
croc_long <- clusters %>% pivot_longer(2:4, names_to = "cluster", values_to = "prob")
head(croc_long)
croc_long$cluster
# manual relevel of the sampling sites (to avoid alphabetical ordering)
croc_long$Morph_Species <- fct_relevel(as.factor(croc_long$Morph_Species), "CA", "CM", "HY", "MCA", "UK")
croc_long

# set up custom facet strips
facetstrips <- strip_nested(
  text_x = elem_list_text(size = c(12, 4)),
  by_layer_x = TRUE, clip = "off"
)

setwd(out_dir) # set to out directory 

pdf(file="DAPC2_barchart_MorphSp.pdf", width=10, height=10) # open the pdf plotting device
ggplot(croc_long, aes(factor(Sample), prob, fill = factor(cluster))) +
  geom_col(color = "gray", size = 0.01) +
  facet_nested(~ Morph_Species,
               switch = "x",
               nest_line = element_line(size = 1, lineend = "round"),
               scales = "free", space = "free", strip = facetstrips,
  ) +
  labs(x = "Individuals", y = "Membership probability") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  #scale_fill_paletteer_d("ghibli::PonyoMedium", guide = "none") +
  theme(
    panel.spacing.x = unit(0.18, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )
dev.off() # close the pdf plotting device
#######################################
## Different way of plotting##

#tidy the data for plotting. 
dapc_data_df <-
  # as_tibble() converts the ind.coord matrix to a special data frame called a tibble. 
  as_tibble(dapc2$ind.coord, rownames = "individual") %>%
  # mutate changes or adds columns to a data frame. Here, we're adding the population and group assignment columns to the data frame
  mutate(population = pops$Morph_Species,
         group = dapc2$grp)

dapc_data_df

#plot the data. Can color the points according to your pre-defined populations and the dapc groups to see if it conforms to your hypothesis.
dapc_plot <-
  ggplot(dapc_data_df, aes(
    x = LD1,
    y = LD2,
    fill = population
  )) +
  geom_point(shape = 21, size = 3) 

dapc_plot
################################################################################
############################ Isolation by distance #############################
################################################################################
## if isolation by distance (IBD) exists in our dataset, programs that seek to 
## cluster individuals but do not model continuous spatial structure can be 
## positively misled by IBD. This can result in overestimating the number of 
## population clusters, potentially identifying discrete population structure 
## when no such structure exists. The method conStruct can explicitly model both 
## processes, but can be finicky to run and have long run times. Do quick test 
## and visualization of isolation by distance to try to determine if IBD is
## misleading our population structure analyses.
################################################################################

gendata_all<-read.vcfR(path_vcf) # read in all of the genetic data
gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
gendata_names <- indNames(gendata) # get the sample names
gendata_names

## make sure there aren't any individuals that don't have coordinates
gendata_names[which(!gendata_names %in% coords[,"Sample"])]
coords

## match up the coordinates to the order of the individuals from genetic data
match_coords<-match(gendata_names, coords[,"Sample"])
coords<-coords[match_coords,]
coords

################################################################################
## Identify individuals missing geographic data 
missing.ind <- coords[rowSums(is.na(coords)) > 0, ]   
missing.ind <- missing.ind$Sample 
missing.ind

# Create new dataframe without individuals missing geographic data
coords_nomissing<- filter(coords.tbl, !(Sample %in% missing.ind))
coords_nomissing <- as.data.frame(coords_nomissing)
coords_nomissing

# Drop individuals with missing geographic data
gendata@ind.names
gendata
missing.ind.list <- as.list(missing.ind)
gl2 <- gl.drop.ind(gendata, ind.list=missing.ind.list,recalc=TRUE)

################################################################################
## Run Mantel test (correlates two different distance matrices); need to convert DNA and geographic data into pairwise distances for this test
Dgeo<-earth.dist(coords[,c("Longitude", "Latitude")]) # get the geographic distances
Dgeo<-earth.dist(coords_nomissing[,c("Longitude", "Latitude")]) # get the geographic distances with no missing individuals

Dgen<-dist(gendata) # get the genetic distances
Dgen<-dist(gl2) # get the genetic distances with missing individuals

ibd<-mantel.randtest(Dgen,Dgeo) # run the mantel test
ibd

## Impute data for missing data if error message: Error in testmantel(nrepet, col, as.matrix(m1), as.matrix(m2)) : NA/NaN/Inf in foreign function call (arg 3)
Dgen1 <- impute(Dgen, "random")
ibd<-mantel.randtest(Dgen1,Dgeo) # run the mantel test with imputed data 

## We can then visualize the result of our empirical estimate of isolation by distance 
  ## compared to a permuted null distribution to see how significant our result is. 
  ## We'll also plot out a kernel density plot of genetic vs. geographic distances to 
  ## visualize how these distances are associated.

## PDF of mantel output and a kernel density plot of genetic vs. geograophic distance
Mantel_KD <-paste0("Mantel_KD_", basefile,".pdf")
pdf(file=Mantel_KD, width=8, height=8)
  plot(ibd, main=paste0("mantel p = ", ibd$pvalue)) # plot out the IBD significance
  ## make kernel density plot of genetic and geographic distances
  dens <- kde2d(as.numeric(Dgeo),as.numeric(Dgen), n=300)
  myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
  plot(Dgeo, Dgen, pch=20,cex=.5, xlab="Geographic Distance", ylab="Genetic Distance")
  image(dens, col=transp(myPal(300),.7), add=TRUE)
  abline(lm(as.numeric(Dgen)~as.numeric(Dgeo)))
  lines(loess.smooth(Dgeo, Dgen), col="red")
  title("IBD plot")
dev.off()

## Plot With imputation
Mantel_KD <-paste0("Mantel_KD_withimputation_", basefile,".pdf")
pdf(file=Mantel_KD, width=8, height=8)
  plot(ibd, main=paste0("mantel p = ", ibd$pvalue)) # plot out the IBD significance
  ## make kernel density plot of genetic and geographic distances
  dens <- kde2d(as.numeric(Dgeo),as.numeric(Dgen1), n=300)
  myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
  plot(Dgeo, Dgen1, pch=20,cex=.5)
  image(dens, col=transp(myPal(300),.7), add=TRUE)
  abline(lm(as.numeric(Dgen1)~as.numeric(Dgeo)))
  title("IBD plot")
dev.off()
## Result --> 2-page pdf of these plots ##
  ## Page 1: we can see that we have highly significant isolation by distance, 
    ## with the lowest possible significance value given our number of permutations 
    ## in the Mantel test (999 by default). The diamond-shaped point with the line
    ## indicates our empirical estimate, with the permutations shown as the gray histogram.
  ## Page 2: If we look at the second page, we can see how the geographic and genetic 
    ## distances are correlated. What we see is a positive relationship, but with 
    ## a major disjunction. This type of disjunction indicates the presence of some 
    ## level of not fully continuous spatial genetic structure.

################################################################################
################################## sNMF ########################################
################################################################################
## The snmf function requires a geno file as input, and requires that it has the 
## extension .geno. We want to use only unlinked SNPs here (i.e., 1 SNP per RAD 
## locus, assumed to be unlinked), and the geno file of unlinked snps has the 
## extension ugeno, so we’ll copy the file and give it a new extension
################################################################################

##########################
# Manage an snmf project #
##########################

# All the runs of snmf for a given file are 
# automatically saved into an snmf project directory and a file.
# The name of the snmfProject file is the same name as 
# the name of the input file with a .snmfProject extension 
# ("genotypes.snmfProject").
# The name of the snmfProject directory is the same name as
# the name of the input file with a .snmf extension ("genotypes.snmf/")
# There is only one snmf Project for each input file including all the runs.

# An snmfProject can be load in a different session.
project_path <- paste0(basefile, ".u.snmfProject")
#project_path <- "C:/Users/helen/Dropbox/Dissertation/DATA/3RADmerged/ipyrad_results/noreponly/3RADmerged_noreponly.u.snmfProject"
obj.at = load.snmfProject(project_path) 
################################################################################

## Use a regular expression substitution to generate the new file name
path_geno<-gsub(".ugeno", ".u.geno", path_ugeno)

file.copy(path_ugeno, path_geno) # do the copying with the new name

## Now we're ready to run sNMF. We'll run this using 1 to 10 possible ancestral 
## populations and evaluate the fit of these different numbers of populations 
## (referred to as k values) to the data using the cross entropy criterion.
obj.at <- snmf(input.file = path_geno,  # input file is the .geno format file
               K = 1:8, # we will test for k=1 through 8
               ploidy = 2, 
               entropy = T, # use the cross entropy criterion for assessing the best k 
               repetitions = 10, # Run 10 independent replicate analyses
               CPU = 4, 
               project = "new", tolerance = 0.00001, iterations = 10000)

## Let's make a pdf of the cross-entropy plot ##
snmf_cross_ent <-paste0("snmf_cross_ent_", basefile,".pdf")
pdf(snmf_cross_ent, width = 8, height=5)
  plot(obj.at, col = "lightblue", cex = 1.2, pch = 19)
dev.off()

## As for DAPC, the best fit model is one with 2 populations, as shown by the 
## lowest cross-entropy score. We can also look at a numeric summary of this result:
outstats <- summary(obj.at)
outstats # best fit is 2 populations shown by lowest cross-entropy score
  # the value for which the function plateaus or increases is our estimate of K 

## We can also confirm cross entropy values for K are consistent across runs and get the single best run for K=2
(ce <- cross.entropy(obj.at, K = 2))
ce ## best ce is run 6 at 0.3860456
(best.run <- which.min(ce)) # find the run with the lowest cross validation error

# Entropy value corresponding to best run for best K
e = round(ce[best.run], digits = 4) # rounding to four digits
e

(ce_3 <- cross.entropy(obj.at, K = 3))
(best.run_3 <- which.min(ce_3)) # find the run with the lowest cross validation error

(ce_4 <- cross.entropy(obj.at, K = 4))
(best.run_4 <- which.min(ce_4)) # find the run with the lowest cross validation error

## Then we can get the snmf Q matrix from the best run at the best k, which is 
## a matrix of the proportion of ancestry that each sample derives from each population.
qmatrix <- Q(obj.at, K = 2, run = best.run)
qmatrix_K3 <- Q(obj.at, K = 3, run = best.run_3)
qmatrix_K4 <- Q(obj.at, K = 4, run = best.run_4)

# cluster assignment for each individual
cluster<- apply(qmatrix, 1, which.max) #this corresponds with the 1:24 order 
cluster
cluster3<- apply(qmatrix_K3, 1, which.max) #this corresponds with the 1:24 order 
table(cluster3)

colnames(qmatrix) <- paste0("P", 1:2)
colnames(qmatrix_K3) <- paste0("P", 1:3)
colnames(qmatrix_K4) <- paste0("P", 1:4)
qmatrix

## Convert the qmatrix into a dataframe to order based on ancestry ##
# Meh way
admix<-as.data.frame(qmatrix)
admix_K3<-as.data.frame(qmatrix_K3)
admix_K4<-as.data.frame(qmatrix_K4)

#sort by mpg (ascending) and cyl (descending)
admix_ordered <- admix[order(admix$P1),]
admix_ordered3 <- admix_K3[order(admix_K3$P1, admix_K3$P2),]
admix_ordered3 <- admix_K3[order(admix_K3[,1], -admix_K3[,3] ),]
admix_ordered4 <- admix_K4[order(admix_K4[,1], -admix_K4[,4]),]

## Better way
# convert the q matrix to a data frame
q_df <- qmatrix %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(Sample = pops$Sample,
  )
q_df
q_df$P2
which.max(q_df$P2)

# Order data by decending P1 and then ascending P2 for graphing
q_df <- arrange(q_df, desc(P1)) 
q_df <- mutate(q_df, order = 1:num_ind)
q_df

# Create csv file for qmatrix and metadata
q_mat <- q_df %>% left_join(pops, by = "Sample")
q_mat

qmatrix_csv<-paste0(out_dir,"/", basefile,"_qmatrix.csv")
qmatrix_csv
write.csv(q_mat, file = qmatrix_csv, row.names = FALSE)

## Create qmatrix dataframe for other K values 
# convert the q matrix to a data frame
q_df3 <- qmatrix_K3 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(Sample = pops$Sample,) 
q_df3

# Order data by decending P1 and then ascending P2 for graphing
q_df3 <- arrange(q_df3, desc(P1)) 
q_df3 <- mutate(q_df3, order = 1:num_ind)
q_df3

# Create csv file for qmatrix and metadata
q_mat3 <- q_df3 %>% left_join(pops, by = "Sample")
q_mat3

qmatrix3_csv<-paste0(out_dir,"/", basefile,"_qmatrix3.csv")
write.csv(q_mat3, file = qmatrix3_csv, row.names = FALSE)

q_df4 <- qmatrix_K4 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(Sample = pops$Sample,) 
q_df4

# Order data by decending P1 and then ascending P2 for graphing
q_df4 <- arrange(q_df4, desc(P1)) 
q_df4 <- mutate(q_df4, order = 1:num_ind)
q_df4

##########################################
## PLOTTING RESULTS: Map ##
##########################################
coords
pops

## Plot results onto a map like we did for our DAPC results ##
sNMF_Map_K2 <-paste0("sNMF_Map_K2_", basefile,".pdf")
pdf(file=sNMF_Map_K2, width=10, height=10) # open the pdf plotting device
# Plot out the map
# map("worldHires", "Belize", xlim=c(-89.5,-87.5), ylim=c(15.8,18.5),col="gray90", fill=TRUE)
plot(Belize1)
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix$P1[x],admix$P2[x]), radius=0.05, col=colors_2)}
dev.off() # close the pdf plotting device
## The majority population membership matches up with what we saw from DAPC, 
## but we are now able to see the admixture present where the populations contact 
## each other.

## SNMF plot for K = 3 ##
sNMF_Map_K3 <-paste0("sNMF_Map_K3_", basefile,".pdf")
pdf(file=sNMF_Map_K3, width=10, height=10) # open the pdf plotting device
# Plot out the map
# map("worldHires", "Belize", xlim=c(-89.5,-87.5), ylim=c(15.8,18.5),col="gray90", fill=TRUE)
plot(Belize1)
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix_K3$P1[x],admix_K3$P2[x],admix_K3$P3[x]), radius=0.05, col=colors_3)}
dev.off() # close the pdf plotting device

## SNMF plot for K = 4 ##
sNMF_Map_K4 <-paste0("sNMF_Map_K4_", basefile,".pdf")
pdf(file=sNMF_Map_K4, width=10, height=10) # open the pdf plotting device
# Plot out the map
# map("worldHires", "Belize", xlim=c(-89.5,-87.5), ylim=c(15.8,18.5),col="gray90", fill=TRUE)
plot(Belize1)
# Plot pies at each locality
for (x in 1:nrow(coords)) {floating.pie(coords$Longitude[x],coords$Latitude[x],
                                        c(admix_K4$P1[x],admix_K4$P2[x],admix_K4$P3[x], admix_K4$P4[x]), radius=0.05, col=colors_4)}
dev.off() # close the pdf plotting device

##########################################
## PLOTTING RESULTS: stacked barchart ##
##########################################
## Basic ggplot stacked barchart ##
# Example                                                                                                                                    -93L))
tbl # example data

plot_data <- tbl %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
plot_data

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  theme_classic()

ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  facet_grid(~likely_assignment, scales = 'free', space = 'free')

ggplot() + 
  geom_col(data = q_df, aes(x=order, ))
############################
# Now, you're transforming the data to a "long" format for plotting. The population column names get their own column and the ancestry proportions (q) get their own column.  
q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 
q_df_long

which.max(q_df_long$q)

q_df_ordered <- q_df_long %>%
  # assign the population assignment according to the max q value (ancestry proportion) and include the assignment probability of the population assignment
  group_by(Sample) %>%
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>%
  # arrange the data set by the ancestry coefficients
  arrange(likely_assignment, desc(assignment_prob)) %>%
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  ungroup() %>%
  mutate(Sample = forcats::fct_inorder(factor(Sample)))
q_df_ordered

q_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = order, y = q, fill = pop)) +
  scale_fill_manual(values = colors_2, labels = c("CA", "CM")) +
  #scale_fill_viridis_d() +
  labs(fill = "Struc_pop") +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
  ) 
q_df_ordered

admix_ordered4


## Plot bar chart like standard structure plot 
barplot(t(admix_ordered),col = colors_2, border = NA, space = .2, xlab = "Individuals", ylab = "Admixture coefficients", main = "Ancestry matrix", horiz = FALSE, names.arg = q_df$Sample ,cex.names=0.4, las = 2)
barplot(t(admix_ordered3), col = colors_3, border = NA, space = .2,xlab = "Individuals", ylab = "Admixture coefficients", main = "Ancestry matrix", names.arg = q_df3$Sample,  cex.names=0.4, las = 2)
barplot(t(admix_ordered4), col = colors_4, border = NA, space = .2,xlab = "Individuals", ylab = "Admixture coefficients", main = "Ancestry matrix", names.arg = q_df4$Sample,  cex.names=0.4, las = 2)


## create pdf of barchart
sNMF_boxplot_K2 <-paste0("sNMF_barplot_K2_", basefile,".pdf")
pdf(file=sNMF_boxplot_K2, width=10, height=10) # open the pdf plotting device
barplot(t(admix_ordered), col = colors_2, border = NA, space = .2, xlab = "Individuals", ylab = "Admixture coefficients", main = "Ancestry matrix", horiz = FALSE, names.arg = q_df$Sample ,cex.names=0.4, las = 2) 
dev.off() # close the pdf plotting device

sNMF_boxplot_K3 <-paste0("sNMF_barplot_K3_", basefile,".pdf")
pdf(file=sNMF_boxplot_K3, width=10, height=10) # open the pdf plotting device
barplot(t(admix_ordered3), col = colors_3, border = NA, space = .2,ylab = "Admixture coefficients", main = "Ancestry matrix", cex.names=0.4, las = 2, xaxt="n")
dev.off() # close the pdf plotting device

sNMF_boxplot_K4 <-paste0("sNMF_barplot_K4_", basefile,".pdf")
pdf(file=sNMF_boxplot_K4, width=10, height=10) # open the pdf plotting device
barplot(t(admix_ordered4), col = colors_4, border = NA, space = .2,ylab = "Admixture coefficients", main = "Ancestry matrix", cex.names=0.4, las = 2, xaxt="n")
dev.off() # close the pdf plotting device

#####################################################################################
##Some formatting crap## 
q_df <- q_df %>% dplyr::select(P1, P2, P3, ID, Sample, order, Longitude, Latitude, Monitoring.Unit, Route, Site, Subdivision, Species)
head(q_df)
as_tibble(q_df)
glimpse(q_df)


q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q") 

glimpse(q_df_long)

q_df_prates <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(ID = forcats::fct_inorder(factor(ID)))

q_df_prates

q_df_ordered <- q_df_long %>% 
  # assign the population assignment according to the max q value (ancestry proportion) and include the assignment probability of the population assignment
  group_by(ID) %>%
  mutate(likely_assignment = pop[which.max(q)],
         assignment_prob = max(q)) %>%
  # arrange the data set by the ancestry coefficients
  arrange(likely_assignment, assignment_prob) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  ungroup() %>% 
  mutate(ID = forcats::fct_inorder(factor(ID)))

q_df_ordered

# a custom palette for plotting
q_palette <- c("#fde725", "#35b779", "#440154")

q_df_ordered %>% 
  ggplot() +
  geom_col(aes(x = ID, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette, labels = c("CA", "CM", "HY")) +
  #scale_fill_viridis_d() +
  labs(fill = "Species") +
  theme_minimal() +
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
  )
#####################################################################################
#####################################################################################
## Organized barplot by ancestral coeff ##
# Adding ID names for plotting
individuals = pops$CrocID
q_mat$ID = individuals
q_mat$Sample = pops$Sample
head(q_mat)
q_mat %>% add_column(plot_order = 1:273)
pops

# Arrange order or samples in plot
df <- arrange(coords_md, desc(CrocID))
df <- df %>% add_column(plot_order = 1:273)
df$plot_order
df
admix$Samples
df$Sample

## match up the coordinates to the order of the individuals from genetic data
match_names<-match(admix$Samples, df$Sample)
df_reordered <- df %>% slice(match(admix$Samples, df$Sample))
df_reordered$plot_order

admix$plot_order = df_reordered$plot_order
admix_ord = plyr::arrange(admix, plot_order)

admix = as_tibble(admix) # dplyr likes the tibble format
admix <- q_df
admix_ord$ID = factor(admix_ord$ID, levels = unique(admix_ord$ID)) # transforming ID into factor to keep order in plot
admix_ord = as.data.frame(admix_ord)

# Save qmatrix
write.csv(admix, file = paste0("qmatrix_", basefile, ".csv"), quote = FALSE, row.names = FALSE)
write.csv(q_df, file = paste0("qmatrix_", basefile, ".csv"), quote = FALSE, row.names = FALSE)

# "Melt" dataframe using the gather function, as required by ggplot
admix_melt = gather(admix, key = cluster, value = coeff, 1:2)

# Plot with ggplot
SNMF_plot = admix_melt %>% ggplot(aes(x= order())) +
  geom_bar(aes(x= order, y = coeff, fill = cluster), stat = "identity", position = "fill") +
  #scale_fill_manual("Population", values = myPalette[c(1:bestK)], labels = c("Atlantic Forest", "Amazonia 1", "Amazonia 2")) + 
  #scale_fill_manual("Population", values = myPalette[c(1:bestK)], labels = c("Amazonia", "Atlantic Forest")) + 
  scale_color_manual(values = colors_2, guide = "none") +
  theme_minimal(base_size = 16) +
  labs(y = "Ancestrality coefficient", x = "") +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 5),
        panel.grid.major.x = element_blank(),
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0))

SNMF_plot

#####################
## Diff method ## 
#####################

# plot the ancestry coefficients for the best run and K = 9
library(MetBrewer)

pdf(file="sNMF_barchart.pdf", width=10, height=10) # open the pdf plotting device
bp <- barchart(obj.at, K = 2, run = best.run, sort.by.Q = TRUE,
         border = NA, space = 0, col = colors_2, 
         xlab = "Individuals", ylab = "Ancestry proportions", 
         main = "Ancestry matrix")

axis(1, at = 1:length(bp$order), 
     labels = bp$order, las = 3, 
     cex.axis = .4)
dev.off() # close the pdf plotting device
bp

b <- c(1:length(gendata_names))
b <- as.numeric(b)
df <- cbind(gendata_names, b)
df

match_df<-match(bp$order, df[,b])
match_df
df <- df[match_df,]
df
coords<-coords[match_coords,]
coords
coords1 <- na.omit(coords)
coords1

# An snmfProject can be exported to be imported in another directory
# or in another computer

export.snmfProject("dataset3.u.snmfProject")

########################################################################################
## Genome scan for selection using an ancestral allele frequency differentiation test ##
########################################################################################
## Function for FST statistic computed from the ancestral allele frequencies estimated by snmf.
fst = function(project,run, K, ploidy = 2){
  library(LEA)
  l = dim(G(project, K = K, run = run))[1]
  q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
  if (ploidy == 2) {
    G1.t = G(project, K = K, run = run)[seq(2,l,by = 3),]
    G2.t = G(project, K = K, run = run)[seq(3,l,by = 3),]
    freq = G1.t/2 + G2.t}
  else {
    freq = G(project, K = K, run = run)[seq(2,l,by = 2),]}
  H.s = apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x))
  P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x))
  H.t = P.t*(1-P.t)
  return(1-H.s/H.t)
}

## Compute the FST stat
fst.values = fst(obj.at,run = best.run, K=2)

## Convert the FST values into absolute values of z-scores 
n = dim(qmatrix)[1]
n
fst.values[fst.values<0] = 0.000001
K = 2
z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
z.scores

## Computing a correct set of p-values ##
## The estimated z-scores must be recalibrated before applying a test for 
## neutrality at each locus. Due to the existence of population structure, the 
## null-hypothesis may be misspecified, and this step is often necessary.
## We suggest using an “empirical-null hypothesis” approach, starting with 
## computing a genomic inflation factor, λ, and then dividing the scores by λ
#Compute the GIF
K=2
lambda = median(z.scores^2)/qchisq(1/2, df = K-1)
lambda  # 1.563627

# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(z.scores^2/lambda, df = K-1, lower = FALSE)

#histogram of p-values
hist(adj.p.values, col = "red")
  # Then, we can check the histogram of p-values. Ideally, this histogram should have 
  # a flat shape with a peak close to zero. Having this shape indicates that the 
  # null-hypothesis is correct, i.e., p-values are sampled from a uniform 
  # distribution under the null-hypothesis.

## The genomic inflation factor is acknowledged to be overly conservative. 
## So, it might be better to calibrate the p-values using other methods 
## (eg, using the R package “fdrtool") or following a trial-and-error approach.
## Here we try λ=1.4.

# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(z.scores^2/1.4, df = 1, lower = FALSE)
adj.p.values = pchisq(z.scores^2/lambda, df = 1, lower = FALSE)

#histogram of p-values
hist(adj.p.values, col = "green")
  # plot looks okay

## Control of false discoveries ##

## FDR control: Benjamini-Hochberg at level q
L = nums_snps # 88559 loci --> 3RADmerged_noreponly
q = 0.05
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates = order(adj.p.values)[w]

## our list of candidate loci is recorded in the R object “candidates” ##
plot(-log10(adj.p.values), main="Manhattan plot", xlab = "Locus", cex = .7, col = "grey")
points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")

########################################################################################
## Estimate effective population size using LD ##
########################################################################################
gl.LDNe(gendata, outfile = "CM_genepopLD.txt", )

nes <- gl.LDNe(gendata, "CM_genepopLD.txt", outpath=tempdir(),
               neest.path = "./path_to Ne-21",
               critical=c(0,0.05), singleton.rm=TRUE, mating='random')

