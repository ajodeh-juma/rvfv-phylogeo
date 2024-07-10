#!/usr/bin/env Rscript

#'@author John Juma
#'
#'
# clear the working environment and free up memory
rm(list=ls(all.names = TRUE))
gc()


# set options
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
options(timeout = 1000)

###########################################
######## load/install packages ############
###########################################

library(lubridate)
library(dplyr)
library(countrycode)
library(colorspace)
library(scales)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggridges)
library(treeio)
library(patchwork)
library(sf)
library(tidyverse)
library(maps)
library(scatterpie)
library(ggthemes)
library(ggimage)
library(sp)
library(hrbrthemes)
library(ggrepel)
library(showtext)
# additional fonts
# check the current search path for fonts
font_paths()
# List available font files in the search path
font_files()  
# Add desired font
font_add("Arial Narrow", "/System/Library/Fonts/Supplemental/Arial Narrow.ttf")
# list loaded fonts
font_families()

# get path of the script
if (rstudioapi::isAvailable()){
  if (require('rstudioapi') != TRUE){
    install.packages('rstudioapi', repos = r)
  } else {
    library(rstudioapi)
  }
  wd <- dirname(getActiveDocumentContext()$path)
} else {
  wd <- getwd()
}


# output directories
figuresDir <- file.path(wd, "figures")
if (!dir.exists(figuresDir)) dir.create(figuresDir, recursive = TRUE)

setwd(wd)

# colours
cols_countries <- c("#EEEE00", "#5D5E36", "#732039", "#17664D", "#30CED8",
                    "#4E11A8", "#E88817", "#FF85A9", "#666666", "#91793C",
                    "#DC143C", "#EE00EE", "#F7D2BC", "#FF69B4", "#AFEEEE",
                    "#9A3324", "#6495ED", "#12B385")

cols_lineages <- c("#9370DB", "#30CED8", "#FF0000", "#82A3FF", "#000033",
                   "#567391", "#FF00B6", "#E88817", "#009FFF", "#FE5000",
                   "#00FFBE", "#8B5A2B", "#36312E", "#BB85AB", "#B1CC71")

cols_hosts <- c("#655454", "#B3C4DD", "#D6C499", "#947059", "#BAA844", 
                "#ffb302", "#606024")


country.codes <-  c("BDI", "BFA", "CAF", "EGY", "GIN", "KEN", 
                   "MDG", "MRT", "MYT", "NAM", "RWA", "SAU", 
                   "SDN", "SEN", "TZA", "UGA", "ZAF", "ZWE")
countrynames <- countrycode(country.codes, origin = 'iso3c', 
                            destination = 'country.name')
hosts <- c("buffalo", "cow", "goat", "human", "mosquito", "sheep", 
           "springbok")
hosts <- str_to_title(hosts)
df_countries <- data.frame(country = country.codes,
                           countryname = countrynames,
                           color = cols_countries)
df_lineages <- data.frame(Lineage = toupper(letters[1:15]),
                          color = cols_lineages)
df_hosts <- data.frame(host=hosts,
                       color=cols_hosts)

names(cols_countries) <- as.character(df_countries$countryname)
names(cols_lineages) <- as.character(df_lineages$Lineage)
names(cols_hosts) <- as.character(df_hosts$host)

# function to merge dataframes with root-to-tip, lineages and metadata files
prepare_metadata <- function(metadata_file, lineages_file, tempest_file) {
  
  # read metadata
  meta <- read.csv(metadata_file)
  meta$date <- sapply(strsplit(as.character(meta$taxa), "|", fixed = TRUE), `[`, 4)
  meta$accession <- sapply(strsplit(as.character(meta$taxa), "|", fixed = TRUE), `[`, 1)
  meta$date <- as.Date(as.character(meta$date), format = "%Y-%m-%d")
  meta$year <- sapply(strsplit(as.character(meta$date), "-", fixed = TRUE), `[`, 1)
  meta$sample_date <- decimal_date(meta$date)
  
  # read lineages file and merge with metadata
  lineages <- read.csv(lineages_file, header = T)
  m <- left_join(meta,lineages, by=c("accession"="Query"))


  # read tempest data
  tempest_data <- read.csv(tempest_file, header = T, sep = "\t")
  tempest_data$date2 <- sapply(strsplit(as.character(tempest_data$tip), "|", fixed = TRUE), `[`, 4)
  tempest_data$date2 <- as.Date(as.character(tempest_data$date2), format = "%Y-%m-%d")
  tempest_data$accession <- sapply(strsplit(tempest_data$tip, "|", fixed = TRUE), `[`, 1)
  # merge dataframes
  tempest_data <- left_join(tempest_data, m, by=c("accession"))
  meta <- tempest_data %>% dplyr::select(c(taxa, accession, Lineage, country, location, host, date.x, date.y, sample_date, distance, date2, year))
  meta$countryname <- countrycode(meta$country, origin = 'iso3c', destination = 'country.name')
  
  # fix host names
  meta$host[meta$host == "bovine"] <- "cow"
  
  # convert host common names to scientific names
  meta$sciname <- NA
  
  meta <- meta %>%
    mutate(sciname = case_when(
      host == 'cow' ~ 'Bos taurus',
      host == 'human' ~ 'Homo sapiens',
      host == 'mosquito' ~ 'Aedes',
      host == 'sheep' ~ 'Ovis aries',
      host == 'bat' ~ 'Chiroptera',
      host == 'goat' ~ 'Capra',
      host == 'buffalo' ~ 'Bubalina',
      host == 'unknown' ~ 'Aedes',
      host == 'springbok' ~ 'Antidorcas',
      TRUE ~ host)
    )
  
  meta$uid <- NA
  meta <- meta %>%
    mutate(uid = case_when(
      sciname == 'Bos taurus' ~ 'ea92020a-298a-4977-b6d3-dccbe816bb5e',
      sciname == 'Homo sapiens' ~ 'c089caae-43ef-4e4e-bf26-973dd4cb65c5',
      sciname == 'Aedes' ~ '37fd37ec-5ac1-4429-bee9-ce66a3a2eca4',
      sciname == 'Ovis aries' ~ 'a9297cbd-10ca-457b-b5d5-a0b038720df7',
      sciname == 'Chiroptera' ~ '384a1c02-9a04-4688-a1e5-eb15966707e3',
      sciname == 'Capra' ~ '52c28a94-15e2-482d-ac4d-f75ce7b09b95',
      sciname == 'Bubalina' ~ '5ede126f-6a18-4a08-898d-e48828e7dcaa',
      sciname == 'Antidorcas' ~ '320dcfd5-c738-484b-ac28-dab0ac07139b',
      TRUE ~ sciname)
    )
  meta$host <- str_to_title(meta$host)
  print(meta$host)
  return(meta)
}

# function to plot temporal data
plot_root_to_tip <- function(tempest_data, root, pos, r2, n, cor, column, segment, cols, breaks) {
  pos <- as.Date(decimal2Date(pos))
  column <- ensym(column)
  min_y = round(min(tempest_data$distance), 4)
  max_y = round(max(tempest_data$distance), 4)
  min_x = as.Date(min(tempest_data$date2))
  max_x = as.Date(max(tempest_data$date2))
  by = (max_y-min_y)/4
  root <- as.Date(decimal2Date(root))
  
  plot <- ggplot(tempest_data, aes(date2, distance)) +
    geom_point(aes(fill=!!column), size=10, shape=21, stroke=0.3, alpha=0.9) +
    geom_smooth(method=lm, se=T, alpha=0.5, color="#000000", fullrange=TRUE, linewidth=0.5) +
    scale_x_date(date_labels = "%Y", breaks=breaks, expand = c(0,0), 
                 limits = c(root, max_x)) +
    scale_y_continuous(expand=c(0,0), limits=c(0, max_y)) +
    scale_fill_manual(values=cols, name=column) +
    ggtitle("", subtitle = segment) +
    ylab("Root-to-tip Distance") +
    xlab("Date") +
    # annotate(geom="text", x=pos, y=by*1,
    #          label=paste("~italic(R)^{2} == ", r2), parse = TRUE,
    #          color="#000000", size=16/.pt, family="Arial Narrow", hjust = 0, fontface = "bold") +
    annotate(geom="text", x=pos, y=(by*1)+(by/4), label=paste("n = ", n),
             color="#000000", size=16/.pt, family="Arial Narrow", hjust = 0, fontface = "bold") +
    # annotate(geom="text", x=pos, y=by*1.5, label=paste("Age = ", year(as.Date(root, format="%Y"))),
    #          color="#000000", size=16/.pt, family="Arial Narrow", hjust = 0) +
    annotate(geom="text", x=pos, y=(by*1.5)+(by/4), label=paste("R = ", cor),
             color="#000000", size=16/.pt, family="Arial Narrow", hjust = 0, fontface = "bold") +
    labs(x = "date (year)", y = "root-to-tip distance") +
    theme_classic(base_size = 18) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      text = element_text(family = "Arial Narrow"),
      axis.title = element_text(size = 25, color="#000000"),
      axis.text = element_text(size = 20, color="#000000"),
      axis.text.x = element_text(size = 20, color="#000000", angle = 0),
      plot.subtitle=element_text(size = 25, color="#000000", face = "bold"),
      legend.position = c(0.3, 0.8),
      legend.background = element_blank(),
      legend.title = element_text(size = 24, color="#000000"),
      legend.text =  element_text(size = 15, color="#000000")
    ) +
    guides(fill = guide_legend(nrow = 5, override.aes = list(size = 8)))
  return(plot)
}


# plot tree
plot_tree <- function(tree, meta, root, cols_lineage, segment){
  
  #'@param tree tree object
  #'@param meta path to the file containing taxa metadata
  #'@param root decimal date value of the TMRCA of the tree (check with tempest)
  
  # make dataframe for clade nodes
  clades.df <- data.frame(clade=unique(meta$Lineage), node=NA)
  
  # find the most recent common ancestor for each clade
  for (i in 1:length(clades.df$clade)) {
    clades.df$node[i] <- MRCA(
      tree,
      meta$taxa[meta$Lineage == clades.df$clade[i]]
    )
  }
  
  # get mrsd
  mrsd <- meta[rev(order(as.Date(meta$date2, format = "%Y-%m-%d"))),]$date2[1]
  print(mrsd)
  
  # plot preliminary tree
  gg.tree <- ggtree(tree,
                    mrsd = mrsd,
                    as.Date = FALSE) %<+% meta +
    geom_highlight(data = clades.df,
                   aes(node=node, fill=clade),
                   alpha=1,
                   align="right",
                   extend=0.1,
                   show.legend=FALSE) +
    geom_tree(linewidth=0.5) +
    geom_tiplab(aes(label=accession), size=4)
  
  # order clades dataframe to match the tree
  clades.df <- clades.df[match(gg.tree$data %>%
                                 dplyr::filter(isTip == "TRUE") %>%
                                 arrange(y) %>%
                                 pull(Lineage) %>%
                                 unique(),
                               clades.df$clade),]
  
  # add a column with alternating binary value
  clades.df$highlight <- rep(c(0,1), length.out=length(clades.df$clade))
  print(clades.df)
  
  # highlight alternating colors for the clades
  pt1 <- ggtree(tree,
               mrsd = mrsd,
               as.Date = FALSE,
               size = 0.3) %<+% meta +
    geom_highlight(data = clades.df,
                   aes(node=node, fill=as.factor(highlight)),
                   alpha = 1,
                   align = "right",
                   extend = 0.02,
                   show.legend = FALSE
    ) +
    # geom_cladelab(data = clades.df,
    #               mapping = aes(node=node, label=clade, color=clade),
    #               fontsize=5,
    #               family="Arial Narrow",
    #               barsize=4.0,
    #               align=TRUE,
    #               offset=0.04,
    #               offset.text=0.01) +
    labs(x="time (year)") +
    geom_tree(linewidth = 0.5) +
    expand_limits(y = 5) +
    ggtitle("", subtitle = segment) +
    geom_tippoint(aes(color = Lineage), alpha = 1.0, size=8) +
    scale_color_manual(name = "Lineage",
                       values = cols_lineage,
                       breaks = toupper(letters[1:15]),
                       guide = guide_legend(override.aes = list(linetype = 0,
                                                                shape = 16,
                                                                color = cols_lineage))) +
    scale_fill_manual(values = c("#F5F5F5", "#ECECEC")) +
    theme_tree2() +
    theme(text = element_text(family = "Arial Narrow"),
          plot.tag = element_text(size = 30),
          axis.text.x = element_text(size = 20, colour = "#000000"),
          axis.title.x = element_text(size = 25, colour = "#000000"),
          legend.title = element_text(size = 24),
          legend.background = element_blank(),
          legend.position = c(0.2, 0.8),
          legend.key = element_blank(),
          legend.text =  element_text(size = 15, color="#000000"),
          plot.subtitle=element_text(size = 25, color="#000000", face = "bold")) + 
    guides(colour = guide_legend(nrow = 5, override.aes = list(size = 9)))
  pt1
  
  pt2 <- ggtree(tree, 
               layout = "circular", 
               mrsd = mrsd, 
               as.Date = FALSE, 
               size=0.3,
               color = "#000000") %<+% meta +
    geom_highlight(data = clades.df,
                   aes(node=node, fill=as.factor(highlight)),
                   alpha=1,
                   align="right",
                   extend=0.01,
                   show.legend=FALSE) +
    geom_tree(linewidth=0.5) +
    ggtitle("", subtitle = segment) +
    geom_tippoint(aes(color = Lineage), alpha = 1.0, size=13) +
    scale_color_manual(name = "Lineage",
                       values = cols_lineage,
                       breaks = sort(unique(meta$Lineage))) +
    scale_fill_manual(values = c("#F5F5F5", "#ECECEC")) +
    theme(text = element_text(family = "Arial Narrow"),
          plot.tag = element_text(size = 30),
          legend.title = element_text(size = 24),
          legend.background = element_blank(),
          legend.position = c(1.1, 0.4),
          legend.key = element_blank(),
          legend.text =  element_text(size = 15, color="#000000"),
          plot.subtitle=element_text(size = 25, color="#000000", face = "bold")) + 
    guides(colour = guide_legend(nrow = 5, override.aes = list(size = 12)))
  trees_list <- list(pt1,pt2,clades.df)
  return(trees_list)
}


# get snps
snps_data <- function(snps_fn){
  SNPs <- data.table::fread(snps_fn)
  gapChar <- "N"
  SNP <- t(SNPs)
  SNP[1,]
  lSNP <- apply(SNP, 1, function(x) {
    x != SNP[1,] & x != gapChar & SNP[1,] != gapChar
  })
  lSNP <- as.data.frame(lSNP)
  lSNP$pos <- as.numeric(rownames(lSNP))
  lSNP <- tidyr::gather(lSNP, name, value, -pos)
  SNP_data <- lSNP[lSNP$value, c("name", "pos")]
  SNP_data$name <- chartr(".", "|", SNP_data$name)
  return(SNP_data)
}


plot_tree_snps <- function(snps_data, tree, column, cols){
  column <- ensym(column)
  p <- tree + geom_facet(panel = "SNPs", 
                         data = snps_data, geom = geom_point,
                         mapping=aes(x = pos, color = !!column), 
                         size=4.5, shape = "•") +
    xlab("             date (year)                                          position (bp)\n") +
    scale_color_manual(name=column,
                       values=cols) +
    theme(text = element_text(family = "Arial Narrow"),
          axis.text.x = element_text(color="#000000", size=20),
          legend.position = c(0.1, 0.6),
          strip.background = element_rect(colour=NA,fill="#BEBEBE"),
          strip.text.x = element_text(size=16, face="bold"),
          legend.background = element_rect(color = "#000000",
                                           fill = "transparent",
                                           linewidth = 5, linetype = "blank"))
  return(p)
}
############################ Medium segment #####################################

# metadata_fn <- file.path(wd, "M-global.filtered.txt")
metadata_fn <- file.path(wd, "M-global.geolocations.csv")
lineages_fn <- file.path(wd, "M-global.lineages.csv")
tempest_fn <- file.path(wd, "M-global.root-to-tip.tsv")
timetree_fn <- file.path(wd, "M-global.timetree.nexus")

meta.medium <- prepare_metadata(metadata_file = metadata_fn,
                                lineages_file = lineages_fn,
                                tempest_file = tempest_fn)

p1 <- plot_root_to_tip(tempest_data = meta.medium,
                       root = 1899.84,
                       r2 = 0.63,
                       n = 237,
                       cor = 0.795,
                       pos = 1910,
                       column = "Lineage",
                       segment = "RVFV-M temporal signal",
                       cols = cols_lineages,
                       breaks = '15 year')
p1

# read global M segment nexus tree from TreeTime
tree_treetime <- read.nexus(timetree_fn)
meta.medium <- meta.medium[match(tree_treetime$tip.label, meta.medium$taxa), ]
all(tree_treetime$tip.label == meta.medium$taxa)

trees.m <- plot_tree(tree=tree_treetime,
                  meta = meta.medium,
                  root = 1899.8452,
                  cols_lineage = cols_lineages,
                  segment = "RVFV-M maximum likelihood tree")
trees.m[[1]]
trees.m[[2]]

# read global M segment MCC tree from BEAST
tree_medium_mcc <- "M-global.mcc.tree"
tree_medium_mcc <- read.beast(tree_medium_mcc)
meta.medium <- meta.medium[match(tree_medium_mcc@phylo$tip.label, 
                                 meta.medium$taxa),]
all(tree_medium_mcc@phylo$tip.label == meta.medium$taxa)
trees.mcc <- plot_tree(tree=tree_medium_mcc,
                       meta = meta.medium,
                       root = 1918.346,
                       cols_lineage = cols_lineages,
                       segment = "RVFV-M maximum clade credibility tree")
trees.mcc[[1]]

# read and process snps data
snps_fn <- "M-global.alignment.snps.csv"
m_snp_data <- snps_data(snps_fn = snps_fn)

m.global.mcc.snps.plot <- plot_tree_snps(snps_data = m_snp_data,
                                         tree = trees.mcc[[1]],
                                         column = "Lineage",
                                         cols = cols_lineages)
m.global.mcc.snps.plot <- m.global.mcc.snps.plot + 
  guides(colour = guide_legend(nrow = 15, override.aes = list(size = 12)))

m.global.mcc.snps.plot


# read log file
logfile <- file.path(wd, "M-global.log")
log_df <- read.csv(logfile, skip = 4, header = T, sep = "\t")
age.root <- log_df %>% dplyr::select(c("age.root.")) %>% dplyr::rename(Age="age.root.")
age.root$Lineage <- "RVFV"
age.linA <- log_df %>% dplyr::select(c("age.LineageA.")) %>% dplyr::rename(Age="age.LineageA.")
age.linA$Lineage <-"A"
age.linC <- log_df %>% dplyr::select(c("age.LineageC.")) %>% dplyr::rename(Age="age.LineageC.")
age.linC$Lineage <- "C"
age.linH <- log_df %>% dplyr::select(c("age.LineageH.")) %>% dplyr::rename(Age="age.LineageH.")
age.linH$Lineage <- "H"

d <- bind_rows(age.root, age.linA, age.linC, age.linH)

custom <- c("#4d4d4d", "#9370DB", "#FF0000", "#E88817")
names(custom) <- unique(d$Lineage)
colors <- data.frame("Lineage"=c(unique(d$Lineage)), "colour"=custom)
d <- left_join(d, colors)

d$Age <- date_decimal(d$Age)
d$Age <- as.Date(cut(d$Age, breaks = "days", 
                     start.on.monday = FALSE),
                 format = "%Y-%m-%d")

mean.ages <- data.frame(Age=c(mean(age.root$Age), 
                              mean(age.linA$Age), 
                              mean(age.linC$Age),
                              mean(age.linH$Age)),
                        Lineage=c("RVFV", "A", "C", "H"))
mean.ages$Age <- date_decimal(mean.ages$Age)
mean.ages$Age <- as.Date(cut(mean.ages$Age, breaks = "days", 
                             start.on.monday = FALSE),
                         format = "%Y-%m-%d")

names(custom) <- unique(mean.ages$Lineage)
colors <- data.frame("Lineage"=c(unique(mean.ages$Lineage)), "colour"=custom)
mean.ages <- left_join(mean.ages, colors)

lineages.age.plot <- ggplot(d, aes(x=Age, fill = Lineage)) +
  scale_x_date(date_labels = "%Y", breaks='15 year',
               limits = c(as.Date(as.character(1850), "%Y"),
                          as.Date(as.character(2022), "%Y"))) +
  scale_fill_manual(name="Lineage", values=custom) +
  geom_density(alpha=1.0, show.legend = T) +
  geom_vline(data = mean.ages[mean.ages$Lineage =="RVFV",], 
             aes(xintercept = Age), colour = custom[1], 
                 linewidth=0.3, show.legend = FALSE) +
  geom_vline(data = mean.ages[mean.ages$Lineage =="A",], 
             aes(xintercept = Age), colour = custom[2], 
             linewidth=0.3, show.legend = FALSE) +
  geom_vline(data = mean.ages[mean.ages$Lineage =="C",], 
             aes(xintercept = Age), colour = custom[3], 
             linewidth=0.3, show.legend = FALSE) +
  geom_vline(data = mean.ages[mean.ages$Lineage =="H",], 
             aes(xintercept = Age), colour = custom[4], 
             linewidth=0.3, show.legend = FALSE) +
  annotate(geom="text", x=as.Date(as.character(1918), "%Y"), y=0.0015, 
           label="1918", hjust=0, fontface="bold",
           color="#4d4d4d", size=8) +
  annotate(geom="text", x=as.Date(as.character(1971), "%Y"), y=0.0019, 
           label="1971",hjust=0, fontface="bold",
           color="#9370DB", size=8) +
  annotate(geom="text", x=as.Date(as.character(1970), "%Y"), y=0.0017, 
           label="1970",hjust=0, fontface="bold",
           color="#FF0000", size=8) +
  annotate(geom="text", x=as.Date(as.character(2007), "%Y"), y=0.0021, 
           label="2007",hjust=0, fontface="bold",
           color="#E88817", size=8) +
  labs(y="node age density", x = "age (year)") +
  ggtitle("", subtitle = "TMRCA estimates") +
  theme_classic(base_size = 18) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title = element_text(size=25, color="#000000"),
    axis.text = element_text(size=20, color="#000000"),
    plot.subtitle=element_text(size=25, color="#000000", face = "bold"),
    legend.position = c(0.15, 0.4),
    legend.background = element_blank(),
    legend.title = element_text(size = 24, color="#000000"),
    legend.text =  element_text(size = 15, color="#000000")
  ) +
  guides(fill = guide_legend(nrow = 4, override.aes = list(size = 9)))
lineages.age.plot

# . plot skygrid population history of the virus
plot_pop <- function(file, age, lower_tmrca, upper_tmrca, youngest_date, 
                     subtitle, col, segment) {
  
  # Read model output from BEAST
  df <- read.csv(file, skip = 1, header = T, sep = "\t")
  df <- df %>% dplyr::rename(time=1,
                             mean=2,
                             median=3,
                             upper=4,
                             lower=5)
  base_n <- dim(df)[1]
  BaseDateVector <- rep(youngest_date, base_n) # base number from BEAST
  print(BaseDateVector)
  decimalDates <- (BaseDateVector - df$time)
  calendarDates = format(date_decimal(df$time), "%Y") # convert to calendar dates
  
  # add calender and decimal dates difference to dataframe
  dfCalendar = data.frame(df, calendarDates)
  dfCalendar <- data.frame(dfCalendar, decimalDates)
  
  # convert to class Date
  dfCalendar$calendarDates <- as.Date(dfCalendar$calendarDates, "%Y")
  
  # define TMRCA from BEAST tracer summary
  tmrca <- as.Date(age, "%Y")
  minTMRCA <- as.Date(lower_tmrca, "%Y")
  maxTMRCA <- as.Date(upper_tmrca, "%Y")
  
  tmrca_year <- as.numeric(format(as.POSIXct(tmrca, format = "%Y-%m-%d %H:%M:%S"), format = "%Y"))
  min_tmrca_year <- as.numeric(format(as.POSIXct(minTMRCA, format = "%Y-%m-%d %H:%M:%S"), format = "%Y"))
  max_tmrca_year <- as.numeric(format(as.POSIXct(maxTMRCA, format = "%Y-%m-%d %H:%M:%S"), format = "%Y"))
  
  
  # define range of median values
  min_lower <- min(dfCalendar$lower)
  max_upper <- max(dfCalendar$upper)
  
  # plot
  
  # adjust these values according to your data
  n = length(seq(from=log(min_lower), to=log(max_upper)))
  by = 2
  breaks = seq(from=log(min_lower), to=log(max_upper), by=by)
  print(breaks)
  
  if (segment == 'L'){
    labels = c(0.1, 1, 10, 100, 1000, 10000)
  }
  if (segment == 'M'){
    labels = c(0.1, 1, 10, 100, 1000)
  }
  if (segment == 'nss'){
    labels = c(0.1, 1, 10)
  }
  
  
  p <- ggplot(data = dfCalendar, mapping = aes(x=calendarDates)) +
    geom_line(aes(y = log(mean)), colour = col , alpha = 0.6) +
    geom_line(aes(y = log(median)), colour = col, alpha = 0.8, linetype="dashed") +
    geom_ribbon(aes(ymin = log(lower), ymax = log(upper)), fill = col, alpha = 0.2) +
    scale_x_date(date_labels = "%Y", date_breaks = "10 year") +
    scale_y_continuous(breaks = seq(from=log(min_lower), to=log(max_upper), by=by),
                       labels = labels) +
    labs(y=expression(paste(italic("N"[e]),tau," ","(years)")), x = "year") +
    ggtitle("", subtitle = subtitle) +
    theme_classic(base_size = 18) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      text = element_text(family = "Arial Narrow"),
      axis.title = element_text(size=25, color="#000000"),
      axis.text = element_text(size=20, color="#000000"),
      plot.subtitle=element_text(size=25, color="#000000", face = "bold")
    )
}

pop.M.global <- plot_pop(file = file.path(wd, "M-global.skygrid.tsv"),
                  age = "1918.246",
                  lower_tmrca = "1902.8385",
                  upper_tmrca = "1930.7063",
                  youngest_date = 2022.695,
                  col = "#FFA500",
                  subtitle = "RVFV demographic history",
                  segment = "M")
pop.M.global


country.codes <- c("AGO", "BDI", "BEN", "BFA", "BWA", "CAF", "CIV",
                   "CMR", "COG", "CPV", "DJI", "DZA", "EGY", "ERI", "ESH",
                   "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GNQ", "KEN",
                   "LBR", "LBY", "LSO", "MAR", "MDG", "MLI", "MOZ", "MRT",
                   "MUS", "MWI", "MYT", "NAM", "NER", "NGA", "RWA", "SDN",
                   "SEN", "SLE", "SOM", "SSD", "SWZ", "SYC", "TCD", "TGO",
                   "TUN", "TZA", "UGA", "ZAF", "ZMB", "ZWE", "SAU", "YEM",
                   "OMN")


country.names <- countrycode(country.codes,
                             origin = 'iso3c',
                             destination = 'country.name')
country.names
country.names <- append(country.names, c("Somaliland", "Ivory Coast", "Algeria",
                                         "Cameroon", "Cape Verde", 
                                         "Central African Republic",
                                         "Democratic Republic of the Congo",
                                         "Equatorial Guinea", "Gabon", "Gambia",
                                         "Guinea-Bissau", "Morocco"
))
country.names <- unique(country.names)

data <- meta.medium
per.country.count <- meta.medium %>% dplyr::count(countryname) %>% dplyr::rename(country="countryname", count="n")
per.country.lineage.count <- meta.medium %>% dplyr::count(countryname,Lineage) %>% dplyr::rename(country="countryname", count="n")
per.country.lineage.count <- per.country.lineage.count[order(per.country.lineage.count$Lineage),]


df.wide <- per.country.lineage.count %>%
  mutate(names = ifelse(str_detect(Lineage, "[A-Z]"), Lineage, NA)) %>%
  fill(names) %>%
  filter(!is.na(count)) %>%
  pivot_wider(names_from = Lineage,
              values_from = count)
df.wide[is.na(df.wide)] <- 0
wld <- map_data('world')
coordinates <- wld %>% dplyr::group_by(region) %>% dplyr::summarise(long=mean(long), lat=mean(lat))
df.lineages.per.country <- left_join(df.wide, coordinates, by=c("country"="region"))
df.lineage.totals <- data.frame(Lineage=LETTERS[1:15], Total=colSums(df.lineages.per.country[ , c(LETTERS[1:15])], na.rm=TRUE))
rownames(df.lineage.totals) <- NULL
df.lineages.per.country <- left_join(df.lineages.per.country, df.lineage.totals, by=c("names"="Lineage"))
df.lineages.per.country$radius <- (rowSums(df.lineages.per.country[ , c(LETTERS[1:15])], na.rm=TRUE)/df.lineages.per.country$Total)
df.lineages.per.country$radius
df.lineages.per.country$radius <- log10(df.lineages.per.country$radius)+3
df.lineages.per.country$radius


spdf <- spData::world %>% 
  dplyr::filter((continent == 'Africa' ) | 
                  (continent == 'Asia' & name_long %in% 
                     c('Saudi Arabia', 'Yemen', 'Oman')))
global <- dplyr::left_join(spdf, per.country.count, 
                           by = c("name_long" = "country")) %>% 
  dplyr::rename("name"="name_long")

scatterpie.plot.m <- ggplot(spdf) +
  geom_sf(aes(geometry = geom), show.legend = "point") +
  coord_sf(default_crs = 4326, 
           xlim = c(-17.520278, 63.47519), 
           ylim = c(-34.65302, 37.34698)) +
  geom_scatterpie(aes(x=long, y=lat, group=country, r=radius),
                  data=df.lineages.per.country, 
                  cols=LETTERS[1:15], color=NA, alpha=.8) +
  geom_scatterpie_legend(df.lineages.per.country$radius, n=3, 
                         x=60, y=37) +
  ggtitle("", subtitle = "RVFV lineage occurrence") +
  theme_map() +
  scale_fill_manual(name="Lineage",
                    values = df_lineages$color,
                    breaks = df_lineages$Lineage,
                    guide = guide_legend(override.aes = list(linetype = 0,
                                                             shape = 16,
                                                             color = df_lineages$color))) +
  theme(
    plot.tag = element_text(size = 30, face = "bold"),
    text = element_text(family = "Arial Narrow"),
    legend.title = element_text(size = 20, colour = "#000000"),
    legend.background = element_blank(),
    legend.position = c(-0.02, 0.1),
    legend.text = element_text(size = 16, colour = "#000000"),
    plot.subtitle = element_text(color = "#000000", size = 25, face = "bold")) +
  guides(fill = guide_legend(nrow = 5,
                             override.aes = list(size = 9, shape=1),
                             title = "Lineage"))
scatterpie.plot.m


genomes.plot <- ggplot(global) +
  geom_sf(aes(geometry = geom, fill = count)) +
  ggrepel::geom_label_repel(
    data = subset(global, count > 0),
    aes(label = name, geometry = geom, size = 5.5),
    stat = "sf_coordinates",
    min.segment.length = 0,
    colour = "#000000",
    fontface = "bold",
    family = "Arial Narrow",
    segment.colour = "#FFFFFF",
    show.legend = F) +
  scale_fill_viridis_c(alpha = 0.8, option="plasma", na.value = "#FFFFFF") +
  geom_sf_label(aes(label = count), alpha=0.8, size=5, 
                label.padding = unit(0.1, "lines")) +
  ggtitle("", subtitle = "RVFV genomes") +
  coord_sf(datum = NA)+
  theme_map()+
  theme(
    plot.tag = element_text(size = 30, face = "bold"),
    text = element_text(family = "Arial Narrow"),
    legend.title = element_text(size = 20),
    legend.background = element_blank(),
    legend.position = c(-0.01, 0.18),
    legend.text = element_text(size = 15),
    plot.subtitle = element_text(size = 25, color = "#000000", face = "bold")) + 
  guides(fill = guide_legend(title = "Genome sequences", 
                             title.position = "top", 
                             title.theme = element_text(size = 20, 
                                                        colour = "#000000", 
                                                        angle = 0))
  )

genomes.plot


figure.1.plots <- ((((genomes.plot + theme(axis.title = element_blank(),
                                           axis.text = element_blank(),
                                           legend.text = element_text(size = 24, colour = "#000000"))) / (scatterpie.plot.m + theme(axis.title = element_blank(),
                                                                                                      axis.text.y = element_blank(),
                                                                                                      axis.text.x = element_text(size=25, colour = "#000000", face="bold", angle=0),
                                                                                                      legend.position = "none"))) |
  (p1 + theme(axis.text.y = element_text(size = 20, colour = "#000000", face = "bold"),
              axis.text.x = element_text(size = 20, colour = "#000000", face = "bold", angle = 90),
              legend.position = "none",
              axis.title.x = element_blank()) +
     guides(fill = guide_legend(nrow = 5,
                                override.aes = list(size = 8),
                                title = "Lineage"))) / 
    (lineages.age.plot + 
       theme(axis.text.y = element_text(size = 20, colour = "#000000", face = "bold"),
             axis.text.x = element_text(size = 20, colour = "#000000", face = "bold", angle = 90),
             axis.title.x = element_blank(),
             legend.position = "none")) /
    (pop.M.global + 
       theme(axis.text.y = element_text(size = 20, colour = "#000000", face = "bold"),
             axis.text.x = element_text(size = 20, colour = "#000000", face = "bold", angle = 90)
       ))) | 
    (m.global.mcc.snps.plot + theme(legend.position = "right",
                                    axis.text.x = element_text(size = 20, colour = "#000000", face = "bold", angle = 90)
    ))) + 
  plot_layout(guides = "keep",
              widths = c(100, 100, 120),
              heights = c(40, 75, 75),
              byrow = T) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 30, colour = "#000000", face = "bold"),
        plot.subtitle = element_text(size = 26, colour = "#000000", face = "bold"),
        strip.text.x = element_text(size=25, colour = "#000000", face="bold"),
        legend.title = element_text(size = 24, colour = "#000000", face = "bold"))


ggsave(file.path(figuresDir, "Figure1.pdf"),
       figure.1.plots,
       width = 28, height = 18, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDir, "Figure1.png"),
       figure.1.plots,
       width = 28, height = 18, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = "png")

ggsave(file.path(figuresDir, "Figure1.tiff"), 
       figure.1.plots,
       width = 28, height = 18, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = "tiff")





############################ Large segment #####################################

metadata_fn <- file.path(wd, "L-global.geolocations.csv")
lineages_fn <- file.path(wd, "L-global.lineages.csv")
tempest_fn <- file.path(wd, "L-global.root-to-tip.tsv")
timetree_fn <- file.path(wd, "L-global.timetree.nexus")

meta.large <- prepare_metadata(metadata_file = metadata_fn,
                               lineages_file = lineages_fn,
                               tempest_file = tempest_fn)

# fix wrongly assigned sequences
meta.large$Lineage[meta.large$accession == 'DQ375434'] <- 'J' # wrongly assigned as I
meta.large$Lineage[meta.large$accession == 'MG659874'] <- 'H' # wrongly assigned as C
meta.large$Lineage[meta.large$accession == 'MG659825'] <- 'H' # wrongly assigned as C
meta.large$Lineage[meta.large$accession == 'MG659886'] <- 'H' # wrongly assigned as C
meta.large$Lineage[meta.large$accession == 'KY126678'] <- 'H' # wrongly assigned as C
meta.large$Lineage[meta.large$accession == 'DQ375425'] <- 'D' # wrongly assigned as C
meta.large$Lineage[meta.large$accession == 'DQ375433'] <- 'O' # wrongly assigned as L


p3 <- plot_root_to_tip(tempest_data = meta.large,
                       root = 1896.25,
                       r2 = 0.94,
                       n = 236,
                       cor = 0.97,
                       pos = 1905,
                       column = "Lineage",
                       segment = "RVFV-L temporal signal",
                       cols = cols_lineages,
                       breaks = '14 year')
p3

# read nexus tree
tree_treetime <- read.nexus(timetree_fn)
meta.large <- meta.large[match(tree_treetime$tip.label, meta.large$taxa), ]
all(tree_treetime$tip.label == meta.large$taxa)

trees.l <- plot_tree(tree=tree_treetime,
                meta = meta.large,
                root = 1896.25,
                cols_lineage = cols_lineages,
                segment = "RVFV-L maximum likelihood tree")

# read and process snps data
snps_fn <- file.path(wd, "L-global.alignment.snps.csv")
l_snp_data <- snps_data(snps_fn = snps_fn)
l.global.treetime.snps.plot <- plot_tree_snps(snps_data = l_snp_data,
                                              tree = trees.l[[1]],
                                              column = "Lineage",
                                              cols = cols_lineages)
l.global.treetime.snps.plot <- l.global.treetime.snps.plot + 
  guides(colour = guide_legend(nrow = 15, override.aes = list(size = 12)))
l.global.treetime.snps.plot

ggsave("~/projects/GenPath_Africa/progress/plots/202405-Figure1.pdf",
       l.global.treetime.snps.plot,
       width = 20, height = 15, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave("~/projects/GenPath_Africa/progress/plots/202405-Figure1.png",
       l.global.treetime.snps.plot,
       width = 20, height = 15, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white")

########################## Small segment - NP ##################################

metadata_fn <- file.path(wd, "S-global.geolocations.csv")
lineages_fn <- file.path(wd, "S-global.lineages.csv")
tempest_fn <- file.path(wd, "np-global.root-to-tip.tsv")
timetree_fn <- file.path(wd, "np-global.timetree.nexus")

meta.small <- prepare_metadata(metadata_file = metadata_fn,
                               lineages_file = lineages_fn,
                               tempest_file = tempest_fn)

# fix wrongly assigned sequences
meta.small$Lineage[meta.small$accession == 'KY196500'] <- 'H' # wrongly assigned as C
meta.small$Lineage[meta.small$accession == 'DQ380167'] <- 'G' # wrongly assigned as E
meta.small$Lineage[meta.small$accession == 'DQ380165'] <- 'G' # wrongly assigned as E
meta.small$Lineage[meta.small$accession == 'DQ380168'] <- 'G' # wrongly assigned as E
meta.small$Lineage[meta.small$accession == 'DQ380161'] <- 'G' # wrongly assigned as E
meta.small$Lineage[meta.small$accession == 'DQ380160'] <- 'G' # wrongly assigned as E

p5 <- plot_root_to_tip(tempest_data = meta.small,
                       root = 1903.33,
                       r2 = 0.72,
                       n = 247,
                       cor = 0.85,
                       pos = 1912,
                       column = "Lineage",
                       segment = "RVFV-NP temporal signal",
                       cols = cols_lineages,
                       breaks = '17 year')
p5

# read nexus tree
tree_treetime <- read.nexus(timetree_fn)
meta.small <- meta.small[match(tree_treetime$tip.label, meta.small$taxa), ]
all(tree_treetime$tip.label == meta.small$taxa)

trees.s.np <- plot_tree(tree=tree_treetime,
                        meta = meta.small,
                        root = 1903.3254,
                        cols_lineage = cols_lineages,
                        segment = "RVFV-S-NP maximum likelihood tree")

######################### Small segment - NSS ##################################

tempest_fn <- file.path(wd, "nss-global.root-to-tip.tsv")
timetree_fn <- file.path(wd, "nss-global.timetree.nexus")

meta.small <- prepare_metadata(metadata_file = metadata_fn,
                               lineages_file = lineages_fn,
                               tempest_file = tempest_fn)

# fix wrongly assigned sequences
meta.small$Lineage[meta.small$accession == 'KY196500'] <- 'H' # wrongly assigned as C
meta.small$Lineage[meta.small$accession == 'DQ380167'] <- 'G' # wrongly assigned as E
meta.small$Lineage[meta.small$accession == 'DQ380165'] <- 'G' # wrongly assigned as E
meta.small$Lineage[meta.small$accession == 'DQ380168'] <- 'G' # wrongly assigned as E
meta.small$Lineage[meta.small$accession == 'DQ380161'] <- 'G' # wrongly assigned as E
meta.small$Lineage[meta.small$accession == 'DQ380160'] <- 'G' # wrongly assigned as E

p7 <- plot_root_to_tip(tempest_data = meta.small,
                       root = 1874.13,
                       r2 = 0.35,
                       n = 247,
                       cor = 0.59,
                       pos = 1890,
                       column = "Lineage",
                       segment = "RVFV-NSS temporal signal",
                       cols = cols_lineages,
                       breaks = '21 year')
p7

# read nexus tree
tree_treetime <- read.nexus(timetree_fn)
meta.small <- meta.small[match(tree_treetime$tip.label, meta.small$taxa), ]
all(tree_treetime$tip.label == meta.small$taxa)

trees.s.nss <- plot_tree(tree = tree_treetime,
                         meta = meta.small,
                         root = 1874.1335,
                         cols_lineage = cols_lineages,
                         segment = "RVFV-S-NSS maximum likelihood tree")

t1 <- (((trees.l[[1]] + theme(legend.position = "none", 
                              axis.text.x = element_text(angle = 90, colour = "#000000", face = "bold"),
                              plot.subtitle = element_blank())) + 
         inset_element((p3 + theme(legend.position = "none",
                                   axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90, colour = "#000000", face = "bold"),
                                   axis.title.y = element_text(size = 22, colour = "#000000"),
                                   axis.text.y = element_text(size = 20, colour = "#000000"))), 
                       left = 0.02,
                       bottom = 0.6,
                       right = 0.50,
                       top = 1.0))) +
  plot_layout(tag_level = "new") +
  labs(tag = "A")
t1

t2 <- (((trees.m[[1]] + theme(plot.subtitle = element_blank(),
                              axis.text.x = element_text(angle = 90, colour = "#000000", face = "bold"),
                              legend.position = "right") +
           guides(colour = guide_legend(nrow = 15,
                                      override.aes = list(size = 20),
                                      title = "Lineage"))) + 
         inset_element((p1 + theme(legend.position = "none",
                                   axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90, colour = "#000000", face = "bold"),
                                   axis.title.y = element_text(size = 22, colour = "#000000"),
                                   axis.text.y = element_text(size = 20, colour = "#000000"))), 
                       left = 0.02,
                       bottom = 0.6,
                       right = 0.50,
                       top = 1.0))) +
  plot_layout(tag_level = "new") +
  labs(tag = "B")
t2

t3 <- (((trees.s.np[[1]] + theme(legend.position = "none",
                                 axis.text.x = element_text(angle = 90, colour = "#000000", face = "bold"),
                                 plot.subtitle = element_blank())) + 
         inset_element((p5 + theme(legend.position = "none",
                                   axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90, colour = "#000000", face = "bold"),
                                   axis.title.y = element_text(size = 22, colour = "#000000"),
                                   axis.text.y = element_text(size = 20, colour = "#000000"))), 
                       left = 0.02,
                       bottom = 0.6,
                       right = 0.50,
                       top = 1.0))) +
  plot_layout(tag_level = "new") +
  labs(tag = "C")

t4 <- (((trees.s.nss[[1]] + theme(legend.position = "none",
                                  axis.text.x = element_text(angle = 90, colour = "#000000", face = "bold"),
                                  plot.subtitle = element_blank())) + 
         inset_element((p7 + theme(legend.position = "none",
                                   axis.title.x = element_blank(),
                                   axis.text.x = element_text(angle = 90, colour = "#000000", face = "bold"),
                                   axis.title.y = element_text(size = 22, colour = "#000000"),
                                   axis.text.y = element_text(size = 20, colour = "#000000"))), 
                       left = 0.02,
                       bottom = 0.6,
                       right = 0.48,
                       top = 1.0))) +
  plot_layout(tag_level = "new") +
  labs(tag = "D")


sf1 <- (((t1) / 
    (t2)) |
  ((t3 + theme(axis.title.y = element_blank())) / 
     (t4 + theme(axis.title.y = element_blank())))) +
  plot_layout(guides = "collect", 
              widths = c(20, 20), 
              heights = c(16, 16)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 34, colour = "#000000", face = "bold"),
        axis.title = element_text(size = 28, colour = "#000000", face = "bold"),
        axis.text.x = element_text(size = 26, colour = "#000000", face = "bold"),
        legend.title = element_text(size = 25, colour = "#000000"),
        legend.text = element_text(size = 24, colour = "#000000"))
  
ggsave(file.path(figuresDir, "FigureS1.pdf"),
       sf1,
       width = 30, height = 28, units = "in",
       limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDir, "FigureS1.png"),
       sf1,
       width = 30, height = 24, units = "in",
       limitsize = FALSE,
       dpi = 300, bg="white", device = "png")

ggsave(file.path(figuresDir, "FigureS1.tiff"),
       sf1,
       width = 30, height = 24, units = "in",
       limitsize = FALSE,
       dpi = 300, bg="white", device = "tiff")


########################### Lineage C mcc trees ################################

lineageC.cols <- c("#FF0000", "#228B22", "#6495ED", "#CC7B00", "#EE82EE")
plot_mcc_lineage_tree <- function(tree, metadata, cols, subtitle, column){
  column <- ensym(column)
  mrsd <- metadata[rev(order(as.Date(metadata$date2, format = "%Y-%m-%d"))),]$date2[1]
  
  p <- ggtree(tree,
         mrsd = mrsd, 
         as.Date = FALSE) %<+% metadata + 
    geom_tippoint(aes(color = !!column), alpha=0.8, size=7.5) +
    geom_nodelab(aes(x=branch, label=round(posterior, 2),
                     subset=(posterior > 0.9)), vjust=-.5, size=6, fontface="bold", family="Arial Narrow") +
    geom_range("height_0.95_HPD", color="#7F7F7F", size=2, alpha=.5, branch.length = "height_median") + 
    # geom_text2(aes(subset=!isTip, label=round(height, digits = 2), hjust=-.3, vjust=1.2), size=3) + 
    scale_color_manual(name=column, values=cols) +
    geom_tree(linewidth = 0.5) +
    labs(subtitle = subtitle) +
    theme_tree2() +
    theme(text = element_text(family = "Arial Narrow"),
          plot.tag = element_text(size = 30),
          axis.text.x = element_text(size = 25, colour = "#000000", face = "bold"),
          axis.title.x = element_text(size = 26, colour = "#000000"),
          legend.title = element_text(size = 25, colour = "#000000"),
          legend.background = element_blank(),
          legend.position = c(0.2, 0.8),
          legend.key = element_blank(),
          legend.text =  element_text(size = 24, colour = "#000000"),
          legend.key.height= unit(0.4, 'cm'),
          legend.key.width= unit(0.2, 'cm'),
          plot.subtitle=element_text(size = 24, colour = "#000000", face = "bold")) + 
    guides(colour = guide_legend(title = str_to_title(column),
                                 nrow = 5, 
                                 override.aes = list(size = 10)))
  return(p)
}


plot_mcc_with_locations <- function(tree, metadata, subtitle){
  startD <- metadata[order(as.Date(metadata$date2, format = "%Y-%m-%d")),]$date2[1]
  mrsd <- metadata[rev(order(as.Date(metadata$date2, format = "%Y-%m-%d"))),]$date2[1]
  subtract_from <- as.numeric(round(Date2decimal(mrsd), 2))
  print(subtract_from)
  tree.plot <- ggtree(tree, mrsd = mrsd, as.Date = F) %<+% 
    metadata +
    geom_tippoint(aes(color = countryname), alpha = 0.8, size=7.5) +
    scale_x_continuous(breaks=seq(year(startD)-10, year(mrsd), 10), minor_breaks=seq(year(startD)-10, year(mrsd), 1)) +
    labs(subtitle = subtitle) +
    geom_range("height_0.95_HPD", color="#7F7F7F", size=2, alpha=.5, branch.length = "height_median") + 
    # geom_text2(aes(subset=!isTip, label=round(height, digits = 2), hjust=-.3, vjust=1.2), size=2, fontface="bold") + 
    # geom_text2(aes(subset=!isTip, 
    #                label=lapply(list(lapply(height_0.95_HPD, round, 2))[[1]], function(x, sub_amnt) sub_amnt - x, sub_amnt = subtract_from),
    #                hjust=1.1, vjust=0), size=2.0, fontface="bold") + 
    geom_nodelab(aes(label=round(as.numeric(posterior), 2), 
                   subset=as.numeric(posterior)> 0.9, 
                   x=branch), vjust=-1, size=6, fontface="bold", family="Arial Narrow") +
    geom_tiplab(aes(label = location), offset=0.8, size = 3.0, align=T, linetype = "dashed", linesize = .3, fontface="bold", family="Arial Narrow") +
    scale_color_manual(name = "country",
                       values = df_countries[df_countries$countryname %in% sort(unique(metadata$countryname)),]$color,
                       breaks = sort(unique(metadata$countryname))) +
    theme_tree2() +
    theme(text = element_text(family = "Arial Narrow"),
          panel.grid.major   = element_line(color="black", size=.2),
          panel.grid.minor   = element_line(color="grey", size=.2),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(size = 26, color = "#000000", face = "bold"),
          legend.text = element_text(size = 25, color = "#000000"),
          legend.title = element_text(size = 26, color = "#000000"),
          plot.subtitle = element_text(size = 30, color = "#000000", face = "bold"),
          legend.position = c(0.2, 0.45),
          legend.background = element_blank()) +
    guides(color = guide_legend(nrow = length(unique(metadata$countryname)), 
                                override.aes = list(size = 15), 
                                title = "country"))
  print(tree.plot)
} 


################################# MCC segment M ################################
metadata_fn <- file.path(wd, "M-global.geolocations.csv")
lineages_fn <- file.path(wd, "M-global.lineages.csv")
tempest_fn <- file.path(wd, "M-lineage-C.root-to-tip.tsv")
mcc.lineageC.medium <- file.path(wd, "M-lineage-C.mcc.tree")

meta.medium.lineageC <- prepare_metadata(metadata_file = metadata_fn,
                                         lineages_file = lineages_fn,
                                         tempest_file = tempest_fn)

p9 <- plot_root_to_tip(tempest_data = meta.medium.lineageC,
                       root = 1968,
                       r2 = 0.91,
                       n = 115,
                       cor = 0.96,
                       pos = 1972,
                       column = "countryname",
                       segment = "RVFV-M Lineage C temporal signal",
                       cols = cols_countries, 
                       breaks = '9 year')
p9 <- p9 + theme(legend.position = c(0.3, 0.78),
                 axis.text.x = element_text(size = 25, colour = "#000000", face = "bold", angle = 90),
                 axis.text.y = element_text(size = 25, colour = "#000000", face = "bold")) +
  guides(fill = guide_legend(nrow = 6,
                             override.aes = list(size = 9),
                             title = "Country"))
p9
# read tree
mcc.tre.lineageC.medium <- read.beast(mcc.lineageC.medium)

# merge with old lineages classification
old_lineages <- read.csv(file.path(wd, "complete-M-Gn-classifier.csv"))
old_lineages <- old_lineages %>% dplyr::rename("Old_Lineage"="Previous_Lineage.Bird.et.al...2008..Bird.et.al..2007.")
merged_lineages <- left_join(meta.medium.lineageC, old_lineages, by=c("accession"="Accession"))

merged_lineages <- merged_lineages %>% dplyr::mutate(Old_Lineage = ifelse(Old_Lineage %in% "", Lineage.x, Old_Lineage))
merged_lineages <- merged_lineages %>% dplyr::mutate(Old_Lineage = ifelse(is.na(Old_Lineage), Lineage.x, Old_Lineage))
merged_lineages <- merged_lineages %>% dplyr::mutate(Strain = ifelse(is.na(Strain), accession, Strain))

meta.medium.lineageC <- meta.medium.lineageC[match(mcc.tre.lineageC.medium@phylo$tip.label, meta.medium.lineageC$taxa), ]
all(mcc.tre.lineageC.medium@phylo$tip.label == meta.medium.lineageC$taxa)

# cols <- cols_countries[names(cols_countries) %in% sort(unique(meta.medium.lineageC$countryname))]
mrsd <- meta.medium.lineageC[rev(order(as.Date(meta.medium.lineageC$date2, format = "%Y-%m-%d"))),]$date2[1]

p10 <- ggtree(mcc.tre.lineageC.medium,
            mrsd = mrsd,
            as.Date = TRUE,
            size = 0.3) %<+% merged_lineages +
  geom_nodelab(aes(label=node), vjust=-.5, size=2.5) +
  geom_nodelab(aes(x=branch, label=round(posterior,2), subset=(posterior > 0.9)), vjust=-.5, size=6) +
  geom_cladelabel(node=179, label="C.1.1", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.3, fontsize = 6.0) +
  geom_cladelabel(node=188, label="C.1.2", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.3, fontsize = 6.0) +
  geom_cladelabel(node=169, label="C.2.1", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.3, fontsize = 6.0) +
  geom_cladelabel(node=125, label="C.2.2", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.08, fontsize = 6.0) +
  # geom_tiplab(aes(label = paste0("italic('", accession, "|", country, "|", year, "')"), offset=0.06), parse = TRUE, size = 6.0, align=F) +
  geom_tippoint(aes(color = Old_Lineage), alpha = 1.0, size=1) +
  scale_x_date(date_breaks= "14 year", 
               date_labels = "%Y") +
  scale_color_manual(name = "Old Lineage",
                     values = c("#FF0000", "#228B22", "#FFDAB9", "#CD00CD"),
                     breaks = unique(sort(merged_lineages$Old_Lineage))) +
  theme_tree2() +
  theme(text = element_text(family = "Arial Narrow"),
        axis.text.x = element_text(size = 18, color = "#000000"),
        legend.text = element_text(size = 14, color = "#000000"),
        legend.title = element_text(size = 15, color = "#000000"),
        plot.title = element_text(hjust = 0.5),
        legend.background = element_blank()) +
  guides(fill = guide_legend(nrow = 4, 
                             override.aes = list(size = 9), 
                             title = "Country"))
p10

# subset nodes that have > 1 posterior probability
c11_tree <- tree_subset(mcc.tre.lineageC.medium, 179, levels_back = 0)
c11_sdn_tree <- tree_subset(mcc.tre.lineageC.medium, 223, levels_back = 0)
c12_tree <- tree_subset(mcc.tre.lineageC.medium, 188, levels_back = 0)
c21_tree <- tree_subset(mcc.tre.lineageC.medium, 169, levels_back = 0)
c22_tree <- tree_subset(mcc.tre.lineageC.medium, 125, levels_back = 0)

merged_lineages$sublineage[merged_lineages$taxa %in% c11_tree@phylo$tip.label] <- "C.1.1"
merged_lineages$sublineage[merged_lineages$taxa %in% c11_sdn_tree@phylo$tip.label] <- "C.1.1"
merged_lineages$sublineage[merged_lineages$taxa %in% c12_tree@phylo$tip.label] <- "C.1.2"
merged_lineages$sublineage[merged_lineages$taxa %in% c21_tree@phylo$tip.label] <- "C.2.1"
merged_lineages$sublineage[merged_lineages$taxa %in% c22_tree@phylo$tip.label] <- "C.2.2"
merged_lineages <- merged_lineages %>% dplyr::mutate(sublineage = ifelse(is.na(sublineage), Lineage.x, sublineage))

merged_lineages <- merged_lineages %>% 
  dplyr::select(c("taxa", "accession", "Lineage.x", "country", "countryname",
                  "host", "distance", "date2", "year", "uid",
                  "Old_Lineage", "sublineage")) %>% 
  dplyr::rename("Lineage"="Lineage.x")

# read and and plot tree with snp data
snps_fn <- file.path(wd, "M-lineage-C.snps.csv")
m.lineage.c.snps <- snps_data(snps_fn = snps_fn)

p12 <- plot_mcc_lineage_tree(tree = mcc.tre.lineageC.medium,
                             metadata = merged_lineages,
                             cols = lineageC.cols,
                             subtitle = "RVFV-M Lineage C maximum clade credibility tree",
                             column = "sublineage")
p12 

m.lineage.c.snps.plot <- plot_tree_snps(snps_data = m.lineage.c.snps,
                                        tree = p12,
                                        column = "sublineage",
                                        cols = lineageC.cols)
m.lineage.c.snps.plot <- m.lineage.c.snps.plot +
  theme(strip.text.x = element_text(size=25, face="bold"))

mcc_m_loc <- plot_mcc_with_locations(tree = mcc.tre.lineageC.medium,
                        metadata = meta.medium.lineageC,
                        subtitle = "RVFV-M")
ggsave(file.path(figuresDir, "FigureS3.pdf"),
       mcc_m_loc,
       width = 48, height = 50, units = "in", 
       limitsize = FALSE, scale = 0.3,
       dpi = 300, bg="white", device = cairo_pdf)

###################### MCC tree S-NSS ########################
metadata_fn <- file.path(wd, "S-global.geolocations.csv")
lineages_fn <- file.path(wd, "S-global.lineages.csv")
tempest_fn <- file.path(wd, "nss-lineage-C.root-to-tip.tsv")
mcc.lineageC.nss <- file.path(wd, "nss-lineage-C.mcc.tree")

meta.nss.lineageC <- prepare_metadata(metadata_file = metadata_fn,
                                      lineages_file = lineages_fn,
                                      tempest_file = tempest_fn)

# fix wrongly assigned sequences
meta.nss.lineageC$Lineage[meta.nss.lineageC$accession == 'KY196500'] <- 'H' # wrongly assigned as C
meta.nss.lineageC$Lineage[meta.nss.lineageC$accession == 'DQ380167'] <- 'G' # wrongly assigned as E
meta.nss.lineageC$Lineage[meta.nss.lineageC$accession == 'DQ380165'] <- 'G' # wrongly assigned as E
meta.nss.lineageC$Lineage[meta.nss.lineageC$accession == 'DQ380168'] <- 'G' # wrongly assigned as E
meta.nss.lineageC$Lineage[meta.nss.lineageC$accession == 'DQ380161'] <- 'G' # wrongly assigned as E
meta.nss.lineageC$Lineage[meta.nss.lineageC$accession == 'DQ380160'] <- 'G' # wrongly assigned as E

p13 <- plot_root_to_tip(tempest_data = meta.nss.lineageC, 
                        root = 1980,
                        r2 = 0.81,
                        n = 141,
                        cor = 0.90,
                        pos = 1983,
                        column = "countryname",
                        segment = "RVFV-NSS Lineage C temporal signal",
                        cols = cols_countries,
                        breaks = '8 year')

p13 <- p13 + theme(legend.position = c(0.3, 0.78),
                   axis.text.x = element_text(size = 25, angle = 90),
                   axis.text.y = element_text(size = 25)) +
  guides(fill = guide_legend(nrow = 6,
                             override.aes = list(size = 9),
                             title = "Country")) 
p13

# read tree
mcc.tre.lineageC.nss <- read.beast(mcc.lineageC.nss)

# merge with old lineages classification
# medium.metadata <- read.csv("~/projects/RVFV/continuous/segments/M/complete/global/merged-sequences/M-global.csv", header = T)
# merged_lineages_nss <- left_join(meta.nss.lineageC, old_lineages, by=c("accession"="Accession"))
meta.nss.lineageC <- meta.nss.lineageC[match(mcc.tre.lineageC.nss@phylo$tip.label, meta.nss.lineageC$taxa), ]
all(mcc.tre.lineageC.nss@phylo$tip.label == meta.nss.lineageC$taxa)

#cols <- cols_countries[names(cols_countries) %in% sort(unique(meta.nss.lineageC$countryname))]
mrsd <- meta.nss.lineageC[rev(order(as.Date(meta.nss.lineageC$date2, format = "%Y-%m-%d"))),]$date2[1]

p.nss <- ggtree(mcc.tre.lineageC.nss,
              mrsd = mrsd,
              as.Date = TRUE,
              size = 0.3) %<+% meta.nss.lineageC +
  geom_nodelab(aes(label=node), vjust=-.5, size=2.5) +
  geom_nodelab(aes(x=branch, label=round(posterior == 1, 1), subset=(posterior == 1)), vjust=-.5, size=4) +
  geom_cladelabel(node=235, label="C.1.1", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.3, fontsize = 6.0) +
  geom_cladelabel(node=240, label="C.1.2", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.08, fontsize = 6.0) +
  geom_cladelabel(node=155, label="C.2.1", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.3, fontsize = 6.0) +
  geom_cladelabel(node=174, label="C.2.2", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.3, fontsize = 6.0) +
  geom_tiplab(aes(label = paste0("italic('", accession, "|", country, "|", year, "')"), offset=0.06), parse = TRUE, size = 6.0, align=F) +
  geom_tippoint(aes(color = Lineage), alpha = 1.0, size=1) +
  scale_x_date(date_breaks= "14 year", 
               date_labels = "%Y") +
  theme_tree2() +
  theme(text = element_text(family = "Arial Narrow"),
        axis.text.x = element_text(size = 18, color = "#000000"),
        legend.text = element_text(size = 14, color = "#000000"),
        legend.title = element_text(size = 15, color = "#000000"),
        plot.title = element_text(hjust = 0.5),
        legend.background = element_blank()) +
  guides(fill = guide_legend(nrow = 4, 
                             override.aes = list(size = 9), 
                             title = "Country"))
p.nss

# subset nodes that have > 1 posterior probability
c11_tree <- tree_subset(mcc.tre.lineageC.nss, 235, levels_back = 0)
c11_sdn_tree <- tree_subset(mcc.tre.lineageC.nss, 232, levels_back = 0)
c12_tree <- tree_subset(mcc.tre.lineageC.nss, 240, levels_back = 0)
c21_tree <- tree_subset(mcc.tre.lineageC.nss, 155, levels_back = 0)
c22_tree <- tree_subset(mcc.tre.lineageC.nss, 174, levels_back = 0)

meta.nss.lineageC$sublineage[meta.nss.lineageC$taxa %in% c11_tree@phylo$tip.label] <- "C.1.1"
meta.nss.lineageC$sublineage[meta.nss.lineageC$taxa %in% c11_sdn_tree@phylo$tip.label] <- "C.1.1"
meta.nss.lineageC$sublineage[meta.nss.lineageC$taxa %in% c12_tree@phylo$tip.label] <- "C.1.2"
meta.nss.lineageC$sublineage[meta.nss.lineageC$taxa %in% c21_tree@phylo$tip.label] <- "C.2.1"
meta.nss.lineageC$sublineage[meta.nss.lineageC$taxa %in% c22_tree@phylo$tip.label] <- "C.2.2"
meta.nss.lineageC <- meta.nss.lineageC %>% dplyr::mutate(sublineage = ifelse(is.na(sublineage), Lineage, sublineage))

# read and process snps data
snps_fn <- file.path(wd, "nss-lineage-C.snps.csv")
nss.lineage.c.snps <- snps_data(snps_fn = snps_fn)

p14 <- plot_mcc_lineage_tree(tree = mcc.tre.lineageC.nss,
                             metadata = meta.nss.lineageC,
                             cols = lineageC.cols,
                             subtitle = "RVFV-NSS Lineage C maximum clade credibility",
                             column = "sublineage")
p14

nss.lineage.c.snps.plot <- plot_tree_snps(snps_data = nss.lineage.c.snps,
                                        tree = p14,
                                        column = "sublineage",
                                        cols = lineageC.cols)
nss.lineage.c.snps.plot <- nss.lineage.c.snps.plot + 
  theme(strip.text.x = element_text(size=25, face="bold"))
nss.lineage.c.snps.plot

mcc_nss_loc <- plot_mcc_with_locations(tree = mcc.tre.lineageC.nss,
                        metadata = meta.nss.lineageC,
                        subtitle = "RVFV-NSS")

###################### MCC tree L ########################
metadata_fn <- file.path(wd, "L-global.geolocations.csv")
lineages_fn <- file.path(wd, "L-global.lineages.csv")
tempest_fn <- file.path(wd, "L-lineage-C.root-to-tip.tsv")
mcc.lineageC.large <- file.path(wd, "L-lineage-C.mcc.tree")

meta.large.lineageC <- prepare_metadata(metadata_file = metadata_fn,
                                        lineages_file = lineages_fn,
                                        tempest_file = tempest_fn)
# fix wrongly assigned sequences
meta.large.lineageC$Lineage[meta.large.lineageC$accession == 'DQ375434'] <- 'J' # wrongly assigned as I
meta.large.lineageC$Lineage[meta.large.lineageC$accession == 'MG659874'] <- 'H' # wrongly assigned as C
meta.large.lineageC$Lineage[meta.large.lineageC$accession == 'MG659825'] <- 'H' # wrongly assigned as C
meta.large.lineageC$Lineage[meta.large.lineageC$accession == 'MG659886'] <- 'H' # wrongly assigned as C
meta.large.lineageC$Lineage[meta.large.lineageC$accession == 'KY126678'] <- 'H' # wrongly assigned as C
meta.large.lineageC$Lineage[meta.large.lineageC$accession == 'DQ375425'] <- 'D' # wrongly assigned as C
meta.large.lineageC$Lineage[meta.large.lineageC$accession == 'DQ375433'] <- 'O' # wrongly assigned as L

p15 <- plot_root_to_tip(tempest_data = meta.large.lineageC,
                        root = 1964,
                        r2 = 0.91,
                        n = 114,
                        cor = 0.95,
                        pos = 1968,
                        column = "countryname",
                        segment = "RVFV-L Lineage C temporal signal",
                        cols = cols_countries,
                        breaks = '8 year')
p15 <- p15 + theme(legend.position = c(0.3, 0.78),
                   axis.text.x = element_text(size = 25, angle = 90),
                   axis.text.y = element_text(size = 25)) +
  guides(fill = guide_legend(nrow = 5,
                             override.aes = list(size = 9),
                             title = "Country"))
p15

# read tree
mcc.tre.lineageC.large <- read.beast(mcc.lineageC.large)

# merge with old lineages classification
# medium.metadata <- read.csv("~/projects/RVFV/continuous/segments/M/complete/global/merged-sequences/M-global.csv", header = T)
# merged_lineages_large <- left_join(meta.large.lineageC, old_lineages, by=c("accession"="Accession"))
meta.large.lineageC <- meta.large.lineageC[match(mcc.tre.lineageC.large@phylo$tip.label, meta.large.lineageC$taxa), ]
all(mcc.tre.lineageC.large@phylo$tip.label == meta.large.lineageC$taxa)

# cols <- cols_countries[names(cols_countries) %in% sort(unique(meta.large.lineageC$countryname))]
mrsd <- meta.large.lineageC[rev(order(as.Date(meta.large.lineageC$date2, format = "%Y-%m-%d"))),]$date2[1]


p.large <- ggtree(mcc.tre.lineageC.large,
                  mrsd = mrsd,
                  as.Date = TRUE,
                  size = 0.3) %<+% meta.large.lineageC +
  geom_nodelab(aes(label=node), vjust=-.5, size=2.5) +
  geom_nodelab(aes(x=branch, label=round(posterior == 1, 1), subset=(posterior == 1)), vjust=-.5, size=4) +
  geom_cladelabel(node=177, label="C.1.1", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.3, fontsize = 6.0) +
  geom_cladelabel(node=183, label="C.1.2", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.08, fontsize = 6.0) +
  geom_cladelabel(node=171, label="C.2.1", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.3, fontsize = 6.0) +
  geom_cladelabel(node=124, label="C.2.2", color="#FF0000", barsize = 3.5, align=TRUE, offset=0.3, fontsize = 6.0) +
  # geom_tiplab(aes(label = paste0("italic('", accession, "|", country, "|", year, "')"), offset=0.06), parse = TRUE, size = 6.0, align=F) +
  geom_tippoint(aes(color = Lineage), alpha = 1.0, size=1) +
  scale_x_date(date_breaks= "14 year", 
               date_labels = "%Y") +
  # scale_color_manual(name = "Old Lineage",
  #                    values = c("#FF0000", "#228B22", "#FFDAB9", "#CD00CD"),
  #                    breaks = unique(sort(merged_lineages$Old_Lineage))) +
  theme_tree2() +
  theme(text = element_text(family = "Arial Narrow"),
        axis.text.x = element_text(size = 18, color = "#000000"),
        legend.text = element_text(size = 14, color = "#000000"),
        legend.title = element_text(size = 15, color = "#000000"),
        plot.title = element_text(hjust = 0.5),
        legend.background = element_blank()) +
  guides(fill = guide_legend(nrow = 4, 
                             override.aes = list(size = 9), 
                             title = "Country"))
p.large

# subset nodes that have > 1 posterior probability
c11_tree <- tree_subset(mcc.tre.lineageC.large, 177, levels_back = 0)
c11_sdn_tree <- tree_subset(mcc.tre.lineageC.large, 175, levels_back = 0)
c12_tree <- tree_subset(mcc.tre.lineageC.large, 183, levels_back = 0)
c21_tree <- tree_subset(mcc.tre.lineageC.large, 171, levels_back = 0)
c22_tree <- tree_subset(mcc.tre.lineageC.large, 124, levels_back = 0)

meta.large.lineageC$sublineage[meta.large.lineageC$taxa %in% c11_tree@phylo$tip.label] <- "C.1.1"
meta.large.lineageC$sublineage[meta.large.lineageC$taxa %in% c11_sdn_tree@phylo$tip.label] <- "C.1.1"
meta.large.lineageC$sublineage[meta.large.lineageC$taxa %in% c12_tree@phylo$tip.label] <- "C.1.2"
meta.large.lineageC$sublineage[meta.large.lineageC$taxa %in% c21_tree@phylo$tip.label] <- "C.2.1"
meta.large.lineageC$sublineage[meta.large.lineageC$taxa %in% c22_tree@phylo$tip.label] <- "C.2.2"
meta.large.lineageC <- meta.large.lineageC %>% dplyr::mutate(sublineage = ifelse(is.na(sublineage), Lineage, sublineage))


# read and process snps data
snps_fn <- file.path(wd, "L-lineage-C.snps.csv")
large.lineage.c.snps <- snps_data(snps_fn = snps_fn)

p16 <- plot_mcc_lineage_tree(tree = mcc.tre.lineageC.large,
                             metadata = meta.large.lineageC,
                             cols = lineageC.cols,
                             subtitle = "RVFV-L Lineage C maximum clade credibility tree",
                             column = "sublineage")
p16
l.lineage.c.snps.plot <- plot_tree_snps(snps_data = large.lineage.c.snps,
                                            tree = p16,
                                            column = "sublineage",
                                            cols = lineageC.cols)
l.lineage.c.snps.plot <- l.lineage.c.snps.plot +
  theme(strip.text.x = element_text(size=25, face="bold"))

mcc_l_loc <- plot_mcc_with_locations(tree = mcc.tre.lineageC.large,
                        metadata = meta.large.lineageC,
                        subtitle = "RVFV-L")

tres <- ((l.lineage.c.snps.plot + theme(legend.position = "none",
                                        axis.text.x = element_text(angle = 90, 
                                                                   size = 25,
                                                                   color = "#000000",
                                                                   face = "bold"),
                                        plot.subtitle = element_blank())) |
           (m.lineage.c.snps.plot + 
              theme(axis.text.x = element_text(angle = 90, 
                                               size = 25,
                                               color = "#000000",
                                               face = "bold"))) |
           (nss.lineage.c.snps.plot + 
              theme(legend.position = "none",
                    axis.text.x = element_text(angle = 90, 
                                               size = 25,
                                               color = "#000000",
                                               face = "bold"),
                    plot.subtitle = element_blank()))) +
  plot_layout(guides = "collect",
              widths = c(90, 90, 90))

# plot occurence of sublineage
plot_occurence_sublineage <- function(data, segment){
  data <- data %>% dplyr::group_by(sublineage, countryname, year) %>% 
    dplyr::summarise(count=dplyr::n())
  data$year <- as.Date(as.character(data$year), format = "%Y")
  
  p1 <- ggplot(data=data) +
    geom_point(aes(x=year, fill=sublineage, 
                   y=reorder(countryname, count)), 
               position = position_jitter(width=0.2, height=0.2), 
               shape=21, stroke=0.2, col='grey50', size=12, alpha=0.8)+
    scale_x_date(date_labels = "%Y",date_breaks = "8 year") +
    scale_fill_manual(values=lineageC.cols, name='Sublineage') +
    ylab('country') + 
    xlab('year') +
    ggtitle("", subtitle = segment) +
    theme_bw() +
    theme(
      text = element_text(family = "Arial Narrow"),
      plot.tag = element_text(size = 30, color = "#000000", face = "bold"),
      plot.subtitle=element_text(size=25, color = "#000000", face = "bold"),
      axis.title = element_blank(),
      axis.text = element_text(size=18, color = "#000000", face = "bold"),
      legend.position="right",
      legend.text = element_text(size=16, color = "#000000")) +
    guides(fill = guide_legend(override.aes = list(size=10)))
  print(p1)
  return(p1)
}

occurence.L <- plot_occurence_sublineage(data = meta.large.lineageC,
                                         segment = "RVFV-L lineage C sublineages occurence")
occurence.M <- plot_occurence_sublineage(data = merged_lineages,
                                         segment = "RVFV-M lineage C sublineages occurence") 
occurence.M <- occurence.M +
  theme(legend.position = "none")
occurence.M
occurence.nss <- plot_occurence_sublineage(data = meta.nss.lineageC,
                                           segment = "RVFV-NSS lineage C sublineages occurence")
occurence.nss <- occurence.nss +
  theme(legend.position = "none")
occurence.nss


# plot subs rates
# read log files
mean_rates <- function(log, segment){
  df <- read.csv(log, skip = 4, header = T, sep = "\t") %>% 
    dplyr::select(c("meanRate"))
  df$segment <- segment
  return(df)
}

log_l <- file.path(wd, "L-lineage-C-phylo.log")
log_m <- file.path(wd, "M-lineage-C-phylo.log")
log_nss <- file.path(wd, "nss-lineage-C-phylo.log")

log.l.df <- mean_rates(log = log_l, segment = "L")
log.m.df <- mean_rates(log = log_m, segment = "M")
log.nss.df <- mean_rates(log = log_nss, segment = "NSS")

d <- bind_rows(log.l.df, log.m.df, log.nss.df)
d$segment <- as.factor(d$segment)

# compute lower and upper whiskers
segs.cols <- c("#8B008B", "#008B00", "#7ad2f6")
subs_plot <- ggplot(d, aes(x = meanRate, y = segment, fill=segment)) + 
  geom_density_ridges(rel_min_height=.01, scale=1.0, alpha=0.8) +
  xlim(0,0.002) +
  scale_fill_manual(name='Segment',
                      values=segs.cols,
                      breaks = unique(d$segment),
                      guide = guide_legend(override.aes = list(linetype = 0,
                                                               shape = 21,
                                                               color = segs.cols))) +
  labs(title = "RVFV Lineage C evolution rate",
       y = "segment",
       x = "subs/site/year") +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(text = element_text(family = "Arial Narrow"),
        panel.border  = element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.title=element_text(size=25, color = "#000000", face = "bold"),
        plot.tag = element_text(size = 30, color = "#000000", face = "bold"),
        axis.text = element_text(size = 20, color = "#000000", face = "bold"),
        axis.title = element_text(size = 25, color = "#000000"),
        legend.text = element_text(size = 20, color = "#000000"),
        legend.title = element_text(size = 20, color = "#000000"),
        legend.position = "none") +
  guides(color = guide_legend(nrow = 3, override.aes = list(size = 14)))
subs_plot

figure.2.plots <- ((tres[[2]] + theme(axis.text.x = element_text(angle = 0))) | (occurence.M / subs_plot)) +
  plot_layout(widths = c(20, 22),
              heights = c(15, 20),
              byrow = T) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 30, colour = "#000000", face = "bold"),
        plot.subtitle = element_text(size = 26, colour = "#000000", face = "bold"),
        axis.text = element_text(size = 24, colour = "#000000", face = "bold"),
        strip.text.x = element_text(size=24, colour = "#000000", face="bold"),
        legend.title = element_text(size = 22, colour = "#000000", face = "bold"))

ggsave(file.path(figuresDir, "Figure2.pdf"), 
       figure.2.plots,
       width = 24, height = 16, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDir, "Figure2.png"), 
       figure.2.plots,
       width = 24, height = 16, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = "png")

ggsave(file.path(figuresDir, "Figure2.tiff"), 
       figure.2.plots,
       width = 24, height = 16, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = "tiff")



sf2 <- 
  ((((p9 + theme(legend.position = "none")) + 
    (p13 + theme(axis.title.y = element_blank(),
                 legend.position = "none")) + 
    (p15 + theme(axis.title.y = element_blank(),
                 legend.position = "right") +
       guides(fill = guide_legend(title = "country",
                                    nrow = 12, 
                                    override.aes = list(size = 15))))) +
     plot_layout(guides = "collect")) / 
  ((p14 + theme(axis.text.y = element_blank(),
                legend.position = "none")) | 
     (p16 + theme(axis.text.y = element_blank(),
                  legend.position = "right"))) / 
  ((occurence.nss + theme(legend.position = "none")) | 
     (occurence.L + theme(legend.position = "none")))) +
  plot_layout(widths = c(8, 7.5, 7.5),
              heights = c(5, 8, 5),
              byrow = T) +
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 30, colour = "#000000", face = "bold"),
        plot.subtitle = element_text(size = 24, colour = "#000000", face = "bold"),
        axis.text = element_text(size = 23, colour = "#000000", face = "bold"),
        axis.text.x = element_text(size = 23, colour = "#000000", angle = 0),
        legend.title = element_text(size = 22, colour = "#000000"))

ggsave(file.path(figuresDir, "FigureS2.pdf"), 
       sf2,
       width = 30, height = 24, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDir, "FigureS2.png"), 
       sf2,
       width = 30, height = 24, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = "png")

ggsave(file.path(figuresDir, "FigureS2.tiff"), 
       sf2,
       width = 30, height = 24, units = "in", 
       limitsize = FALSE,
       dpi = 300, bg="white", device = "tiff")

ggsave(file.path(figuresDir, "FigureS3.pdf"),
       mcc_m_loc,
       width = 48, height = 50, units = "in", 
       limitsize = FALSE, scale = 0.3,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDir, "FigureS3.png"),
       mcc_m_loc,
       width = 48, height = 50, units = "in", 
       limitsize = FALSE, scale = 0.3,
       dpi = 300, bg="white", device = "png")

ggsave(file.path(figuresDir, "FigureS4.pdf"),
       mcc_nss_loc,
       width = 48, height = 50, units = "in", 
       limitsize = FALSE, scale = 0.3,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDir, "FigureS4.png"),
       mcc_nss_loc,
       width = 48, height = 50, units = "in", 
       limitsize = FALSE, scale = 0.3,
       dpi = 300, bg="white", device = "png")

ggsave(file.path(figuresDir, "FigureS5.pdf"),
       mcc_l_loc,
       width = 48, height = 50, units = "in", 
       limitsize = FALSE, scale = 0.3,
       dpi = 300, bg="white", device = cairo_pdf)

ggsave(file.path(figuresDir, "FigureS5.png"),
       mcc_l_loc,
       width = 48, height = 50, units = "in", 
       limitsize = FALSE, scale = 0.3,
       dpi = 300, bg="white", device = "png")
