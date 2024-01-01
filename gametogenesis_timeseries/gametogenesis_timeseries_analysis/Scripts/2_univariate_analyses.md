Univariate analysis of heatwave time series biological data
================
Ariana S Huffmyer, E5 RoL Team
20220830

# Set Up

``` r
## install packages if you dont already have them
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
if (!require("lme4")) install.packages("lme4")
if (!require("lmerTest")) install.packages("lmerTest")
if (!require("car")) install.packages("car")
if (!require("effects")) install.packages("effects")
if (!require("ggfortify")) install.packages("ggfortify")
if (!require("cowplot")) install.packages("cowplot")
if (!require("vegan")) install.packages("vegan")
if (!require("corrr")) install.packages("corrr")
if (!require("ggcorrplot")) install.packages("ggcorrplot")
if (!require("GGally")) install.packages("GGally")
if (!require("broom")) install.packages("broom")
if (!require("cowplot")) install.packages("cowplot")

# load packages
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(car)
library(effects)
library(ggfortify)
library(cowplot)
library(vegan)
library(corrr)
library(ggcorrplot)
library(GGally)
library(broom)
library(cowplot)
```

# Load dataframe

Load in master dataframe generated from 1\_assemble\_data.Rmd.

``` r
master<-read.csv("gametogenesis_timeseries_analysis/Output/master_timeseries.csv")
```

# Univariate Responses

Plot univariate response plots and then generate panel of all plots at
the end of this section. Individual plots will not be displayed.

### Plot: Symbiont Densities

View by treatment and timepoint.

``` r
symbplot_cm2<-master %>%
  filter(!is.na(cells.cm2)) %>%
  filter(!is.na(timepoint))%>%
  select(colony_id, timepoint, cells.cm2)%>%
  
  ggplot(., aes(x = timepoint, y = cells.cm2, fill=timepoint)) +
    geom_boxplot(color="black", fill="white", alpha=0.2) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
  scale_x_discrete(limits = c("JAN","FEB","MARCH", "APRIL", "MAY", "JUNE", "JULY")) +
    xlab("") + 
    ylab(expression(bold(paste("Symbiont Cells cm"^-2))))+
    theme_classic() + 
    theme(legend.position = "none",
      axis.title=element_text(face="bold", size=30),
      axis.text=element_text(size=18, color="black")
      )
```

### Plot: Host Protein

View by treatment and timepoint.

``` r
host_prot_cm2<-master %>%
  filter(!is.na(host_prot_mg.cm2)) %>%
  filter(!is.na(timepoint))%>%
  select(colony_id, timepoint, host_prot_mg.cm2)%>%
  
  ggplot(., aes(x = timepoint, y = host_prot_mg.cm2, fill = timepoint)) +
    geom_boxplot(color="black", fill="white", alpha=0.2) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) + 
    ylim(0, 3.5)+
    scale_x_discrete(limits = c("JAN","FEB","MARCH", "APRIL", "MAY", "JUNE", "JULY")) +
    xlab("") + 
    ylab(expression(bold(paste("Host Protein ( mg cm"^-2, ")"))))+
    theme_classic() + 
    theme(legend.position = "none",
      axis.title=element_text(face="bold", size=30),
      axis.text=element_text(size=18, color="black")
      )
```

### Plot: Holobiont protein

View by treatment and timepoint.

``` r
holobiont_prot_cm2<-master %>%
  filter(!is.na(holobiont_prot_mg.cm2)) %>%
  filter(!is.na(timepoint))%>%
  select(colony_id, timepoint, holobiont_prot_mg.cm2)%>%
  
  ggplot(., aes(x = timepoint, y = holobiont_prot_mg.cm2, fill = timepoint)) +
    geom_boxplot(color="black", fill="white", alpha=0.2) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) +
    ylim(0, 3.5)+
    scale_x_discrete(limits = c("JAN","FEB","MARCH", "APRIL", "MAY", "JUNE", "JULY")) +
    xlab("") +
    ylab(expression(bold(paste("Holobiont Protein (mg cm"^-2, ")"))))+
    theme_classic() + 
    theme(legend.position = "none",
      axis.title=element_text(face="bold", size=30),
      axis.text=element_text(size=18, color="black")
      )
```

### Plot: CHL-a

View by treatment and timepoint.

``` r
chla_cm2<-master %>%
  filter(!is.na(chla.ug.cm2)) %>%
  filter(!is.na(timepoint))%>%
  select(colony_id, timepoint, chla.ug.cm2)%>%
  
  ggplot(., aes(x = timepoint, y = chla.ug.cm2, fill = timepoint)) +
    geom_boxplot(color="black", fill="white", alpha=0.2) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) +
    ylim(0, 12)+
    scale_x_discrete(limits = c("JAN","FEB","MARCH", "APRIL", "MAY", "JUNE", "JULY")) +
    xlab("") + 
    ylab(expression(bold(paste("Chlorophyll a (", mu, "g cm"^-2, ")"))))+
    theme_classic() + 
    theme(legend.position = "none",
      axis.title=element_text(face="bold", size=30),
      axis.text=element_text(size=18, color="black")
      )
```

### Plot: CHL-c

View by treatment and timepoint.

``` r
chlc_prot_cm2<-master %>%
  filter(!is.na(chlc2.ug.cm2)) %>%
  filter(!is.na(timepoint))%>%
  select(colony_id, timepoint, chlc2.ug.cm2)%>%
  
  ggplot(., aes(x = timepoint, y = chlc2.ug.cm2, fill = timepoint)) +
    geom_boxplot(color="black", fill="white", alpha=0.2) +
    geom_point(pch = 21, size=2, position = position_jitterdodge(0.2)) +
    ylim(0, 12)+
    scale_x_discrete(limits = c("JAN","FEB","MARCH", "APRIL", "MAY", "JUNE", "JULY")) +
    xlab("") + 
    ylab(expression(bold(paste("Chlorophyll c2 (", mu, "g cm"^-2, ")"))))+
    theme_classic() + 
    theme(legend.position = "none",
      axis.title=element_text(face="bold", size=30),
      axis.text=element_text(size=18, color="black")
      )
```

Join all plots for each normalization together.

``` r
phys_sym_figure<-plot_grid(symbplot_cm2, chla_cm2, chlc_prot_cm2, ncol=3, nrow=1, labels = c('A', 'B', 'C'),label_size = 36,rel_heights= 1, rel_widths = 1, label_y=1, align="h")

ggsave(filename="gametogenesis_timeseries_analysis/Figures/Phys_Sym_Figure.pdf", plot=phys_sym_figure, dpi=300, width=24, height=10, units="in")

phys_host_figure<-plot_grid( host_prot_cm2, holobiont_prot_cm2, ncol=2, nrow=1, labels = c('A', 'B'),label_size = 36,rel_heights= 1, rel_widths = 1, label_y=1, align="h")

ggsave(filename="gametogenesis_timeseries_analysis/Figures/Phys_Host_Figure.pdf", plot=phys_host_figure, dpi=300, width=24, height=10, units="in")
```
