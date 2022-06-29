# ------------------------------------------------------------------------------
# 
# Simulation Multigroup Analysis using PLS-PM
#
# (c) Michael Klesel, Florian Schuberth, Bjoern Niehaves, Joerg Henseler
#
# Klesel, M., Schuberth, F., Henseler, J., and Niehaves, B. (forthcoming) 
# “Multigroup Analysis in Information Systems Research using PLS-PM: 
#  A Systematic Investigation of Approaches,“ The DATA BASE for 
#  Advances in Information Systems
# 
# ------------------------------------------------------------------------------

# Prepare workspace ------------------------------------------------------------
rm(list = ls())
options(scipen=999)
graphics.off()


# Libraries --------------------------------------------------------------------
library("tidyverse")     
library("scales")
library("svglite")


library("renv")
# renv::init()

# serialize current state
# renv::snapshot()

# Restore from renv
# renv::restore()


# User functions ---------------------------------------------------------------

# Add Confidence Intervals
addCI <- function(data, runs, colum, all){
  # Add Confidence Intervals
  data <- data %>% mutate_at(colum, 
                             .f = list(lowerBound = ~ .-sqrt(.*(1-.)/runs*qnorm(0.975)),
                                       upperBound = ~ .+sqrt(.*(1-.)/runs*qnorm(0.975))))

  return(data)
}

# http://jfly.iam.u-tokyo.ac.jp/color
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
myBlack   <- "#000000"
myOrange  <- "#E69F00"
myBlue    <- "#0072B2"
myGreen   <- "#009E73"


data.plot <- readRDS("data/Results.rds")


#  -----------------------------------------------------------------------------
# SINGLE PATH
#  -----------------------------------------------------------------------------

# Type I Error -----------------------------------------------------------------

data.plot1 <- data.plot %>% 
  filter(differences == "none",                       # homogenous groups only
         normalData == F,                             # non-normal distributed data only
         correction == "holm" | is.na(correction) |   # use holm for multiple groups
           (test == "NBT" & correction == "none"),      # Include Henseler
         type_ci == "CI_bca" | is.na(type_ci),        # bca                              
         what == "path",                              # path only
         (is.na(distance) | distance == "dG"),        # distance
         ssizediff == "moderately unequal")           # moderately unequal sample sizes

# Long Format
data.plot1 <- data.plot1 %>% 
  pivot_longer(cols = c(`1%_perCentage`,`5%_perCentage`,`10%_perCentage`), 
               values_to = "RR", names_to = "alpha") %>%
  mutate(alpha = factor(alpha, 
                        levels = c("10%_perCentage", 
                                   "5%_perCentage", 
                                   "1%_perCentage"))) %>%
  mutate(alpha = forcats::fct_recode(alpha, 
                                     `10%` = "10%_perCentage",
                                     `5%`  = "5%_perCentage", 
                                     `1%`  = "1%_perCentage")) %>%
  addCI(data = ., runs = data.plot1$runs[1], colum = "RR") %>% 
  filter(alpha == "5%")

p1 <- ggplot(data.plot1 , aes(y = RR, x = ssize, 
                              group = test)) +
  facet_grid(cols = vars(test), rows = vars(groups))+
  geom_point()+
  geom_line() +
  geom_hline(yintercept=c(0.05), linetype="dashed", color = "black") +
  geom_line(aes(y=upperBound), linetype='dotted')+
  geom_line(aes(y=lowerBound), linetype='dotted')+
  geom_ribbon(aes(ymin = lowerBound, ymax = upperBound), alpha = 0.15)+
  scale_y_continuous(name="Rejection rates", 
                     limits=c(0, 1), 
                     breaks=seq(0,1,.2),
                     labels = scales::percent) + 
  scale_x_discrete(name = "Total sample size") +
  theme_bw(base_size = 14) +
  theme(legend.position="bottom")
p1

# Figure 2
ggsave(p1, device = "svg", filename = "images/Figure02.svg", width = 30, height = 16)


# Power 
data.plot2 <- data.plot %>% 
  filter(differences != "none",                            # homogenous groups only
         normalData == F,                             # non-normal distributed data only
         correction == "holm" | is.na(correction) |   # use holm for multiple groups
           (test == "NBT" & correction == "none"),      # Include Henseler
         type_ci == "CI_bca" | is.na(type_ci),        # bca                              # path only
         what == "path",
         (is.na(distance) | distance == "dG"),
         ssizediff == "moderately unequal")   

data.plot2 <- data.plot2 %>% 
  pivot_longer(cols = c(`1%_perCentage`,`5%_perCentage`,`10%_perCentage`), 
               values_to = "RR", names_to = "alpha") %>%
  rename("Parameter difference" = differences) %>%
  filter(alpha == "5%_perCentage")

p2 <- ggplot(data.plot2 , aes(y = RR, x = ssize, 
                              group = interaction(groups, `Parameter difference`),
                              color = `Parameter difference`)) +
  facet_grid(cols = vars(test), rows = vars(groups))+
  geom_point(aes(shape = `Parameter difference`, color = `Parameter difference`))+
  # Define shapes
  scale_shape_manual(values=c(0, 8, 10, 17))+
  # Define colors
  scale_color_manual(values=c(myOrange,myGreen,myBlue, myBlack))+
  geom_line()+
  # geom_line(aes(linetype = differences)) +
  # Define line-types
  # https://stackoverflow.com/questions/52885265/change-line-width-in-ggplot-not-size
  # Define different linetypes
  # The first numeral is units of dash length, 
  # the second units in the gap in hexadecimal. 
  # scale_linetype_manual(values=c("81", "82","83","84"))+
  geom_hline(yintercept=c(0.8), linetype="dotted", color = "black") +
  scale_y_continuous(name="Rejection rates", 
                     limits=c(0, 1), 
                     breaks=seq(0,1,.2),
                     labels = scales::percent) + 
  scale_x_discrete(name = "Total sample size") + 
  theme_bw(base_size = 14) +
  theme(legend.position="bottom") 
p2

# Figure 3
ggsave(p2, device = "svg", filename = "images/Figure03.svg", width = 18, height = 10)


#  -----------------------------------------------------------------------------
#  COMPLETE MODEL
#  -----------------------------------------------------------------------------

# Type I Error -----------------------------------------------------------------

data.plot3 <- data.plot %>% 
  filter(differences == "none",                       # homogenous groups only
         # test != "OTG",                             # exclude OTG
         test != "NBT",                               # exclude NBT
         normalData == F,                             # non-normal distributed data only
         correction == "holm" | is.na(correction),    # use holm for multiple groups
         type_ci == "CI_bca" | is.na(type_ci),        # bca                              
         what == "complete",                          # path only
         (is.na(distance) | distance == "dG"),        # distance
         ssizediff == "moderately unequal")           # moderately unequal sample sizes

# Long Format
data.plot3 <- data.plot3 %>% 
  pivot_longer(cols = c(`1%_perCentage`,`5%_perCentage`,`10%_perCentage`), 
               values_to = "RR", names_to = "alpha") %>%
  mutate(alpha = factor(alpha, 
                        levels = c("10%_perCentage", 
                                   "5%_perCentage", 
                                   "1%_perCentage"))) %>%
  mutate(alpha = forcats::fct_recode(alpha, 
                                     `10%` = "10%_perCentage",
                                     `5%`  = "5%_perCentage", 
                                     `1%`  = "1%_perCentage")) %>%
  addCI(data = ., runs = data.plot3$runs[1], colum = "RR") %>%
  filter(alpha == "5%")



p3 <- ggplot(data.plot3 , aes(y = RR, x = ssize, group = test)) +
  facet_grid(cols = vars(test), rows = vars(groups))+
  geom_point()+
  geom_line() +
  geom_hline(yintercept=c(0.05), linetype="dashed", color = "black") +
  geom_line(aes(y=upperBound), linetype='dotted')+
  geom_line(aes(y=lowerBound), linetype='dotted')+
  geom_ribbon(aes(ymin = lowerBound, ymax = upperBound), alpha = 0.15)+
  scale_y_continuous(name="Rejection rates", 
                     limits=c(0, 1), 
                     breaks=seq(0,1,.2),
                     labels = scales::percent) + 
  scale_x_discrete(name = "Total sample size") +
  theme_bw(base_size = 18) +
  theme(legend.position="bottom")
p3

# Figure 4
ggsave(p3, device = "svg", filename = "images/Figure04.svg", width = 14, height = 10)


# Power 
data.plot4 <- data.plot %>% 
  filter(differences != "none",                       # homogenous groups only
         # test != "OTG",                             # exclude OTG
         test != "NBT",                               # exclude NBT
         normalData == F,                             # non-normal distributed data only
         correction == "holm" | is.na(correction),    # use holm for multiple groups
         type_ci == "CI_bca" | is.na(type_ci),        # bca
         distance == "dG" | is.na(distance),          # use geodesic distance only
         what == "complete",
         ssizediff == "moderately unequal")          # severe path differenc

data.plot4 <- data.plot4 %>% 
  pivot_longer(cols = c(`1%_perCentage`,`5%_perCentage`,`10%_perCentage`), 
               values_to = "RR", names_to = "alpha") %>%
  filter(alpha == "5%_perCentage")

data.plot4 <- data.plot4 %>% rename("Structural model difference" = differences)


p4 <- ggplot(data.plot4 , aes(y = RR, x = ssize, 
                              group = interaction(groups, `Structural model difference`),
                              color = `Structural model difference`)) +
  facet_grid(cols = vars(test), rows = vars(groups))+
  geom_point(aes(shape = `Structural model difference`, color = `Structural model difference`))+
  # Define shapes
  scale_shape_manual(values=c(0, 8, 10, 17))+
  # Define colors
  scale_color_manual(values=c(myOrange,myGreen,myBlue, myBlack))+
  geom_line()+
  # geom_line(aes(linetype = `Structural model difference`)) +
  # Define line-types
  # https://stackoverflow.com/questions/52885265/change-line-width-in-ggplot-not-size
  # Define different linetypes
  # The first numeral is units of dash length, 
  # the second units in the gap in hexadecimal. 
  # scale_linetype_manual(values=c("81", "82","83","84"))+
  geom_hline(yintercept=c(0.8), linetype="dotted", color = "black") +
  scale_y_continuous(name="Rejection rates", 
                     limits=c(0, 1), 
                     breaks=seq(0,1,.2),
                     labels = scales::percent) + 
  scale_x_discrete(name = "Total sample size") + 
  theme_bw(base_size = 18) +
  theme(legend.position="bottom")
p4

# Figure 6
ggsave(p4, device = "svg", filename = "images/Figure06.svg", width = 18, height = 10)



#  -----------------------------------------------------------------------------
#  MULTIPLE COMPARISON ISSUE
#  -----------------------------------------------------------------------------

data.plot5 <- data.plot %>% 
  filter(differences == "none" | differences == "medium-large",                            
         test == "PTE" | test == "PTU" | test == "NPT",
         normalData == F,   # non-normal distributed data only
         correction != "fdr" & correction != "BY",
         #  correction == "none" | correction == "bonferroni" | 
         # correction == "hochberg" | correction == "holm" correction == "BH", 
         comparison == "overall",                     # overall comparision only
         ssizediff == "moderately unequal", 
         groups == "3 groups") 

data.plot5 <- data.plot5 %>%
  mutate(correction = case_when(
    correction == "none"       ~ "none",
    correction == "bonferroni" ~ "Bonferroni",
    correction == "holm"       ~ "Holm",
    correction == "hommel"     ~ "Hommel",
    correction == "hochberg"   ~ "Hochberg",
    correction == "BH"         ~ "Benjamini & Hochberg",
    TRUE ~ "ERROR"
  )) %>% dplyr::rename(., c("Adjustment" = "correction"))

data.plot5$differences <- plyr::revalue(data.plot5$differences, c(none = "No group difference"))
data.plot5$differences <- plyr::revalue(data.plot5$differences, c("medium-large" = "Medium-large difference"))

data.plot5.test <- data.plot5 %>% 
  pivot_wider(names_from = Adjustment, values_from = `5%_perCentage`)

# fit <- aov(as.numeric(data.plot5$`5%_perCentage`) ~ as.factor(data.plot5$Adjustment))
# TukeyHSD(fit)

p5 <- ggplot(data.plot5 , aes(y = `5%_perCentage`, x = ssize, 
                              group = interaction(differences, Adjustment),
                              shape = Adjustment)) +
  facet_grid(cols = vars(test), rows = vars(differences))+
  geom_point(aes(shape = Adjustment, color = Adjustment))+
  geom_line()+
  geom_hline(yintercept=c(0.05), linetype="dotted", color = "black") +
  scale_y_continuous(name="Rejection rates", 
                     limits=c(0, 1), 
                     breaks=seq(0,1,.2),
                     labels = scales::percent) + 
  scale_x_discrete(name = "Total sample size") + 
  theme_bw(base_size = 18) +
  theme(legend.position="bottom")

p5

# Image 7
ggsave(p5, device = "svg", filename = "images/Figure07.svg", width = 10, height = 7)


#  -----------------------------------------------------------------------------
#  NON-NORMAL DATA, DATA DISTRIBUTION
#  -----------------------------------------------------------------------------

data.plot6 <- data.plot %>% 
  filter(differences == "small-medium",                            
         test == "PTE" | test == "PTU" | test == "NPT" | test == "NDT",
         # normalData == F,                             # non-normal distributed data only
         correction == "holm" | is.na(correction),
         distance == "dG" | is.na(distance),
         comparison == "overall",                     # overall comparision only
         what == "complete",
         ssizediff != "moderately unequal", 
         groups == "3 groups") 

# rename cols
data.plot6 <- data.plot6 %>% dplyr::rename(., c("Sample size distribution" = "ssizediff", 
                                                "Data distribution" = "normalData" ))

# change values
data.plot6$`Data distribution` <- recode_factor(data.plot6$`Data distribution`,
                                                "FALSE" = "non-normal", 
                                                "TRUE" = "normal")

p6 <- ggplot(data.plot6 , aes(y = `5%_perCentage`, x = ssize, 
                              group = interaction(`Data distribution`, `Sample size distribution`),
                              linetype = `Data distribution`,
                              color = `Sample size distribution`)) +
  facet_grid(cols = vars(test), rows = vars(differences))+
  geom_point(aes(shape = `Sample size distribution`))+
  geom_line()+
  scale_color_manual(values=c(myOrange,myBlue))+
  geom_hline(yintercept=c(0.8), linetype="dotted", color = "black") +
  scale_y_continuous(name="Rejection rates", 
                     limits=c(0, 1), 
                     breaks=seq(0,1,.1),
                     labels = scales::percent) + 
  scale_x_discrete(name = "Total sample size") + 
  theme_bw(base_size = 18) +
  theme(legend.position="bottom")
p6

# Figure 8
ggsave(p6, device = "svg", filename = "images/Figure08.svg", width = 12, height = 8)


#  -----------------------------------------------------------------------------
#  Confidence Intervals - NOT PUBLISHED
#  -----------------------------------------------------------------------------

data.plot7 <- data.plot %>% 
  filter(              
    normalData == F,                             
    !is.na(type_ci),
    comparison == "Eta3 ~ Eta1",                     
    ssizediff == "moderately unequal", 
    groups == "2 groups") 

p7 <- ggplot(data.plot7 , aes(y = `5%_perCentage`, x = ssize, 
                              group = type_ci, color = type_ci)) +
  facet_grid(cols = vars(test), rows = vars(differences))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept=c(0.8), linetype="dotted", color = "black") +
  scale_y_continuous(name="Rejection rates", 
                     limits=c(0, 1), 
                     breaks=seq(0,1,.1),
                     labels = scales::percent) + 
  scale_x_discrete(name = "Total sample size") + 
  theme_bw(base_size = 14) +
  theme(legend.position="bottom")
p7

#  -----------------------------------------------------------------------------
#  5000 BS runs
#  -----------------------------------------------------------------------------


# 5000 Permutations ------------------------------------------------------------
data.5000 <- readRDS("data/ResultsWith5000BS.rds")


# Type I Error -----------------------------------------------------------------

data.plot8 <- data.5000  %>% 
  filter(differences == "none",                       # homogenous groups only
         test == "NPT",                               # exclude NPT
         normalData == F,                             # non-normal distributed data only
         correction == "holm" | is.na(correction),    # use holm for multiple groups
         type_ci == "CI_bca" | is.na(type_ci),        # bca
         distance == "dG" | is.na(distance),          # use geodesic distance only
         comparison == "overall",                     # overall comparision only
         ssizediff == "moderately unequal ")           # severe path difference
data.plot8 <- data.plot8 %>% rename("RR" = `5%_perCentage`)

# Add CI
data.plot8 <- data.plot8 %>% addCI(data = ., runs = data.plot8$runs[1], colum = "RR")

p8 <- ggplot(data.plot8 , aes(y = RR, x = ssize, group = test)) +
  facet_grid(cols = vars(perm))+
  geom_point()+
  geom_line() +
  geom_hline(yintercept=c(0.05), linetype="dashed", color = "black") +
  geom_line(aes(y=upperBound), linetype='dotted')+
  geom_line(aes(y=lowerBound), linetype='dotted')+
  geom_ribbon(aes(ymin = lowerBound, ymax = upperBound), alpha = 0.15)+
  scale_y_continuous(name="Rejection rates", 
                     limits=c(0, .2), 
                     breaks=seq(0,1,.05),
                     labels = scales::percent) + 
  scale_x_discrete(name = "Total sample size") +
  theme_bw(base_size = 24) +
  theme(legend.position="bottom")
p8

# Figure 5
ggsave(p8, device = "svg", filename = "images/Figure05.svg", width = 14, height = 10)


