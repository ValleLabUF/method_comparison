###########################
### Hard-Clustering Sim ###
###########################

### Assess model accuracy
dat<- read.csv("Modeled HC Sim Tracks w Behav.csv", as.is = T)

## Overall
true.b<- dat$true.behav[-1] %>% factor(levels = c("Resting","Exploratory","Transit")) %>%
  as.numeric()
model.b<- dat$behav[-1] %>% factor(levels = c("Resting","Exploratory","Transit")) %>%
  as.numeric()

(which(true.b == model.b) %>% length()) / length(true.b)
# 99.0% accuracy when including all different behaviors together


## For 'Resting' behavior
true.b_rest<- which(true.b == 1)
model.b_rest<- which(model.b == 1)
(which(true.b_rest %in% model.b_rest) %>% length()) / length(true.b_rest)
# 99.3% accuracy for 'Resting'


## For 'Exploratory' behavior
true.b_exp<- which(true.b == 2)
model.b_exp<- which(model.b == 2)
(which(true.b_exp %in% model.b_exp) %>% length()) / length(true.b_exp)
# 98.3% accuracy for 'Exploratory'


## For 'Transit' behavior
true.b_transit<- which(true.b == 3)
model.b_transit<- which(model.b == 3)
(which(true.b_transit %in% model.b_transit) %>% length()) / length(true.b_transit)
# 100.0% accuracy for 'Exploratory'





############################
### Mixed-Membership Sim ###
############################

### Assess model accuracy at coarse-scale

dat<- read.csv("Modeled MM Sim Tracks w Behav_multinom.csv", as.is = T)
dat$behav_fine<- gsub('Exploratory', 'ARS', dat$behav_fine)
dat$behav_coarse<- gsub('Exploratory', 'ARS', dat$behav_coarse)
dat$behav<- gsub('Exploratory', 'ARS', dat$behav)

## Overall
true.b.coarse<- dat$behav_coarse[-1] %>% factor(levels = c("Resting","ARS","Transit")) %>%
  as.numeric()
model.b<- dat$behav[-1] %>% factor(levels = c("Resting","ARS","Transit")) %>%
  as.numeric()

(which(true.b.coarse == model.b) %>% length()) / length(true.b.coarse)
# 98.9% accuracy when including all different behaviors together at coarse scale


## For 'Resting' behavior
true.b.coarse_rest<- which(true.b.coarse == 1)
model.b_rest<- which(model.b == 1)
(which(true.b.coarse_rest %in% model.b_rest) %>% length()) / length(true.b.coarse_rest)
# 99.5% accuracy for 'Resting' at coarse scale

## For 'ARS' behavior
true.b.coarse_ars<- which(true.b.coarse == 2)
model.b_ars<- which(model.b == 2)
(which(true.b.coarse_ars %in% model.b_ars) %>% length()) / length(true.b.coarse_ars)
# 97.5% accuracy for 'ARS' at coarse scale

## For 'Transit' behavior
true.b.coarse_transit<- which(true.b.coarse == 3)
model.b_transit<- which(model.b == 3)
(which(true.b.coarse_transit %in% model.b_transit) %>% length()) / length(true.b.coarse_transit)
# 99.8% accuracy for 'Transit' at coarse scale





### Assess model accuracy at fine-scale

## Overall
true.b.fine<- dat$behav_fine[-1] %>% factor(levels = c("Resting","ARS","Transit")) %>%
  as.numeric()

(which(true.b.fine == model.b) %>% length()) / length(true.b.fine)
# 83.0% accuracy when including all different behaviors together at fine scale


## For 'Resting' behavior
true.b.fine_rest<- which(true.b.fine == 1)
model.b_rest<- which(model.b == 1)
(which(true.b.fine_rest %in% model.b_rest) %>% length()) / length(true.b.fine_rest)
# 91.0% accuracy for 'Resting' at fine scale

## For 'ARS' behavior
true.b.fine_ars<- which(true.b.fine == 2)
model.b_ars<- which(model.b == 2)
(which(true.b.fine_ars %in% model.b_ars) %>% length()) / length(true.b.fine_ars)
# 81.4% accuracy for 'ARS' at fine scale

## For 'Transit' behavior
true.b.fine_transit<- which(true.b.fine == 3)
model.b_transit<- which(model.b == 3)
(which(true.b.fine_transit %in% model.b_transit) %>% length()) / length(true.b.fine_transit)
# 70.3% accuracy for 'Transit' at fine scale






####################################
### Mixed-Membership Sim 2 Behav ###
####################################

### Assess model accuracy

dat<- read.csv("Modeled MM Sim Tracks w 2 Behav.csv", as.is = T)

## Overall
true.b.coarse<- dat$behav_coarse[-1] %>% factor(levels = c("Encamped","Transit")) %>%
  as.numeric()
model.b<- dat$behav[-1] %>% factor(levels = c("Encamped","Transit")) %>%
  as.numeric()

(which(true.b.coarse == model.b) %>% length()) / length(true.b.coarse)
# 99.0% accuracy when including all different behaviors together at coarse scale



## For 'Encamped' behavior
true.b.coarse_enc<- which(true.b.coarse == 1)
model.b_enc<- which(model.b == 1)
(which(true.b.coarse_enc %in% model.b_enc) %>% length()) / length(true.b.coarse_enc)
# 99.3% accuracy for 'Encamped' at coarse scale

## For 'Transit' behavior
true.b.coarse_transit<- which(true.b.coarse == 2)
model.b_transit<- which(model.b == 2)
(which(true.b.coarse_transit %in% model.b_transit) %>% length()) / length(true.b.coarse_transit)
# 98.0% accuracy for 'Transit' at coarse scale

