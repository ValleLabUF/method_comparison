
############################
### Mixed-Membership Sim ###
############################

library(tidyverse)

### Assess model accuracy at coarse-scale

dat<- read.csv("Modeled MM Sim Tracks w Behav.csv", as.is = T)

## Overall
dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_coarse == behav) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 95.0% to 99.7%


## For 'Resting' behavior
rest.size<- dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_coarse == 1) %>% 
  tally()

dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_coarse == 1 & behav == 1) %>% 
  tally() %>% 
  left_join(., rest.size, by = "id") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 93.5% to 99.8%

## For 'ARS' behavior
ars.size<- dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_coarse == 2) %>% 
  tally()

dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_coarse == 2 & behav == 2) %>% 
  tally() %>% 
  left_join(., ars.size, by = "id") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 95.4% to 100.0%

## For 'Transit' behavior
transit.size<- dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_coarse == 3) %>% 
  tally()

dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_coarse == 3 & behav == 3) %>% 
  tally() %>% 
  left_join(., transit.size, by = "id") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 95.0% to 99.7%





### Assess model accuracy at fine-scale

## Overall
dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_fine == behav) %>% 
  tally() %>% 
  mutate(acc = n/track_length) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 78.1% to 84.6%


## For 'Resting' behavior
rest.size2<- dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_fine == 1) %>% 
  tally()

dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_fine == 1 & behav == 1) %>% 
  tally() %>% 
  left_join(., rest.size2, by = "id") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 81.9% to 96.0%

## For 'ARS' behavior
ars.size2<- dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_fine == 2) %>% 
  tally()

dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_fine == 2 & behav == 2) %>% 
  tally() %>% 
  left_join(., ars.size2, by = "id") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 62.6% to 85.8%

## For 'Transit' behavior
transit.size2<- dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_fine == 3) %>% 
  tally()

dat %>% 
  group_by(track_length, id) %>% 
  slice(., 2:n()) %>% 
  filter(., behav_fine == 3 & behav == 3) %>% 
  tally() %>% 
  left_join(., transit.size2, by = "id") %>% 
  mutate(acc = n.x/n.y) %>% 
  summarise(min=min(acc), max=max(acc), mean=mean(acc))
# accuracy ranges from 47.2% to 84.4%


