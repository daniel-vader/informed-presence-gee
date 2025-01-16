# Pull and set up data
library(tidyverse)
library(lubridate)

data.rootdir <- "Z:/18-015433_CMV in SOT_Downes/05 Data/Derived Data/"

# Load data
d.demo <- read.csv(paste0(data.rootdir, 
                         "aim1_demos_testing_outcomes_is_20210525.csv")) %>%
  mutate(failure01 = ifelse(is.na(failure), 0, ifelse(failure == 1, 1, failure)))

# Summarize exclusions for paper
d2 <- d.demo %>% group_by(transplant_id) %>% 
  summarize(fu = max(follow_day),
            cmv_test = max(ifelse(follow_day > 13, blood_cmv_test, 0)),
            dead = max(death_during_followup),
            fail = max(d.demo$failure01),
            sot = sot_type_collapse[1],
            dr = cmv_donor_recip[1])
d3 <- d2 %>% filter(fu > 13, cmv_test==1)
d4 <- d2 %>% filter(cmv_test==0, fu > 13)

table(d2$cmv_test, d2$fu > 13) # checking follow up and any testing
table(d3$fu < 365) # Checking loss-to-follow up in eligible period
table(d4$sot) # Checking organ for patients who were dropped due to no testing
table(d4$dr) # Checking risk group for patients who were dropped due to no testing

# Include end time and failure state at every observation line
endt <- d.demo %>% 
  group_by(transplant_id) %>%
  slice_max(follow_day) %>%
  select(transplant_id, follow_day, failure) %>%
  rename(endt = follow_day, endf = failure)

d.demo <- d.demo %>% left_join(endt, by="transplant_id")

# Save data
saveRDS(d.demo, 
        paste0("Z:/18-015433_CMV in SOT_Downes/05 Data/dan_analytic_data/1_adat_", 
               Sys.Date(), ".rds"))
