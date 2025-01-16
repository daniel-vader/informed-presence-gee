# 2021-05-10 Irregular measurement with CMV in SOT data
# Execute primary analyses

################################################################################
# Load and setup data
################################################################################
library(dplyr)
library(geepack)
library(gee)
library(splines)

source("H:/projects/cmv_sot/1-5_functions.R")

################################################################################
# Baseline Table 
################################################################################
# basedat <- loaddat(sel="cmv") %>%
#   group_by(id) %>%
#   summarize(drany = drany[1])
# 
# basedat2 <- readRDS("Z:/18-015433_CMV in SOT_Downes/05 Data/dan_analytic_data/1_adat_2021-06-02.rds")

################################################################################
# Regular GEE 
################################################################################

### steroids model #############################################################
cmvdat <- loaddat(sel="cmv")


# independence model
steroid.fit <- geeglm(cmv ~ steroid_2wk + 
                    #drany + drd + 
                    liver + lung + 
                    ns(time, 3),
                  family=poisson(link="log"),
                  data=cmvdat,
                  id=id,
                  corstr="ind") %>% summary()

################################################################################
# Time elapsed since last obs 
# De Bruijne et al (2001)
################################################################################
# independence model
steroid.fit.el <- geeglm(cmv ~ steroid_2wk + 
                        #drany + drd + 
                        liver + lung + 
                        time_el +
                        ns(time, 3),
                      family=poisson(link="log"),
                      data=cmvdat,
                      id=id,
                      corstr="ind")

steroid.fit.el <- summary(steroid.fit.el)

################################################################################
# number of measurements
# author (year)
# This approach is problematic in a longitudinal context because it "looks ahead"
################################################################################
# independence model
steroid.fit.num <- geeglm(cmv ~ steroid_2wk + 
                           #drany + drd + 
                           liver + lung + 
                           obsrate +
                           ns(time, 3),
                         family=poisson(link="log"),
                         data=cmvdat,
                         id=id,
                         corstr="ind")

steroid.fit.num <- summary(steroid.fit.num)


################################################################################
# Inverse intensity weighting
# E. M. Pullenayegum tutorial: 
# https://cran.r-project.org/web/packages/IrregLong/vignettes/Irreglong-vignette.html
################################################################################

# Load suggested packages 
library(IrregLong)
library(MEMSS)
library(survival)
library(geepack)
library(data.table)
#library(GLMMadaptive)
library(geesmv)

# Inverse intensity weighting ##################################################

# Order by transplant id and day. (Should already be ordered, but running
# just in case)
cmvdat <- loaddat(sel="cmv")
cmvdat <- cmvdat[order(cmvdat$id, cmvdat$time),] %>% as.data.frame()

# Calculate inverse intenisty weights as a function of test results
iiwts.cmv <- iiw.weights(Surv(time.lag, time, blood_cmv_test)~cmv.lag + 
                       antiviral.lag + steroid.lag +
                       drany + age + sex +
                       liver + lung + cluster(id),
                     id="id",
                     time="time",
                     event="blood_cmv_test",
                     data=cmvdat,
                     invariant=c("id", "drany", "age", "sex", "liver", "lung"),
                     lagvars=c("time", "cmv",  "steroid", "log_oimm_count", 
                               "antiviral"),
                     lagfirst=c(0,0,0,0,0), 
                     maxfu=365,
                     first=F
)
# Export intensity model coefficients
im.coef <- data.frame(exp(iiwts.cmv$m$coefficients), 
                 exp(confint(iiwts.cmv$m)[,1]),
                 exp(confint(iiwts.cmv$m)[,2]))
names(im.coef) <- c("HR", "LCI", "UCI")
write.csv(im.coef, "H:/projects/cmv_sot/tables/intensity-model-2022-03-25.csv")

# Fit GEE model
iiw.fit.cmv <- iiwgee(
  # Outcome model
  cmv ~ steroid_2wk +
      # drany + 
      liver + lung + 
      ns(time,3),
                  
  # IIW model
  Surv(time.lag, time, blood_cmv_test)~cmv.lag + 
    antiviral.lag + steroid.lag +
    drany + age + sex +
    liver + lung + cluster(id),
  
  formulanull=NULL,
  family=poisson(link="log"),
  
  id="id",
  time="time",
  event="blood_cmv_test",
  invariant=c("id", "drany", "age", "sex", "liver", "lung"),
  lagvars=c("time", "cmv",  "steroid", "antiviral"),
  lagfirst=c(0,0,0,0,0), 
  maxfu=365,
  first=F,
  data=cmvdat
  )

iiw.fit.cmv.s <- iiw.fit.cmv$geefit %>% summary()

# Multiple outputation gee (uses weights from IIW) #############################
# Load modified mo function (handles NAs)
source("H:/projects/cmv_sot/1-5-1_mo-castaway.R")

# sample(1:100000, 1) # randomly select seed (6961)
set.seed(6961)
stercmv.mo <- mo.castaway(100, moreg, data=cmvdat, iiwts.cmv$iiw.weight, 
             singleobs = F, id="id.n", time="time", 
             keep.first = F, var=T, model="stercmv")


# Bootstrapped IIW GEE #########################################################
# sample(1:100000, 1) # 64267
set.seed(64267)
sid <- unique(cmvdat$id.n)
iiw.boot.fit.cmv <- NA
for(i in 1:1000){
  # Create bootstrap sample
  sid.b <- sample(sid, replace=T)
  for(k in 1:length(sid)){
    tdat <- cmvdat %>% filter(id.n == sid.b[k]) %>%
      mutate(id2 = k)
    if(k == 1){
      cmvdat2 <- tdat
    } else{
      cmvdat2 <- cmvdat2 %>% add_row(tdat)
    }
  }
  cmvdat2$id <- as.factor(cmvdat2$id2)
  # Fit model
  tfit <- iiwgee(
    # Outcome model
    cmv ~ steroid_2wk + 
      #drany + 
      liver + lung + 
      ns(time,3),
    
    # IIW model
    Surv(time.lag, time, blood_cmv_test)~cmv.lag + 
      antiviral.lag + steroid.lag +
      drany + age + sex +
      liver + lung + cluster(id),
    
    formulanull=NULL,
    family=poisson(link="log"),
    
    id="id",
    time="time",
    event="blood_cmv_test",
    invariant=c("id", "drany", "drd", "age", "liver", "lung", "cmv_2wk"),
    lagvars=c("time", "cmv",  "steroid", "log_oimm_count", 
              "antiviral"),
    lagfirst=c(0,0,0,0,0), 
    maxfu=365,
    first=F,
    data=cmvdat2
  )
  if(i == 1){
    iiw.boot.fit.cmv <- tfit$geefit$coefficients[2]
  } else {
    iiw.boot.fit.cmv <- c(iiw.boot.fit.cmv, tfit$geefit$coefficients[2])
  }
}
iiw.boot.fit.cmv <- list(main=iiw.fit.cmv$geefit$coefficients[2], bootstraps=iiw.boot.fit.cmv)




################################################################################
# Store model coefficients
################################################################################
est.eff <- tibble(model = character(), rel=character(),
                  est = numeric(), se = numeric(), 
                  OR=numeric(), LCI=numeric(), UCI=numeric())

est.eff <- storec(est.eff, steroid.fit, "GEE", "steroids -> CMV")
est.eff <- storec(est.eff, steroid.fit.el, "TSLO", "steroids -> CMV")
est.eff <- storec(est.eff, steroid.fit.num, "N Obs", "steroids -> CMV")
est.eff <- storec(est.eff, iiw.fit.cmv.s, "IIW", "steroids -> CMV")
est.eff <- storec(est.eff, iiw.boot.fit.cmv, "IIW-boot", "steroids -> CMV")
est.eff <- storec(est.eff, stercmv.mo, "MO", "steroids -> CMV")

saveRDS(est.eff, paste0("H:/projects/cmv_sot/output/effects_IRR_", Sys.Date(), ".rds"))

t <- readRDS("H:/projects/cmv_sot/output/effects_IRR_2022-03-28.rds")
write.csv(t, "H:/projects/cmv_sot/tables/effects_IRR_2022-03-28.csv")


################################################################################
# Sensitivity analysis: Exclude Low Risk patients?
################################################################################