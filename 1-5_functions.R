# Analysis and figure functions

# Load data for steroids-CMV analysis
# Event is measurement of CMV status. 
loaddat <- function(sel = "cmv"){
  d <- readRDS("Z:/18-015433_CMV in SOT_Downes/05 Data/dan_analytic_data/1_adat_2021-06-02.rds")
  d <- d %>%
    # Create steroid indicator for if any steroids administered in past 2 weeks
    group_by(transplant_id) %>%
    mutate(steroid_2wk = lag(zoo::rollapply(any_Steroids, width=14, FUN=max, 
                                            na.rm=T, fill=NA, align="right", 
                                            partial=T)),
           creat_locf = zoo::na.locf(creatinine_Value, na.rm=F),
           creat_locf_lag = lag(creat_locf),
           creat_gt1 = ifelse(creat_locf > 1, 1, 0),
           creat_gt1_lag = ifelse(creat_locf_lag > 1, 1, 0),
           test_2wk = zoo::rollapply(blood_cmv_test, width=14, FUN=sum, na.rm=T,
                                     fill=NA, align="right", partial=T)
    ) %>% 
    ungroup() %>%
    # Only look at days 14+
    filter(follow_day > 13) %>%
    
    # Set up other analytic variables
    mutate(id = as.factor(as.numeric(transplant_id)),
           id.n = as.numeric(id),
           time = follow_day - 13,
           cmv = ifelse(time == 1, 0, blood_test_result),
           steroid = any_Steroids,
           # Donor recipient status; D-/R- ref
           drany = ifelse(cmv_donor_recip_sens > 0, 1, 0), 
           # Organ transplanted; heart/kidney ref
           liver = ifelse(sot_type_collapse == "Liver", 1, 0),
           lung = ifelse(sot_type_collapse == "Lung", 1, 0),
           oimm = ifelse(other_immunosupp_count > 0, 1, 0),
           log_oimm_count = log(other_immunosupp_count + 1), # Add 1 so lowest value is transformed to 0
           creat_test = creatinine_Tested,
           creat = creatinine_Value,
           antiviral = any_antiviral
           #drstat = factor(cmv_donor_recip_sens, levels=c(0,1,2), labels=c("D-/R-", "Any D/R+", "D+/R-"))
    ) 
  
  # Restrict to cmv visits
  if(sel == "cmv"){
    d <- d %>%
      filter(blood_cmv_test == 1, 
             #!is.na(cmv_donor_recip_sens)
             ) 
    
  # Restrict to creatinine visits
  } else if(sel == "creatinine"){
    d <- d %>%
      filter(creat_test == 1, 
             blood_cmv_test == 1,
             !is.na(cmv_donor_recip_sens)) %>% # missing for 1 subject
      group_by(id) %>%
      filter(sum(creat_test) > 1) %>%
      ungroup()
  }
  d <- d %>% 
    group_by(id) %>%
    mutate(time_el = time - lag(time, default=0),
           totalobs = sum(!is.na(time)),
           obsrate = totalobs/endt) %>%
    ungroup() %>%
    select(id, id.n, time, time_el, totalobs, obsrate, blood_cmv_test, cmv, cmv_locf, 
                    steroid, steroid_2wk,
                    drany, liver, lung, 
                    age, oimm, log_oimm_count, antiviral,
                    creat_test, creat, creat_locf, creat_locf_lag, creat_gt1_lag,
                    creat_gt1, test_2wk, endt, sot_type_collapse, sot_type,
                    race_cat2, sex)
  return(d)
}


moreg <- function(data, model){
  # steroids -> cmv model 
  if(model == "stercmv"){
    # Fits gee model and pulls estimate + SE
    est <- summary(geeglm(cmv ~ steroid_2wk +
                            #drany +
                            liver + lung +
                            ns(time,3),
                          data=data,
                          id=id.n,
                          family=poisson("log")))$coefficients[,1:2]
    
    # If more than 1 obs per subject, calculate the modified GEE variance
    # estimator proposed by Mancl and DeRouen (2001).
    if(max(table(data$id.n))>1){
      nse <- try(geesmv::GEE.var.md(
        cmv ~ steroid_2wk +
          # drany +
          liver + lung +
          ns(time,3),
        data=data,
        id=id.n,
        family=poisson)$cov.beta)
      #print("SANDWICH!")
    }
    if(class(nse) == "try-error"){
      est <- NA
    } else {
      est[,2] <- nse
      est <- data.matrix(est)
    }
    return(est)
  }
}


# Store coefficients in data frame with info
storec <- function(d, m, model, rel){
  if(model == "MO"){
    d2 <- d %>% add_row(model = model, 
                        rel=rel, 
                        est = m$est[2],
                        se = m$se[2],
                        OR = exp(est),
                        LCI = exp(est - 1.96*se),
                        UCI = exp(est + 1.96*se)) 
  } else if(model == "IIW-boot"){
    d2 <- d %>% add_row(model = model, 
                        rel=rel, 
                        est = m[[1]],
                        se = NA,
                        OR = exp(est),
                        LCI = exp(quantile(m[[2]], 0.025)),
                        UCI = exp(quantile(m[[2]], 0.975)))   
  } else{
    d2 <- d %>% add_row(model = model, 
                        rel=rel, 
                        est = m$coefficients$Estimate[2],
                        se = m$coefficients$Std.err[2],
                        OR = exp(est),
                        LCI = exp(est - 1.96*se),
                        UCI = exp(est + 1.96*se)) 
  }
  return(d2)

}
