############################################################################################################
### R codes for "Floods and cause-specific mortality in the United States: A triply robust approach"
### Dr. Lingzhi Chu, Yale School of Public Health
### 6/11/2024
############################################################################################################


require(fixest); require(splines); library(dplyr)

### step 1: calculate the propensity scores
ps_model = feglm(treat.status ~ flood.ma.surrounding + prcp.surrounding + ns(temp.ma, df=5) + prcp.ma + snow.ma | county + year^month + state^year, data = model.data, family="binomial", cluster = "county", combine.quick = FALSE) # feglm considers 1 as success and predicts the probability of being 1.
model.data$ps <- stats::predict(ps_model, newdata=model.data)
model.data <- model.data %>% 
  mutate(iptw = case_when(treat.status==1 ~ 1/ps,
                          treat.status==0 ~ 1/(1-ps)))
# model data stored mortality rates, exposure data, covariates, and treatment/control status
# For a specific county, flood.ma.surrounding and prcp.surrounding are the moving averages of monthly flood days and monthly precipitation over the time window of interest, averaged over the other counties in the same state.
# temp.ma, prcp.ma, and snow.ma are the moving averages of county-specific monthly temperature, precipitation, and snow.



### step 2: counterfactual estimation
## (1) modeling in the control group
# [Inverse Probability of Treatment Weighting Using the Propensity Score]
control.data <- model.data %>% 
  filter(treat.status==0)
model.step2 = feglm(death.count ~ flood.ma.surrounding + prcp.surrounding + ns(temp.ma, df=5) + prcp.ma + snow.ma | county + year^month + state^year, data = control.data, family="quasipoisson", offset = ~log(pop), weights = control.data$iptw, cluster = "county", combine.quick = FALSE)
# death counts are age-standardized observed counts.

## (2) predict counterfactual outcomes for the treatment group
treat.data <- model.data %>% 
  filter(treat.status==1)
treat.data$death.count.pred <- stats::predict(model.step2, newdata = treat.data)



### step 3: flood association estimation in the treatment group (confounder adjustment)
# [Inverse Probability of Treatment Weighting Using the Propensity Score]
treat.data$residual <- treat.data$death.count - treat.data$death.count.pred
# residual counts not explained = observed counts minus predicted counts
model.step3 <- lm(residual/pop~flood.ma.county + flood.ma.surrounding + prcp.surrounding + ns(temp.ma, df=5) + prcp.ma + snow.ma, data=treat.data, weights = iptw, na.action="na.exclude")
model.summary <- summary(model.step3)
result <- model.summary$coefficients
