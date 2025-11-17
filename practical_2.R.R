## ----child="practicals/LM_ex_priors.qmd"--------------------------------------

## -----------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false

library(webexercises)



## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     


## -----------------------------------------------------------------------------
#| code-fold: show
beta = c(2,0.5)
sd_error = 0.1

n = 100
x = rnorm(n)
y = beta[1] + beta[2] * x + rnorm(n, sd = sd_error)

df = data.frame(y = y, x = x)  



## -----------------------------------------------------------------------------
# Model components
cmp =  ~ -1 + beta_0(1) + beta_1(x, model = "linear")
# Linear predictor
formula = y ~ Intercept + beta_1
# Observational model likelihood
lik =  bru_obs(formula = y ~.,
            family = "gaussian",
            data = df)
# Fit the Model
fit.lm = bru(cmp, lik)


## -----------------------------------------------------------------------------
#| eval: true
#| purl: true 

inla.priors.used(fit.lm)

#change the priors on the intercept and slope to be more informative
cmp2 =  ~-1 +  beta_0(1, prec.linear = 0.1) + beta_1(x, model = "linear",prec.linear = 0.1)

fit.lm2 = bru(cmp2, lik)

inla.priors.used(fit.lm1)
## -----------------------------------------------------------------------------

# First we define the logGamma (0.01,0.01) prior 

prec.tau <- list(prec = list(prior = "loggamma",   # prior name
                             param = c(0.01, 0.01))) # prior values

lik2 =  bru_obs(formula = y ~.,
                family = "gaussian",
                data = df,
                control.family = list(hyper = prec.tau))

fit.lm2 = bru(cmp2, lik2) 



## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-align: center
#| fig-height: 4
#| code-fold: show

plot(fit.lm, "beta_0")
plot(fit.lm2, "beta_1")

summary(fit.lm2)

plot(fit.lm2, "Precision for the Gaussian observations")


## ----child="practicals/LMM_fish_example.qmd"----------------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     



## -----------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false
#| eval: false

# library(FSAdata)
# df = PygmyWFBC
# # subset for 2001 and remove nets with zero observations
# df = df %>% filter(year==2001 & !is.na(sex) & net_no %in% c(1:5,7:8)) %>% dplyr::select(net_no,wt,tl,sex)
# write.csv(df,file = "datasets/PygmyWFBC.csv",row.names = F)
# 


## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
#| warning: false
#| message: false

PygmyWFBC <- read.csv("PygmyWFBC.csv")

ggplot(PygmyWFBC, aes(x = factor(net_no), y = wt,fill = sex)) + 
  geom_boxplot() + 
  labs(y="Weight (g)",x = "Net no.")



## -----------------------------------------------------------------------------
PygmyWFBC$sex_M <- ifelse(PygmyWFBC$sex=="F",0,1)


## -----------------------------------------------------------------------------
cmp =  ~ -1 + sex_M + 
  beta_0(1)  + 
  beta_1(tl, model = "linear") +  
  net_eff(net_no, model = "iid")

lik =  bru_obs(formula = wt ~ .,
            family = "gaussian",
            data = PygmyWFBC)

fit = bru(cmp, lik)

summary(fit)




## -----------------------------------------------------------------------------
plot(fit,"net_eff")



## -----------------------------------------------------------------------------

sampvars <-  predict(fit,PygmyWFBC, ~ {
   tau_e <- Precision_for_the_Gaussian_observations
   tau_u <- Precision_for_net_eff
   list(sigma_u = 1/tau_u,
        sigma_e = 1/tau_e)
   },
   n.samples = 1000
  )

names(sampvars) = c("Error variance","Between-net Variance")

sampvars

ICC <-  predict(fit,PygmyWFBC, ~ {
  tau_e <- Precision_for_the_Gaussian_observations
  tau_u <- Precision_for_net_eff
  sigma_u = 1/tau_u
  sigma_e = 1/tau_e
  list(ICC = sigma_u / (sigma_u + sigma_e))
       },    
n.samples = 1000
)

