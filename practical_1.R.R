## ----child="practicals/LM_ex.qmd"---------------------------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false
#| code-summary: "Load libraries"

library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     
# load some libraries to generate nice plots
library(scico)


## -----------------------------------------------------------------------------
#| code-fold: show
#| code-summary: "Simulate Data from a LM"

beta = c(2,0.5)
sd_error = 0.1

n = 100
x = rnorm(n)
y = beta[1] + beta[2] * x + rnorm(n, sd = sd_error)

df = data.frame(y = y, x = x)  











## -----------------------------------------------------------------------------
#| code-summary: "Fit LM in `inlabru`"
fit.lm = bru(cmp, lik)


## -----------------------------------------------------------------------------
#| code-summary: "Model summaries"
#| collapse: true
summary(fit.lm)


## -----------------------------------------------------------------------------
new_data = data.frame(x = c(df$x, runif(10)),
                      y = c(df$y, rep(NA,10)))
pred = predict(fit.lm, new_data, ~ beta_0 + beta_1,
               n.samples = 1000)





## -----------------------------------------------------------------------------
#| code-fold: true
#| fig-cap: Data and 95% credible intervals
#| echo: false
#| message: false
#| warning: false
#| fig-align: center
#| fig-width: 4
#| fig-height: 4

pred %>% ggplot() + 
  geom_point(aes(x,y), alpha = 0.3) +
  geom_line(aes(x,mean)) +
  geom_line(aes(x, q0.025), linetype = "dashed")+
  geom_line(aes(x, q0.975), linetype = "dashed")+
  xlab("Covariate") + ylab("Observations")



## ----child="practicals/LMM_ex.qmd"--------------------------------------------

## -----------------------------------------------------------------------------
#| code-summary: "Simulate data from a LMM"
#| 
set.seed(12)
beta = c(1.5,1)
sd_error = 1
tau_group = 1

n = 100
n.groups = 5
x = rnorm(n)
v = rnorm(n.groups, sd = tau_group^{-1/2})
y = beta[1] + beta[2] * x + rnorm(n, sd = sd_error) +
  rep(v, each = 20)

df = data.frame(y = y, x = x, j = rep(1:5, each = 20))  


## ----plot_data_lmm------------------------------------------------------------
#| code-fold: true
#| fig-cap: Data for the linear mixed model example with 5 groups
#| fig-align: center
#| fig-width: 4
#| fig-height: 4

ggplot(df) +
  geom_point(aes(x = x, colour = factor(j), y = y)) +
  theme_classic() +
  scale_colour_discrete("Group")



## ----define_components_lmm----------------------------------------------------
# Define model components
cmp =  ~ -1 + beta_0(1) + beta_1(x, model = "linear") +
  u(j, model = "iid")


## ----define_likelihood_lmm----------------------------------------------------
# Construct likelihood
lik =  bru_obs(formula = y ~.,
            family = "gaussian",
            data = df)


## -----------------------------------------------------------------------------
#| collapse: true
#| code-summary: "Fit a LMM in inlabru"
fit = bru(cmp, lik)
summary(fit)


## -----------------------------------------------------------------------------
#| code-fold: true
#| code-summary: "LMM fitted values"
#| fig-align: center
#| fig-width: 4
#| fig-height: 4

# New data
xpred = seq(range(x)[1], range(x)[2], length.out = 100)
j = 1:n.groups
pred_data = expand.grid(x = xpred, j = j)
pred = predict(fit, pred_data, formula = ~ beta_0 + beta_1 + u) 


pred %>%
  ggplot(aes(x=x,y=mean,color=factor(j)))+
  geom_line()+
  geom_ribbon(aes(x,ymin = q0.025, ymax= q0.975,fill=factor(j)), alpha = 0.5) + 
  geom_point(data=df,aes(x=x,y=y,colour=factor(j)))+
  facet_wrap(~j)




## ----child="practicals/GLM_ex.qmd"--------------------------------------------

## -----------------------------------------------------------------------------
#| code-summary: "Simulate Data from a GLM"
set.seed(123)
n = 100
beta = c(1,1)
x = rnorm(n)
lambda = exp(beta[1] + beta[2] * x)
y = rpois(n, lambda  = lambda)
df = data.frame(y = y, x = x)  






## -----------------------------------------------------------------------------
#| code-summary: "GLM likelihood"

lik =  bru_obs(formula = eta,
            family = "poisson",
            data = df)


## -----------------------------------------------------------------------------
#| code-summary: "Fit a GLM"
fit_glm = bru(cmp, lik)


## -----------------------------------------------------------------------------
#| code-summary: "GLM summaries"
summary(fit_glm)


## ----get_predictions_glm------------------------------------------------------
#| code-summary: "Predcited values for Poisson GLM"

# Define new data, set to NA the values for prediction

new_data = data.frame(x = c(df$x, runif(10)),
                      y = c(df$y, rep(NA,10)))

# Define predictor formula
pred_fml <- ~ exp(beta_0 + beta_1)

# Generate predictions
pred_glm <- predict(fit_glm, new_data, pred_fml)


## -----------------------------------------------------------------------------
#| echo: false
#| code-summary: "Plot GLM predicted values"
#| fig-cap: Data and 95% credible intervals
#| fig-align: center
#| fig-width: 4
#| fig-height: 4
#| warning: false

pred_glm %>% ggplot() + 
  geom_point(aes(x,y), alpha = 0.3) +
  geom_line(aes(x,mean)) +
    geom_ribbon(aes(x = x, ymax = q0.975, ymin = q0.025),fill = "tomato", alpha = 0.3)+
  geom_line(aes(x, q0.025), linetype = "dashed")+
  geom_line(aes(x, q0.975), linetype = "dashed")+
  xlab("Covariate") + ylab("Observations (counts)")




## -----------------------------------------------------------------------------
#| code-summary: "GLM Task"
set.seed(123)
n = 100
alpha = c(0.5,1.5)
w = rnorm(n)
psi = plogis(alpha[1] + alpha[2] * w)
y = rbinom(n = n, size = 1, prob =  psi) # set size = 1 to draw binary observations
df_logis = data.frame(y = y, w = w)  


