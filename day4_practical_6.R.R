## ----child="practicals/space_time_ex.qmd"-------------------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(inlabru) 
library(sf)
library(terra)


# load some libraries to generate nice map plots
library(scico)
library(ggplot2)
library(patchwork)
library(mapview)
library(tidyterra)



## -----------------------------------------------------------------------------
#| message: false
#| warning: false
library(sdmTMB)

pcod_df = sdmTMB::pcod  %>% filter(year<=2005)
qcs_grid = sdmTMB::qcs_grid



## -----------------------------------------------------------------------------
#| message: false
#| warning: false
#| 
pcod_sf =   st_as_sf(pcod_df, coords = c("lon","lat"), crs = 4326)
pcod_sf = st_transform(pcod_sf,
          crs = "+proj=utm +zone=9 +datum=WGS84 +no_defs +type=crs +units=km" )


## -----------------------------------------------------------------------------

depth_r <- rast(qcs_grid, type = "xyz")
crs(depth_r) <- crs(pcod_sf)

bru_options_set(control.compute = list(waic = TRUE,dic= TRUE,mlik = TRUE))





## -----------------------------------------------------------------------------
mesh = fm_mesh_2d(loc = pcod_sf,    
                  cutoff = 1,
                  max.edge = c(10,20),     
                  offset = c(5,50),
                  crs = st_crs(pcod_df))   

spde_model =  inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(100, 0.5))


ggplot()+gg(mesh)+ geom_sf(data= pcod_sf, aes(color = factor(present)), size = 0.1)

  

## -----------------------------------------------------------------------------
depth_r$depth_group = inla.group(values(depth_r$depth_scaled))



## -----------------------------------------------------------------------------
pcod_sf = pcod_sf %>%
     mutate(time_idx = match(year, c(2003, 2004, 2005)),
         id = 1:nrow(.)) # Observation id for CV




## -----------------------------------------------------------------------------
# Select the raster of interest
depth_orig = depth_r$depth_group
re <- extend(depth_orig, ext(depth_orig)*1.05)
# Convert to an sf spatial object
re_df <- re %>% stars::st_as_stars() %>%  st_as_sf(na.rm=F)
# fill in missing values using the original raster 
re_df$depth_group =  bru_fill_missing(depth_orig,re_df,re_df$depth_group)
# rasterize
depth_filled <- stars::st_rasterize(re_df) %>% rast()
plot(depth_filled)


cmp_spat=~Intercept(1)+ covariate(depth_filled, model="rw2",scale.model=TRUE)+
  trend(time_idx, model="iid") +
  space(geometry, model=spde_model)

formula_spat=present~Intercept + covariate + trend + space

lik_spat<-bru_obs(formula=formula_spat,
                     family="binomial",
                     data=pcod_sf)
# Fit Model 
fit_spat = bru(cmp_spat,lik_spat)

summary(fit_spat)

## -----------------------------------------------------------------------------
pxl = st_as_sf(data.frame(crds(depth_r)), coords = c("x","y") ,
               crs  = st_crs(pcod_sf))


p_cov=fit_spat$summary.random$covariate%>%
  ggplot()+geom_ribbon(aes(ID,ymin = `0.025quant`, ymax= `0.975quant`), alpha = 0.5, fill="lightblue")+
  geom_line(aes(ID,mean), color="blue")

p_time=fit_spat$summary.random$trend%>%
  ggplot()+geom_errorbar(aes(ID,ymin = `0.025quant`, ymax= `0.975quant`), alpha = 0.5, fill="lightgreen")+
  geom_point(aes(ID,mean))

#space effect

pred_space= predict(fit_spat, pxl,
                       formula=~ space)

p_space_mean<-ggplot()+gg(pred_space,aes(color = mean))
p_space_sd<-ggplot()+gg(pred_space,aes(color = sd))

## -----------------------------------------------------------------------------
# PC prior for AR(1) correlation parameter
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.1)))

cmp_st=~ Intercept(1)+
  covariate(depth_filled, model="rw2",scale.model=TRUE)+
  space_time(geometry, 
             group=time_idx,
             model=spde_model,
             control.group=list(model="ar1", hyper = h.spec))

formula_st=present~Intercept + covariate + space_time
lik_st<-bru_obs(formula=formula_st,
                     family="binomial",
                     data=pcod_sf)
# Fit Model
fit_st = bru(cmp_st,lik_st)

summary(fit_st)
fit_st$summary.hyperpar

fit_st$summary.random$covariate%>%
  ggplot()+geom_ribbon(aes(ID,ymin = `0.025quant`, ymax= `0.975quant`), alpha = 0.5, fill="lightblue")+
  geom_line(aes(ID,mean), color="blue")

inv_logit = function(x){ exp(x) / (1 + exp(x))}

pxl_all = fm_cprod(pxl, data.frame(time_idx = 1:3))

pred_space_time=predict(fit_st,pxl_all, ~data.frame(
  logit_prob=Intercept + covariate + space_time,
  prob= inv_logit(Intercept + covariate + space_time)))

ggplot(pred_space_time$prob)+geom_sf(aes(color = mean))+
  facet_wrap(~time_idx)+ scale_color_scico(direction = -1) +
  cowplot::theme_map()


out= data.frame(Model = c("Model 1", "Model 2"),
                DIC = c(fit_spat$dic$dic, fit_st$dic$dic),
                WAIC = c(fit_spat$waic$waic, fit_st$waic$waic),
                MLIK = c(fit_spat$mlik[1], fit_st$mlik[1]))
## ----child="practicals/LM_resids.qmd"-----------------------------------------

## -----------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false

library(webexercises)



## -----------------------------------------------------------------------------
#| warning: false
#| message: false

library(dplyr)
library(tidyr)
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
res_samples <- predict(
  fit.lm,         # the fitted model
  df,             # the original data set
  ~ data.frame(   
    res = y-(beta_0 + beta_1)  # compute the residuals
  ),
  n.samples = 1000   # draw 1000 samples
)



## -----------------------------------------------------------------------------
#| fig-width: 6
#| fig-height: 4
#| fig-align: center
#| fig-cap: "Bayesian residual plots: the left panel is the residual index plot; the right panel is the plot of the residual versus the covariate x"
#| code-fold: true
#| code-summary: "Residuals checks for Linear Model"

ggplot(res_samples,aes(y=mean,x=1:100))+geom_point() +
ggplot(res_samples,aes(y=mean,x=x))+geom_point()



## -----------------------------------------------------------------------------
#| fig-width: 4
#| fig-height: 4
#| fig-align: center
#| code-fold: true
#| code-summary: "QQPlot for Linear Model"

arrange(res_samples, mean) %>%
  mutate(theortical_quantiles = qnorm(1:100 / (1+100))) %>%
  ggplot(aes(x=theortical_quantiles,y= mean)) + 
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "grey70")+
  geom_abline(intercept = mean(res_samples$mean),
              slope = sd(res_samples$mean)) +
  geom_point() +
  labs(x = "Theoretical Quantiles (Normal)",
       y= "Sample Quantiles (Residuals)") 


## -----------------------------------------------------------------------------
#| warning: false
#| message: false

samples =  generate(fit.lm, df,
  formula = ~ {
    mu <- (beta_0 + beta_1)
    sd <- sqrt(1 / Precision_for_the_Gaussian_observations)
    rnorm(100, mean = mu, sd = sd)
  },
  n.samples = 500
) 


## -----------------------------------------------------------------------------
# Tidy format for plotting
samples_long = data.frame(samples) %>% 
  mutate(id = 1:100) %>% # i-th observation
  pivot_longer(-id)

# compute the mean and quantiles for the samples
draws_summaries = data.frame(mean_samples = apply(samples,1,mean),
q25 = apply(samples,1,function(x)quantile(x,0.025)),  
q975 = apply(samples,1,function(x)quantile(x,0.975)),
observations = df$y)  

p1 = ggplot() + geom_density(data = samples_long, 
                        aes(value, group = name),  color = "#E69F00") +
  geom_density(data = df, aes(y))  +
  xlab("") + ylab("") 

p2 = ggplot(draws_summaries,aes(y=mean_samples,x=observations))+
  geom_errorbar(aes(ymin = q25,
                   ymax = q975), 
               alpha = 0.5, color = "grey50")+
geom_point()+geom_abline(slope = 1,intercept = 0,lty=2)+labs()

p1 +p2



## ----child="practicals/GLM_scores_ex.qmd"-------------------------------------

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
crabs <- read.csv("datasets/crabs.csv")

# conditional means and variances
crabs %>%
  summarise( Mean = mean(satell ),
             Variance = var(satell),
                     .by = color)


## -----------------------------------------------------------------------------

crabs_df = model.matrix( ~  color , crabs) %>%
  as.data.frame() %>%
  select(-1) %>%        # drop intercept
  bind_cols(crabs) %>%  # append to original data
  select(-color)        # remove original color categorical variable



## -----------------------------------------------------------------------------

cmp =  ~ -1 + beta0(1) +  colordarker +
       colorlight + colormedium +
       w(weight, model = "linear")

lik =  bru_obs(formula = satell ~.,
            family = "poisson",
            data = crabs_df)

fit_pois = bru(cmp, lik)

summary(fit_pois)



## -----------------------------------------------------------------------------
bru_options_set(control.compute = list(cpo = TRUE))

fit_pois = bru(cmp, lik)




## -----------------------------------------------------------------------------
#| eval: false

# 
# fit_pois$cpo$pit %>%
#   hist(main = "Histogram of PIT values")
# 
# qqplot(qunif(ppoints(length(fit_pois$cpo$pit))),
#        fit_pois$cpo$pit,
#        main = "Q-Q plot for Unif(0,1)",
#        xlab = "Theoretical Quantiles",
#        ylab = "Sample Quantiles")
# 
# qqline(fit_pois$cpo$pit,
#        distribution = function(p) qunif(p),
#        prob = c(0.1, 0.9))


