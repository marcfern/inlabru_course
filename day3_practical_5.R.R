## ----child="practicals/Areal_ex.qmd"------------------------------------------

## -----------------------------------------------------------------------------
#| echo: false

# load webexercises library for tasks and questions (just for a preview - the practical compiler should take care of this when compiling multiple excercises)
library(webexercises)



## -----------------------------------------------------------------------------
library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)   
library(mapview)

# load some libraries to generate nice map plots
library(scico)


## -----------------------------------------------------------------------------

library(CARBayesdata)

data(pollutionhealthdata)
data(GGHB.IZ)




## -----------------------------------------------------------------------------

resp_cases <- merge(GGHB.IZ %>%
                      mutate(space = 1:dim(GGHB.IZ)[1]),
                             pollutionhealthdata, by = "IZ") %>%
  filter(year == 2007) %>%
    mutate(SMR = observed/expected)

ggplot() + geom_sf(data = resp_cases, aes(fill = SMR)) + scale_fill_scico(direction = -1)




## -----------------------------------------------------------------------------
library(spdep)

W.nb <- poly2nb(GGHB.IZ,queen = TRUE)
R <- nb2mat(W.nb, style = "B", zero.policy = TRUE)

diag = apply(R,1,sum)
Q = -R
diag(Q) = diag


## -----------------------------------------------------------------------------
#| echo: true
#| eval: false

# cmp = ~ Intercept(1) + space(...) + iid(...)
# 
# formula = ...
# 
# 
# lik = bru_obs(formula = formula,
#               family = ...,
#               E = ...,
#               data = ...)
# 
# fit = bru(cmp, lik)
# 


## -----------------------------------------------------------------------------

cmp = ~ Intercept(1) + space(space, model = "besag", graph = Q) + iid(space, model = "iid")

formula = observed ~ Intercept + space + iid

lik = bru_obs(formula = formula, 
              family = "poisson",
              E = expected,
              data = resp_cases)

fit = bru(cmp, lik)



## -----------------------------------------------------------------------------
#| echo: true
#| eval: false

# pred = predict(fit, resp_cases, ~data.frame(log_risk = ...,
#                                              risk = exp(...),
#                                              cases = ...
#                                              ),
#                n.samples = 1000)
# 


## -----------------------------------------------------------------------------

# produce predictions
pred = predict(fit, resp_cases, ~data.frame(log_risk = Intercept + space,
                                             risk = exp(Intercept + space),
                                             cases = expected * exp(Intercept + space)
                                             ),
               n.samples = 1000)
# plot the predictions

p1 = ggplot() + geom_sf(data = pred$log_risk, aes(fill = mean)) + scale_fill_scico(direction = -1) + ggtitle("mean log risk")
p2 = ggplot() + geom_sf(data = pred$log_risk, aes(fill = sd)) + scale_fill_scico(direction = -1) + ggtitle("sd log risk")
p1 + p2

p1 = ggplot() + geom_sf(data = pred$risk, aes(fill = mean)) + scale_fill_scico(direction = -1) + ggtitle("mean  risk")
p2 = ggplot() + geom_sf(data = pred$risk, aes(fill = sd)) + scale_fill_scico(direction = -1) + ggtitle("sd  risk")
p1 + p2

p1 = ggplot() + geom_sf(data = pred$cases, aes(fill = mean)) + scale_fill_scico(direction = -1)+ ggtitle("mean  expected counts")
p2 = ggplot() + geom_sf(data = pred$cases, aes(fill = sd)) + scale_fill_scico(direction = -1)+ ggtitle("sd  expected counts")
p1 + p2



## -----------------------------------------------------------------------------
pred$cases %>% ggplot() + geom_point(aes(observed, mean)) + 
  geom_errorbar(aes(observed, ymin = q0.025, ymax = q0.975)) +
  geom_abline(intercept = 0, slope = 1)



## -----------------------------------------------------------------------------

# simulate 1000 realizations of E_i\lambda_i
expected_counts = generate(fit, resp_cases, 
                           ~ expected * exp(Intercept + space),
                           n.samples = 1000)


# simulate poisson data
aa = rpois(271*1000, lambda = as.vector(expected_counts))
sim_counts = matrix(aa, 271, 1000)

# summarise the samples with posterior means and quantiles
pred_counts = data.frame(observed = resp_cases$observed,
                         m = apply(sim_counts,1,mean),
                         q1 = apply(sim_counts,1,quantile, 0.025),
                         q2 = apply(sim_counts,1,quantile, 0.975),
                         vv = apply(sim_counts,1,var)
                         )



## ----child="practicals/Geostat_ex.qmd"----------------------------------------

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

pcod_df = sdmTMB::pcod %>% filter(year==2003)
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



## -----------------------------------------------------------------------------

ggplot()+ 
  geom_spatraster(data=depth_r$depth)+
  geom_sf(data=pcod_sf,aes(color=factor(present))) +
    scale_color_manual(name="Occupancy status for the Pacific Cod",
                     values = c("black","orange"),
                     labels= c("Absence","Presence"))+
  scale_fill_scico(name = "Depth",
                   palette = "nuuk",
                   na.value = "transparent" ) + xlab("") + ylab("")



## -----------------------------------------------------------------------------

mesh = fm_mesh_2d(loc = pcod_sf,           # Build the mesh
                  max.edge = c(10,20),     # The largest allowed triangle edge length.
                  offset = c(5,50))       # The automatic extension distance
ggplot() + gg(mesh) + geom_sf(data= pcod_sf, aes(color = factor(present)), size = 0.1) + xlab("") + ylab("")




## -----------------------------------------------------------------------------

mesh = fm_mesh_2d(loc = pcod_sf,           # Build the mesh
                  cutoff = 2,
                  max.edge = c(10,20),     # The largest allowed triangle edge length.
                  offset = c(5,50))       # The automatic extension distance
ggplot() + gg(mesh) + geom_sf(data= pcod_sf, aes(color = factor(present)), size = 0.1) + xlab("") + ylab("")



## -----------------------------------------------------------------------------
#| eval: false

# ?fm_mesh_2d


## -----------------------------------------------------------------------------
dens_prior_range = function(rho_0, p_alpha)
{
  # compute the density of the PC prior for the
  # range rho of the Matern field
  # rho_0 and p_alpha are defined such that
  # P(rho<rho_0) = p_alpha
  rho = seq(0, rho_0*10, length.out =100)
  alpha1_tilde = -log(p_alpha) * rho_0
  dens_rho =  alpha1_tilde / rho^2 * exp(-alpha1_tilde / rho)
  return(data.frame(x = rho, y = dens_rho))
}

dens_prior_sd = function(sigma_0, p_sigma)
{
  # compute the density of the PC prior for the
  # sd sigma of the Matern field
  # sigma_0 and p_sigma are defined such that
  # P(sigma>sigma_0) = p_sigma
  sigma = seq(0, sigma_0*10, length.out =100)
  alpha2_tilde = -log(p_sigma)/sigma_0
  dens_sigma = alpha2_tilde* exp(-alpha2_tilde * sigma) 
  return(data.frame(x = sigma, y = dens_sigma))
}


## -----------------------------------------------------------------------------
spde_model1 =  inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(.1, 0.5),
                                  prior.range = c(30, 0.5))
spde_model2 =  inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(10, 0.5),
                                  prior.range = c(1000, 0.5))
spde_model3 =  inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(100, 0.5))


## -----------------------------------------------------------------------------
ggplot() + 
  geom_line(data = dens_prior_range(30,.5), aes(x,y, color = "model1")) +
  geom_line(data = dens_prior_range(1000,.5), aes(x,y, color = "model2")) +
  geom_line(data = dens_prior_range(100,.5), aes(x,y, color = "model3")) 


## -----------------------------------------------------------------------------
ggplot() + 
  geom_line(data = dens_prior_sd(1,.5), aes(x,y, color = "model1")) +
  geom_line(data = dens_prior_sd(10,.5), aes(x,y, color = "model2")) +
  geom_line(data = dens_prior_sd(.1,.5), aes(x,y, color = "model3")) 


## -----------------------------------------------------------------------------

cmp = ~ Intercept(1) + space(geometry, model = spde_model3)


## -----------------------------------------------------------------------------

formula = present ~ Intercept  + space

lik = bru_obs(formula = formula, 
              data = pcod_sf, 
              family = "binomial")



## -----------------------------------------------------------------------------
fit1 = bru(cmp,lik)






## -----------------------------------------------------------------------------
pxl = fm_pixels(mesh)


## -----------------------------------------------------------------------------
preds = predict(fit1, pxl, ~data.frame(spatial = space,
                                      total = Intercept + space))


## -----------------------------------------------------------------------------
ggplot() + geom_sf(data = preds$spatial,aes(color = mean)) + scale_color_scico() + ggtitle("Posterior mean")

ggplot() + geom_sf(data = preds$spatial,aes(color = sd)) + scale_color_scico() + ggtitle("Posterior sd")


## -----------------------------------------------------------------------------
pxl1 = data.frame(crds(depth_r), 
                  as.data.frame(depth_r$depth)) %>% 
       filter(!is.na(depth)) %>%
st_as_sf(coords = c("x","y"))




## -----------------------------------------------------------------------------
# we simulate 4 samples from the 
gens = generate(fit1, pxl1, ~ (Intercept + space),
                n.samples = 4)

pp = cbind(pxl1, gens)

pp %>% select(-depth) %>%
  pivot_longer(-geometry) %>%
    ggplot() + 
      geom_sf(aes(color = value)) +
      facet_wrap(.~name) +
        scale_color_scico(direction = -1) +
        ggtitle("Sample from the fitted model")




## -----------------------------------------------------------------------------

# create the grouped variable
depth_r$depth_group = inla.group(values(depth_r$depth_scaled))

# run the model
cmp = ~ Intercept(1) + space(geometry, model = spde_model3) +
        covariate(depth_r$depth_group, model = "rw2")

formula = present ~ Intercept  + space + covariate

lik = bru_obs(formula = formula, 
              data = pcod_sf, 
              family = "binomial")


fit3 = bru(cmp, lik)

# plot the estimated effect of depth

fit3$summary.random$covariate %>% 
  ggplot() + geom_line(aes(ID,mean)) + 
                                  geom_ribbon(aes(ID, ymin = `0.025quant`, 
                                                      ymax = `0.975quant`), alpha = 0.5)


## -----------------------------------------------------------------------------
#| eval: false

# inv_logit = function(x) (1+exp(-x))^(-1)



## ----child="practicals/PointProcess_ex.qmd"-----------------------------------

## -----------------------------------------------------------------------------
#| echo: false
#| message: false
#| warning: false

# load webexercises library for tasks and questions (just for a preview - the practical compiler should take care of this when compiling multiple excercises)
library(webexercises)



## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     
library(spatstat)
library(sf)
library(scico)
library(spatstat)
library(lubridate)
library(terra)
library(tidyterra)




## -----------------------------------------------------------------------------
#| label: fig-points
#| fig-cap: "Distribution of the observed forest fires caused by lightning in Castilla-La Mancha in 2004"
#| 
data("clmfires")
pp = st_as_sf(as.data.frame(clmfires) %>%
                mutate(x = x, 
                       y = y),
              coords = c("x","y"),
              crs = NA) %>%
  filter(cause == "lightning",
         year(date) == 2004)

poly = as.data.frame(clmfires$window$bdry[[1]]) %>%
  mutate(ID = 1)

region = poly %>% 
  st_as_sf(coords = c("x", "y"), crs = NA) %>% 
  dplyr::group_by(ID) %>% 
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 
  
ggplot() + geom_sf(data = region, alpha = 0) + geom_sf(data = pp)  


## -----------------------------------------------------------------------------
# define integration scheme

ips = st_sf(
geometry = st_sample(region, 1)) # some random location inside the domain
ips$weight = st_area(region) # integration weight is the area of the domain

cmp = ~ 0 + beta_0(1)

formula = geometry ~ beta_0

lik = bru_obs(data = pp,
              family = "cp",
              formula = formula,
              ips = ips)
fit1 = bru(cmp, lik)




## -----------------------------------------------------------------------------
#|label: fig-altitude
#|fig-cap: "Distribution of the observed forest fires and scaled altitude"
#| 
elev_raster = rast(clmfires.extra[[2]]$elevation)
elev_raster = scale(elev_raster)
ggplot() + 
  geom_spatraster(data = elev_raster) + 
  geom_sf(data = pp) +
  scale_fill_scico()



## -----------------------------------------------------------------------------
#|label: fig-int2
#|fig-cap: "Integration scheme."

n.int = 1000
ips = st_sf(geometry = st_sample(region,
            size = n.int,
            type = "regular"))

ips$weight = st_area(region) / n.int
ggplot() + geom_sf(data = ips, aes(color = weight)) + geom_sf(data= region, alpha = 0)



## -----------------------------------------------------------------------------
cmp = ~ Intercept(1) + elev(elev_raster, model = "linear")
formula = geometry ~ Intercept + elev
lik = bru_obs(data = pp,
              family = "cp",
              formula = formula,
              ips = ips)
fit2 = bru(cmp, lik)




## -----------------------------------------------------------------------------
n.int2 = 50

ips2 = st_sf(geometry = st_sample(region,
            size = n.int2,
            type = "regular"))
ips2$weight = st_area(region) / n.int2




## -----------------------------------------------------------------------------
est_grid = st_as_sf(data.frame(crds(elev_raster)), coords = c("x","y"))
est_grid  = st_intersection(est_grid, region)




## -----------------------------------------------------------------------------
N_fires = generate(fit2, ips,
                      formula = ~ {
                        lambda = sum(weight * exp(elev + Intercept))
                        rpois(1, lambda)},
                    n.samples = 2000)

ggplot(data = data.frame(N = as.vector(N_fires))) +
  geom_histogram(aes(x = N),
                 colour = "blue",
                 alpha = 0.5,
                 bins = 20) +
  geom_vline(xintercept = nrow(pp),
             colour = "red") +
  theme_minimal() +
  xlab(expression(Lambda))



## -----------------------------------------------------------------------------
mesh = fm_mesh_2d(boundary = region,
                  max.edge = c(5, 10),
                  cutoff = 4, crs = NA)

ggplot() + gg(mesh) + geom_sf(data = pp)

spde_model =  inla.spde2.pcmatern(mesh,
                                  prior.sigma = c(1, 0.5),
                                  prior.range = c(100, 0.5))



## -----------------------------------------------------------------------------
ips = fm_int(mesh, samplers = region)

ggplot() + geom_sf(data = ips, aes(color = weight)) +
  gg(mesh) +
   scale_color_scico()


## -----------------------------------------------------------------------------
#| echo: true
#| eval: false

# 
# cmp = ~ ...
# 
# formula = geometry ~ ...
# 
# lik = bru_obs("cp",
#               formula = formula,
#               data = pp,
#               ips = ...)
# 
# fit3 = bru(cmp, lik)
# 




## -----------------------------------------------------------------------------
#| echo: false
#| eval: true
#| warning: true

fit3 = bru(cmp, lik)





## -----------------------------------------------------------------------------
#| fig-width: 4.5
#| fig-height: 4.5
#| fig-align: center
# Extend raster ext by 30 % of the original raster so it covers the whole mesh
re <- extend(elev_raster, ext(elev_raster)*1.3)
# Convert to an sf spatial object
re_df <- re %>% stars::st_as_stars() %>%  st_as_sf(na.rm=F)
# fill in missing values using the original raster 
re_df$lyr.1 <- bru_fill_missing(elev_raster,re_df,re_df$lyr.1)
# rasterize
elev_rast_p <- stars::st_rasterize(re_df) %>% rast()
ggplot() + geom_spatraster(data = elev_rast_p) 





## -----------------------------------------------------------------------------
sim_fields = generate(fit3, pxl, ~data.frame(spde = space,
                                       log_int = Intercept + space + elev),
                     n.samples = 4)

cbind(pxl,sapply(sim_fields, function(x) x$spde)) %>%
  pivot_longer(-geometry) %>%
  ggplot() + geom_sf(aes(color = value)) + 
  facet_wrap(.~name) + scale_color_scico() +
  ggtitle("simulated spatial fields")


cbind(pxl,sapply(sim_fields, function(x) x$log_int)) %>%
  pivot_longer(-geometry) %>%
  ggplot() + geom_sf(aes(color = value)) + 
  facet_wrap(.~name) + scale_color_scico() + 
  ggtitle("simulated log intensity")




