## ----child="practicals/Multiple_likelihood.qmd"-------------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(inlabru) 
library(sf)
library(terra)
library(tidyverse)

# load some libraries to generate nice map plots
library(scico)
library(ggplot2)
library(patchwork)



## -----------------------------------------------------------------------------

N = 200
x =  runif(N)
df = data.frame(idx = 1:N,
                x = x)

# simulate data
df = df %>% 
  mutate(y_gaus = rnorm(N, mean = 1 + 1.5 * x), sd = 0.5) %>%
  mutate(y_pois = rpois(N, lambda  = exp( -1 + 1.5 * x))) 

# plot the data
df %>% ggplot() + 
  geom_point(aes(x, y_gaus, color = "Gaussian")) +
  geom_point(aes(x, y_pois, color = "Poisson")) 



## -----------------------------------------------------------------------------
cmp = ~ -1 + 
  Intercept_gaus(1) + 
  Intercept_pois(1) +
  covariate(x, model = "linear") 


## -----------------------------------------------------------------------------
lik_gaus = bru_obs(formula = y_gaus ~ Intercept_gaus + covariate,
                    data = df)

lik_pois = bru_obs(formula = y_pois ~ Intercept_pois + covariate,
                    data = df,
                   family = "poisson")


fit_gaus<-bru(cmp,
         lik_gaus)

fit_pois<-bru(cmp,
         lik_pois)


fit_join<-bru(cmp,
         lik_gaus,
         lik_pois)

ggplot()+
  geom_line(data=fit_gaus$marginals.fixed$covariate,
                   aes(x = x, y = y,color = "gaussian"))+
  geom_line(data=fit_pois$marginals.fixed$covariate,
            aes(x = x, y = y, color = "poisson"))+
  geom_line(data=fit_join$marginals.fixed$covariate,
            aes(x = x, y = y, color = "joint"))

summary(fit)

## -----------------------------------------------------------------------------
rMatern <- function(n, coords, sigma=1, range, 
                    kappa = sqrt(8*nu)/range, 
                    variance = sigma^2, 
                    nu=1) {
  m <- as.matrix(dist(coords))
  m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  diag(m) <- 1
  return(drop(crossprod(chol(variance*m),
                        matrix(rnorm(nrow(coords)*n), ncol=n))))
}


## -----------------------------------------------------------------------------
# Intercept on reparametrized model
beta <- c(-5, 3) 
# Random field marginal variances for omega1 and omega2:
m.var <- c(0.5, 0.4) 
# GRF range parameters for omega1 and omega2:
range <- c(4, 6)
# Copy parameters: reparameterization of coregionalization 
# parameters
lambda <- c(0.7) 
# Standard deviations of error terms
e.sd <- c(0.3, 0.2)



## -----------------------------------------------------------------------------
# define the area of interest
poly_geom = st_polygon(list(cbind(c(0,10,10,0,0), c(0,0,5,5,0)) ))
# Wrap it in an sfc (simple feature collection)
poly_sfc <- st_sfc(poly_geom)
# Now create the sf object
border <- st_sf(id = 1, geometry = poly_sfc)



# how many observation we have
n1 <- 200
n2 <- 150
n_common = 50

# simulate observation locations

loc_common = st_sf(geometry = st_sample(border, n_common))
loc_only1 = st_sf(geometry = st_sample(border, n1-n_common))
loc_only2 = st_sf(geometry = st_sample(border, n2-n_common))



# simulate the two gaussian field at the locations
z1 <- rMatern(1, st_coordinates(rbind( loc_common,loc_only1, loc_only2)), range = range[1],
                  sigma = sqrt(m.var[1]))

z2 <- rMatern(1, st_coordinates(rbind(loc_common, loc_only2)), range = range[2],
                  sigma = sqrt(m.var[2]))


## Create data.frame
loc1 = rbind( loc_common, loc_only1)
loc2 = rbind( loc_common, loc_only2)

df1 =  loc1 %>% mutate(z1 = z1[1:n1])
df2 =  loc2 %>% mutate(z1 = z1[-c(1:(n1-n_common))], z2 =z2)


## create the linear predictors

df1  = df1 %>%
  mutate(eta1 = beta[1] + z1)

df2  = df2 %>%
  mutate(eta2 = beta[2] + lambda * z1 + z2)


# simulate data by addint the obervation noise

df1  = df1 %>%
  mutate(y = rnorm(n1, mean = eta1, sd = e.sd[1]))

df2  = df2 %>%
  mutate(y = rnorm(n2, mean = eta2, sd = e.sd[1]))


## ----out.width="95%"----------------------------------------------------------
p1 = ggplot(data = df1) + geom_sf(aes(color = z1)) 
p2 = ggplot(data = df2) + geom_sf(aes(color = z2)) 
p1+p2+plot_layout(ncol = 1)


## -----------------------------------------------------------------------------
mesh <-  fm_mesh_2d(loc = rbind(loc1, loc2), 
                   boundary = border,
                     max.edge = c(0.5, 1.5), 
                     offset = c(0.1, 2.5), 
                     cutoff = 0.1)



## ----echo = F-----------------------------------------------------------------
ggplot() + 
  gg(mesh) + 
  geom_sf(data =  df1,  size = 2, aes(color = "data 1")) +
    geom_sf(data =  df2, aes(color = "data 2")) + xlab("") + ylab("")


spde<-inla.spde2.pcmatern(mesh,
                prior.range = c(0.5,0.01),
                prior.sigma = c(1,0.01))

## -----------------------------------------------------------------------------
cmp = ~ -1 +  Intercept1(1) + Intercept2(1) +
  omega1(geometry, model = spde) +
  omega1_copy(geometry, copy = "omega1", fixed = FALSE) +
  omega2(geometry, model = spde)

lik1 = bru_obs(formula = y ~ Intercept1 + omega1,
               family = "gaussian",
                data = df1)
lik2 = bru_obs(formula = y ~ Intercept2 + omega1_copy + omega2,
               family = "gaussian",
                data = df2)
res = bru(cmp,lik1,lik2)

summary(fit)

#fixed effects
fixed = data.frame(true = beta, res$summary.fixed[,c(1,3,5)])

#hyperparameters
hyper = data.frame(true = c(1/e.sd^2, range[1], sqrt(m.var[1]),
                            range[2], sqrt(m.var[2]),
                            lambda),
                   res$summary.hyperpar[,c(1,3,5)])


#Compute predictions from the model at the observation points and compare them with the observed values.

pred1=predict(res, df1,formula= ~Intercept1+ omega1)
pred2=predict(res, df2,formula= ~Intercept2+ omega1_copy + omega2)

ggplot()+ geom_point(data = pred1 , aes(y, mean)) +
  geom_errorbar(data = pred1 , aes(y, ymin= q0.025, ymax= q0.975), width=0.1, alpha=0.5)+
  geom_abline(slope=1, intercept=0, color='red') 

ggplot()+ geom_point(data = pred2 , aes(y, mean)) +
  geom_errorbar(data = pred2 , aes(y, ymin= q0.025, ymax= q0.975), width=0.1, alpha=0.5)+
  geom_abline(slope=1, intercept=0, color='red')

#Compute predictions from the model over the area of interest. Plot the posterior mean and the posterior sd.
pxl<-fm_pixels(mesh, mask = border)
pred_pxl1=predict(res, pxl,formula= ~Intercept1+ omega1)
pred_pxl2=predict(res, pxl,formula= ~Intercept2+ omega1_copy + omega2)
#plot mean and sd for the two fields
ggplot() + 
  gg(pred_pxl1, aes(color=mean)) + 
  ggtitle("Posterior mean for eta_1") +  xlab("") + ylab("")
ggplot() +
  gg(pred_pxl1, aes(color=sd)) + 
  ggtitle("Posterior sd for eta_1") +  xlab("") + ylab("")

ggplot() +
  gg(pred_pxl2, aes(color=mean)) + 
  ggtitle("Posterior mean for eta_2") +  xlab("") + ylab("")
ggplot() +
  gg(pred_pxl2, aes(color=sd)) + 
  ggtitle("Posterior sd for eta_2") +  xlab("") + ylab("")


#look at the posterior samples of the two fields
samples = generate(res, pxl,
                   ~ data.frame(omega1 = omega1,
                                omega2 = omega2),
                   n.samples = 5)


omega1 = sapply(samples, function(x) x$omega1)
p1 = cbind(pxl,omega1) %>%
  pivot_longer(-geometry) %>% ggplot() +
  geom_sf(aes(color =value)) + facet_wrap(.~name) + ggtitle("Omega 1")

omega2 = sapply(samples, function(x) x$omega2)
p2 = cbind(pxl,omega2) %>%
  pivot_longer(-geometry) %>% ggplot() +
  geom_sf(aes(color =value)) + facet_wrap(.~name) + ggtitle("Omega 2")


