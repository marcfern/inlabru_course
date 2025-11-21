## ----child="practicals/distance_sampling.qmd"---------------------------------

## -----------------------------------------------------------------------------
#| warning: false
#| message: false


library(dplyr)
library(INLA)
library(ggplot2)
library(patchwork)
library(inlabru)     
library(sf)
# load some libraries to generate nice map plots
library(scico)


## -----------------------------------------------------------------------------

mexdolphin <- mexdolphin_sf

mexdolphin$depth <- mexdolphin$depth %>% mutate(depth=scale(depth)%>%c())

ggplot() + geom_sf(data = mexdolphin$points, color = "red" ) +

  geom_sf(data = mexdolphin$samplers) +

  geom_sf(data = mexdolphin$ppoly, alpha = 0)
          



## -----------------------------------------------------------------------------

boundary0 = fm_nonconvex_hull(mexdolphin$points,convex = -0.1)
mesh_0 = fm_mesh_2d(boundary = boundary0,
                          max.edge = c(30, 150), # The largest allowed triangle edge length.
                          cutoff = 15,
                          crs = fm_crs(mexdolphin$points))
ggplot() + gg(mesh_0)



## -----------------------------------------------------------------------------

mesh_1 = fm_mesh_2d(boundary = mexdolphin$ppoly,
                    max.edge = c(30, 150),
                    cutoff = 15,
                    crs = fm_crs(mexdolphin$points))

ggplot() + gg(mesh_1)





## -----------------------------------------------------------------------------

spde_model <- inla.spde2.pcmatern(mexdolphin$mesh,
  prior.sigma = c(2, 0.01),
  prior.range = c(50, 0.01)
)





## -----------------------------------------------------------------------------

hn <- function(distance, sigma) {
  exp(-0.5 * (distance / sigma)^2)

}



## -----------------------------------------------------------------------------

cmp <- ~ space(main = geometry, model = spde_model) +
  sigma(1,
    prec.linear = 1,
    marginal = bm_marginal(qexp, pexp, dexp, rate = 1 / 8)
  ) +
  Intercept(1)



## -----------------------------------------------------------------------------
#| echo: true
#| eval: false

 eta <- geometry+distance~ space+log(hn(distance,sigma))+Intercept+ log(2)





## -----------------------------------------------------------------------------
# build integration scheme

distance_domain <-  fm_mesh_1d(seq(0, 8,
                              length.out = 30))

ips = fm_int(list(geometry = mexdolphin$mesh,
                  distance = distance_domain),
             samplers = mexdolphin$samplers)



## -----------------------------------------------------------------------------

lik = bru_obs("cp",
              formula = eta,
              data = mexdolphin$points,
              ips = ips)



## -----------------------------------------------------------------------------
fit = bru(cmp, lik)


fit$summary.fixed
fit$summary.hyperpar



## -----------------------------------------------------------------------------

plot( spde.posterior(fit, "space", what = "range")) +
plot( spde.posterior(fit, "space", what = "log.variance"))





## -----------------------------------------------------------------------------

ggplot() + geom_sf(data = pr.int$spatial,aes(color = mean)) + scale_color_scico() + ggtitle("Posterior mean")

ggplot() + geom_sf(data = pr.int$spatial,aes(color = sd)) + scale_color_scico() + ggtitle("Posterior sd")



pxl<-fm_pixels(mexdolphin$mesh, dims=c(200,100), mask=mexdolphin$ppoly)
pr.pxl=predict(fit, pxl, ~data.frame(space=space,
                                     lambda=exp(Intercept+space)))

ggplot() + geom_sf(data = pr.pxl$space,aes(color = mean)) + scale_color_scico() + ggtitle("Posterior mean")
ggplot() + geom_sf(data = pr.pxl$space,aes(color = sd)) + scale_color_scico() + ggtitle("Posterior sd")


ggplot() +
  
  geom_sf(data = pr.pxl$lambda,aes(color = mean)) +
  
  scale_color_scico(palette = "imola") +
  
  ggtitle("Posterior mean")
## -----------------------------------------------------------------------------

distdf <- data.frame(distance = seq(0, 8, length.out = 100))

dfun <- predict(fit, distdf, ~ hn(distance, sigma))

plot(dfun)



## -----------------------------------------------------------------------------

predpts <- fm_int(mexdolphin$mesh, mexdolphin$ppoly)
Lambda <- predict(fit, predpts, ~ sum(weight * exp(space + Intercept)))
Lambda



## -----------------------------------------------------------------------------
Ns <- seq(50, 450, by = 1)

Nest <- predict(fit, predpts,
  ~ data.frame(
    N = Ns,
    density = dpois(
      Ns,
      lambda = sum(weight * exp(space + Intercept))
    )
  ),
  n.samples = 2000
)



## -----------------------------------------------------------------------------

Nest <- dplyr::bind_rows(
  cbind(Nest, Method = "Posterior"),
  data.frame(
    N = Nest$N,
    mean = dpois(Nest$N, lambda = Lambda$mean),
    mean.mc_std_err = 0,
    Method = "Plugin"
  )
)



## -----------------------------------------------------------------------------

ggplot(data = Nest) +
  geom_line(aes(x = N, y = mean, colour = Method)) +
  geom_ribbon(
    aes(
      x = N,
      ymin = mean - 2 * mean.mc_std_err,
      ymax = mean + 2 * mean.mc_std_err,
      fill = Method,
    ),
    alpha = 0.2
  ) +
  geom_line(aes(x = N, y = mean, colour = Method)) +
  ylab("Probability mass function")



