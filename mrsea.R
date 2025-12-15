
library(sf)
library(ggplot2)
library(ggpubr)
library(INLA)
library(INLAspacetime)
library(inlabru)

## data
mrsea <- inlabru::mrsea

gg0 <- ggplot() + theme_minimal() + gg(mrsea$boundary)

## visualize the points over the samples and boundary 
gg0 +
  gg(mrsea$samplers) +
  gg(mrsea$points, size = 0.5) +
  geom_fm(data = mrsea$mesh, fill = 'transparent') +
  facet_wrap(~season) +
  ggtitle("MRSea observation seasons")

## deeper -> unlike to have a point
ggarrange(
    gg0 + geom_sf(aes(color = depth), mrsea$covar),
    gg0 + geom_sf(aes(color = depth), mrsea$samplers),
    gg0 + geom_sf(aes(color = depth), mrsea$points),
    nrow = 3
)

## defining the data model 
dataModel <-  bru_obs(
    geometry + season ~ .,
    family = "cp",
    data = mrsea$points,
    samplers = mrsea$samplers,
    domain = list(
        geometry = fmesher::fm_subdivide(mrsea$mesh, n = 3),
        season = seq_len(4)
    )
)

## depth covariate evaluator function, since eval_spatial() does not eval over POINTS
evDepth <- function(where) {
    i <- st_nearest_feature(where, mrsea$covar)
    mrsea$covar$depth[i]
}

## check
cor(evDepth(mrsea$covar), mrsea$covar$depth)
cor(evDepth(mrsea$samplers), mrsea$samplers$depth)
cor(evDepth(mrsea$points), mrsea$points$depth)

## define the components of model M0
cmp0 <- ~ 0 +
    season(season, model = 'factor_full') +
    Depth(evDepth(.data.), model = 'linear')

## fit M0
fit0 <- bru(cmp0, dataModel)

fit0$summary.fixed
fit0$summary.random$season

## model M1: smoothed depth effect

## setup of knots over depth
depth.bk <- pretty(mrsea$covar$depth, n = 31)
depth.bk
fDepth <- function(where, breaks = depth.bk) {
    n <- length(breaks)-1
    mids <- (breaks[1:n] + breaks[1 + 1:n])/2
    depth <- evDepth(where)
    mids[findInterval(depth, breaks)]    
}

## check
table(fDepth(mrsea$covar))
table(fDepth(mrsea$samplers))
table(fDepth(mrsea$points))

## define M1 components
cmp1 <- ~ 0 +
    season(season, model = 'factor_full') +
    Depth(fDepth(.data.), model = 'rw2', scale.model = TRUE,
          hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.1))))


fit1 <- bru(cmp1, dataModel)

fit1$summary.fixed
fit1$summary.random$season

## plot function for a smoothed effect
plotr <- function(r, add = FALSE, ...) {
    n <- nrow(r)
    if(!add) {
        plot(r$ID, r$mean,
             ylim = range(r[, c(4,6)]),
             type = 'n', ...)
    }
    polygon(r$ID[c(1:n, n:1, 1)],
            c(r[,4], rev(r[,6]), r[1,4]), 
            ...)
    lines(r$ID, r$mean, lty = 2)
}

par(bty = 'n')
plotr(fit1$summary.random$Depth,
      xlab = 'Depth', ylab = 'effect',
      col = gray(0.5, 0.5))
plot(function(x) mean(fit0$summary.random$season$mean) +
                 x * fit0$summary.fixed$mean, -30, 0,
     lwd = 2, add = TRUE)
legend('topleft', c("M0", "M1"), bty = 'n',
       lty = 1:2, lwd = 2:1)

## SPDE model
matern <- inla.spde2.pcmatern(
    mesh = mrsea$mesh,
    prior.sigma = c(0.1, 0.01),
    prior.range = c(10, 0.01),
    constr = TRUE
)

cmps <- update(
    cmp1, .~.+spatial(geometry, model = matern, group = season))

fits <- bru(cmps, dataModel)

grep('constraints', fits$logfile, value = TRUE) ## should be 1 + 1*4

fits$summary.random$season
fits$summary.hyperpar

par(bty = 'n')
plotr(fit1$summary.random$Depth,
      xlab = 'Depth', ylab = 'effect',
      col = gray(0.2, 0.5))
plot(function(x) mean(fit0$summary.random$season$mean) +
                 x * fit0$summary.fixed$mean, -30, 0,
     lwd = 2, add = TRUE)
plotr(fits$summary.random$Depth, add = TRUE,
      xlab = 'Depth', ylab = 'effect',
      col = rgb(0.5, 0.7, 1, 0.5))
legend('topleft', c("M0", "M1", "M2"), bty = 'n',
       lty = c(1,0,0), lwd = c(2,0,0), border = c('transparent', 1, 1), 
       fill = c('transparent', gray(0.2,0.5), rgb(.5,.7,1,0.5)))

## re-define a spacetime model
M3 <- update(
    cmp1, . ~ . +
              spacetime(list(space = geometry,
                             time = season),
                        model = stmodel)
)

## define a non-separable spacetime
stmodel <- stModel.define(
    smesh = mrsea$mesh,
    tmesh = fm_mesh_1d(loc = 1:4),
    model = '121',
    control.priors = list(
        prs = c(0.5, 0.1),
        prt = c(1.0, 0.1),
        psigma = c(1, 0.1)
    ),
    constr = TRUE
)

fitst <- bru(M3, dataModel)

fitst$cpu.used
fitst$summary.random$season
fitst$summary.hyperpar

stpmarginals <- c(
    list(inla.tmarginal(
        function(x) exp(-x/2),
        fitst$internal.marginals.hyperpar[[1]])),
    lapply(fitst$internal.marginals.hyperpar[2:4], function(m)
        inla.tmarginal(exp, m)))
names(stpmarginals) <- c('Depth sigma', 'Practical spatial range',
                         'Practical temporal range', 'spacetime sigma')

par(mfrow = c(2,2), mar = c(4,4,1,1))
for(k in 1:4)
    plot(stpmarginals[[k]], type = 'l',
         xlab = names(stpmarginals)[k], ylab = 'Density')


par(mfrow = c(1,1), mar = c(4,4,1,1), bty = 'n')
plotr(fit1$summary.random$Depth,
      xlab = 'Depth', ylab = 'effect',
      col = gray(0.2, 0.5))
plot(function(x) mean(fit0$summary.random$season$mean) +
                 x * fit0$summary.fixed$mean, -30, 0,
     lwd = 2, add = TRUE)
plotr(fits$summary.random$Depth, add = TRUE,
      xlab = 'Depth', ylab = 'effect',
      col = rgb(0.5, 0.7, 1, 0.5))
plotr(fitst$summary.random$Depth, add = TRUE,
      xlab = 'Depth', ylab = 'effect',
      col = rgb(1.0, 0.7, 0.5, 0.5))
legend('topleft', c("M0", "M1", "M2", "M3"), bty = 'n',
       lty = c(1,0,0,0), lwd = c(2,0,0,0), border = c('transparent', 1, 1,1), 
       fill = c('transparent', gray(0.2,0.5), rgb(.5,.7,1,0.5), rgb(1.0,.7,0.5,0.5)))


## Visualize the spatio temporal effect

nt <- c(mrsea$mesh, 4)
st.mean <- matrix(fitst$summary.random$spacetim$mean, ncol = 4)
dim(st.mean)

bb <- st_bbox(mrsea$boundary)
grid <- fm_evaluator(mesh = mrsea$mesh, xlim = bb[c(1,3)], ylim = bb[c(2,4)])

grid.sf <- st_as_sf(
    data.frame(grid$lattice$loc),
    coords = 1:2,
    crs = st_crs(mrsea$boundary)
)

st.grid <- lapply(1:4, function(j)
                  fm_evaluate(grid, field = st.mean[, j]))

grid.out <- which(sapply(
    st_within(grid.sf, mrsea$boundary), length)==0)

str(st.grid)

library(fields)
par(mfrow = c(2,2))
for(k in 1:4) {
    z <- st.grid[[k]]
    z[grid.out] <- NA
    image.plot(z, asp = 1)
}
