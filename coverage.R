prior.range = c(1, 0.5); prior.sigma = c(1, 0.1)
zlim.sd <- c(0,1.5)
zlim.mean <- c(-2, 2)
xlim.marginal = c(0,6)
#geom 
SM <- c(1.5)
W <- c(3)
R <- 2
N <- c(30)
fr <- c(0.01,0.5,0.8,1)

set.seed(880)
SS <- floor(runif(100, 1, 10000))

tbm.frac <- list()


  for (d in 1:length(SM)) {
    for (w in 1:length(W)) {
      for(r in 1:length(R)) {
        for (n in 1:length(N)) {
          
          smalldist = SM[d]
          width = c(W[w],W[w])
          range = R[r]
          max.edge.length = 0.4
          n = N[n]
          
          x.mid = 10
          xlim.big  = c(0,20)
          xlim.small = c(3,17)
          ylim.big = c(0,20)
          ylim.small = c(3,17)
          
          prior.range = prior.range
          prior.sigma = c(1,0.1)
          sigma.u = 1
          sigma.epsilon = 0.2
          prior.range.st = prior.range
          prior.sigma.st = c(1,0.1)
          
          
          
          poly1 <- local.square.polygon(xlim=c(xlim.big[1]-1, x.mid-smalldist/2), 
                                        ylim=x.mid+width[1]*c(-.5, .5))
          
          poly2 <- local.square.polygon(xlim=c(x.mid+smalldist/2, xlim.big[2]+1), 
                                        ylim=x.mid+width[2]*c(-.5, .5))
          
          poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))
          
          loc1 <- matrix(c(xlim.big[1],ylim.big[1], 
                           xlim.big[2],ylim.big[1], 
                           xlim.big[2],ylim.big[2], 
                           xlim.big[1],ylim.big[2]), 4, 2, byrow = T)
          
          seg <- inla.sp2segment(poly.original)
          # - Transforms a SpatialPolygon to an "inla polygon"
          mesh <- inla.mesh.2d(loc=loc1, 
                               interior = seg, 
                               max.e = max.edge.length, 
                               offset=1)
          
          tl <- length(mesh$graph$tv[,1])
          # - the number of triangles in the mesh
          posTri <- matrix(0, tl, 2)
          
          for (t in 1:tl){
            temp = mesh$loc[mesh$graph$tv[t, ], ]
            posTri[t,] = colMeans(temp)[c(1,2)] 
          }
          
          posTri <- SpatialPoints(posTri)
          # - the positions of the triangle centers
          
          bar.original <- over(poly.original, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar.original <- unlist(bar.original)
          poly.bar.orginal <- inla.barrier.polygon(mesh, barrier.triangles = bar.original)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 1
          bar1 <- over(poly1, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar1 <- unlist(bar1)
          poly.bar1 <- inla.barrier.polygon(mesh, barrier.triangles = bar1)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 2
          bar2 <- over(poly2, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar2 <- unlist(bar2)
          poly.bar2 <- inla.barrier.polygon(mesh, barrier.triangles = bar2)
          
          mat <-  inla.barrier.fem.plus(mesh, list(bar1, bar2)) 
          fem <- mat 
          barrier.triangles <- list(bar1, bar2)
          
          # PLOTS
          poly1_h <- local.square.polygon_T(xlim=c(xlim.small[1], x.mid-smalldist/2), 
                                            ylim=x.mid+width[1]*c(-.5, .5))
          
          poly2_h <- local.square.polygon_T(xlim=c(x.mid+smalldist/2, xlim.small[2]), 
                                            ylim=x.mid+width[2]*c(-.5, .5))
          
          loc2 <- matrix(c(3,3, 17,3, 17,17, 3,17), 4, 2, byrow = T)
          
          locp2 <- Polygon(loc2, hole = FALSE)
          
          poly.water <- SpatialPolygons(list(Polygons(list(locp2, poly1_h, poly2_h), '0')))
          poly.water_sf <- st_as_sf(poly.water)
          
          for (s in 34:length(SS)) {
            set.inla.seed = SS[s]
            set.seed(set.inla.seed)
            loc.data <- spsample(x = poly.water, n = n, type = "random")
            loc.data_sf <- st_as_sf(loc.data)
            loc.data <- loc.data@coords
          
            poly.rect <- SpatialPolygons(list(Polygons(list(poly1_h, poly2_h), '0')))
            poly.rect_sf <- st_as_sf(poly.rect)
            poly.sq <- SpatialPolygons(list(Polygons(list(locp2), '0')))
            poly.sq_sf <- st_as_sf(poly.sq)
          
            sq_bbox <- st_bbox(poly.sq_sf) %>% st_as_sfc()
            poly1_bbox <- st_bbox(poly1) %>% st_as_sfc()
            poly2_bbox <- st_bbox(poly2) %>% st_as_sfc()
          
          
            tbm.frac[[s]] <- list()
          
            for (p in 1:length(fr)) { #
            #fp <- file.path("plots/geom1.2_2/coverage", paste0("sm", SM[d], "w", W[w], "r", R[r], "seed", SS[s], ".n30", "fr",fr[p], "/"))
            #if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
            #gif.path = fp
            #marginals.path = fp
            
                f <- fr[p]
                range.fraction <- c(0.01, f)
             
                tbm.frac[[s]][[p]] <- model.tbm.plus(mesh = mesh, fem = fem, 
                                            barrier.triangles = barrier.triangles,
                                            prior.range = prior.range, prior.sigma = prior.sigma,  
                                            range.fraction = range.fraction,
                                            range = range,
                                            set.inla.seed = set.inla.seed,
                                            loc.data = loc.data,
                                            sigma.u = 1, sigma.epsilon = 0.2,
                                            poly.original = poly.bar.orginal,
                                            prior.range.st = prior.range.st,
                                            prior.sigma.st = prior.sigma.st,
                                            return.list = TRUE)
            
            
            
          }
          
          
        }
        
      }
    }
  }
  
}

#tbm.frac.config1_34.67 <- tbm.frac
#save(tbm.frac.config1_34.67, file = "plots/coverage/config1_34.100.RData")
#save objects as .qs instead of RData that is way too heavy

#####################################################################################################################################
#sm0w1r4seed987.n50
zlim.sd <- c(0,1.4)
zlim.mean <- c(-2, 2.822)


prior.range = c(1, 0.5); prior.sigma = c(1, 0.1)
#geom 
SM <- c(0)
W <- c(1)
R <- 4

N <- c(50)
fr <- c(0.01,0.5,0.8,1)

set.seed(880)
SS <- floor(runif(100, 1, 10000))

SS <- SS[seq(31, 100, 1)]
tbm.frac2 <- list()


  for (d in 1:length(SM)) {
    for (w in 1:length(W)) {
      for(r in 1:length(R)) {
        for (n in 1:length(N)) {
          
          
          
          
          smalldist = SM[d]
          width = c(W[w],W[w])
          range = R[r]
          max.edge.length = 0.4
          n = N[n]
          
          x.mid = 10
          xlim.big  = c(0,20)
          xlim.small = c(3,17)
          ylim.big = c(0,20)
          ylim.small = c(3,17)
          
          prior.range = prior.range
          prior.sigma = c(1,0.1)
          sigma.u = 1
          sigma.epsilon = 0.2
          prior.range.st = prior.range
          prior.sigma.st = c(1,0.1)
          
          
          
          poly1 <- local.square.polygon(xlim=c(xlim.big[1]-1, x.mid-smalldist/2), 
                                        ylim=x.mid+width[1]*c(-.5, .5))
          
          poly2 <- local.square.polygon(xlim=c(x.mid+smalldist/2, xlim.big[2]+1), 
                                        ylim=x.mid+width[2]*c(-.5, .5))
          
          poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))
          
          loc1 <- matrix(c(xlim.big[1],ylim.big[1], 
                           xlim.big[2],ylim.big[1], 
                           xlim.big[2],ylim.big[2], 
                           xlim.big[1],ylim.big[2]), 4, 2, byrow = T)
          
          seg <- inla.sp2segment(poly.original)
          # - Transforms a SpatialPolygon to an "inla polygon"
          mesh <- inla.mesh.2d(loc=loc1, 
                               interior = seg, 
                               max.e = max.edge.length, 
                               offset=1)
          
          tl <- length(mesh$graph$tv[,1])
          # - the number of triangles in the mesh
          posTri <- matrix(0, tl, 2)
          
          for (t in 1:tl){
            temp = mesh$loc[mesh$graph$tv[t, ], ]
            posTri[t,] = colMeans(temp)[c(1,2)] 
          }
          
          posTri <- SpatialPoints(posTri)
          # - the positions of the triangle centers
          
          bar.original <- over(poly.original, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar.original <- unlist(bar.original)
          poly.bar.orginal <- inla.barrier.polygon(mesh, barrier.triangles = bar.original)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 1
          bar1 <- over(poly1, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar1 <- unlist(bar1)
          poly.bar1 <- inla.barrier.polygon(mesh, barrier.triangles = bar1)
          # - the Barrier model's polygon
          # - in most cases this should be the same as poly.original
          
          # BARRIER 2
          bar2 <- over(poly2, SpatialPoints(posTri), returnList=T)
          # - checking which mesh triangles are inside the barrier area
          bar2 <- unlist(bar2)
          poly.bar2 <- inla.barrier.polygon(mesh, barrier.triangles = bar2)
          
          mat <-  inla.barrier.fem.plus(mesh, list(bar1, bar2)) 
          fem <- mat 
          barrier.triangles <- list(bar1, bar2)
          
          # PLOTS
          poly1_h <- local.square.polygon_T(xlim=c(xlim.small[1], x.mid-smalldist/2), 
                                            ylim=x.mid+width[1]*c(-.5, .5))
          
          poly2_h <- local.square.polygon_T(xlim=c(x.mid+smalldist/2, xlim.small[2]), 
                                            ylim=x.mid+width[2]*c(-.5, .5))
          
          loc2 <- matrix(c(3,3, 17,3, 17,17, 3,17), 4, 2, byrow = T)
          
          locp2 <- Polygon(loc2, hole = FALSE)
          
          poly.water <- SpatialPolygons(list(Polygons(list(locp2, poly1_h, poly2_h), '0')))
          poly.water_sf <- st_as_sf(poly.water)
          
          for (s in 1:length(SS)) { #length(SS)) {
            set.inla.seed = SS[s]
            set.seed(set.inla.seed)
            loc.data <- spsample(x = poly.water, n = n, type = "random")
            loc.data_sf <- st_as_sf(loc.data)
            loc.data <- loc.data@coords
          
            poly.rect <- SpatialPolygons(list(Polygons(list(poly1_h, poly2_h), '0')))
            poly.rect_sf <- st_as_sf(poly.rect)
            poly.sq <- SpatialPolygons(list(Polygons(list(locp2), '0')))
            poly.sq_sf <- st_as_sf(poly.sq)
          
            sq_bbox <- st_bbox(poly.sq_sf) %>% st_as_sfc()
            poly1_bbox <- st_bbox(poly1) %>% st_as_sfc()
            poly2_bbox <- st_bbox(poly2) %>% st_as_sfc()
        
          
            tbm.frac2[[s]] <- list()
            for (p in 1:length(fr)) { #
            #fp <- file.path("plots/geom1.2_2/coverage", paste0("sm", SM[d], "w", W[w], "r", R[r], "seed", SS[s], ".n30", "fr",fr[p], "/"))
            #if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
            #gif.path = fp
            #marginals.path = fp
            
                f <- fr[p]
                range.fraction <- c(0.01, f)
            
                tbm.frac2[[s]][[p]] <- model.tbm.plus(mesh = mesh, fem = fem, 
                                                  barrier.triangles = barrier.triangles,
                                                  prior.range = prior.range, prior.sigma = prior.sigma,  
                                                  range.fraction = range.fraction,
                                                  range = range,
                                                  set.inla.seed = set.inla.seed,
                                                  loc.data = loc.data,
                                                  sigma.u = 1, sigma.epsilon = 0.2,
                                                  poly.original = poly.bar.orginal,
                                                  prior.range.st = prior.range.st,
                                                  prior.sigma.st = prior.sigma.st,
                                                  return.list = TRUE)
            
          }
          
          
          }
          
          
        
      }
    }
  }
  
}

#################################################################################################################################
#CONFIG2 FROM 31 TO 100

#tbm.frac.config2 <- tbm.frac2
#tbm.frac.config2_21.30 <- tbm.frac2
#save(tbm.frac.config2_21.30, file = "plots/coverage/config2_21.30.RData")

# Install if needed
#install.packages("qs")
#to actually use qs I need
#install.packages(c("qs", "stringfish", "RcppParallel"), type = "source")
#so I won't do it now and do it later so I can Restart R (in RStudio or terminal) to ensure a clean session.

# Save object
#qs::qsave(tbm.frac.config2_21.30, "plots/coverage/tbm.frac.config2_30.qs")

# Load it later
#my_loaded_object <- qs::qread("my_large_object.qs")

#use rds instead now
saveRDS(tbm.frac.config2_21.30, "plots/coverage/config2_21.30.rds")

#remove heavy object
rm(tbm.frac.config2_21.30)
rm(tbm.frac.config2)
rm(tbm.frac2)
#check
ls()
#free up memory, r doesn't release memory immediately after rm(), so i can trigger garbage collection manually
gc()


tbm.frac.config2_31.45 <- tbm.frac2[(seq(1,15,1))]
tbm.frac.config2_46.60 <- tbm.frac2[(seq(16,30,1))]
tbm.frac.config2_61.75 <- tbm.frac2[(seq(31,45,1))]
tbm.frac.config2_76.90 <- tbm.frac2[(seq(46,60,1))]
tbm.frac.config2_91.100 <- tbm.frac2[(seq(61,70,1))]

saveRDS(tbm.frac.config2_31.45, "plots/coverage/config2_31.45.rds")
saveRDS(tbm.frac.config2_46.60, "plots/coverage/config2_46.60.rds")
saveRDS(tbm.frac.config2_61.75, "plots/coverage/config2_61.75.rds")
saveRDS(tbm.frac.config2_76.90, "plots/coverage/config2_79.90.rds")
saveRDS(tbm.frac.config2_91.100, "plots/coverage/config2_91.100.rds")

#################################################################################################################################
#coverage for config 2 for to 31:100 seed

#set.seed(880); SS <- floor(runif(100, 1, 10000)); SS <- SS[seq(31, 100, 1)]

rule.tbm <- list()
rule.bm <- list()
rule.st <- list()
ci.tbm <- list() #length of confidence interval
ci.bm <- list()
ci.st <- list()

for(s in 1:length(SS)){ #SS from 31 to 100
  rule.tbm[[s]] <- list() # rule to chek if mean is inside ci
  rule.bm[[s]] <- list()
  rule.st[[s]] <- list()
  
  ci.tbm[[s]] <- list() #length of confidence interval
  ci.bm[[s]] <- list()
  ci.st[[s]] <- list()
  
  for (p in 1:length(fr)) {
    #rule
    if(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,3] < tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,1] &&
       tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,1] < tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,5]) {
      rule.tbm[[s]][[p]] <- 1
    } else {
      rule.tbm[[s]][[p]] <- 0
    }
    
    if(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,3] < tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,1] &&
       tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,1] < tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,5]) {
      rule.bm[[s]][[p]] <- 1
    } else {
      rule.bm[[s]][[p]] <- 0
    }
    
    if(tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,3] < tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,1] &&
       tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,1] < tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,5]) {
      rule.st[[s]][[p]] <- 1
    } else {
      rule.st[[s]][[p]] <- 0
    }
    
    #length 
    ci.tbm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,5]) - exp(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,3])
    ci.bm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,5]) - exp(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,3])
    ci.st[[s]][[p]] <- tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,5] - tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,3]
  }
}

#config2.rule_31.100 <- list(rule.tbm=rule.tbm, rule.bm=rule.bm, rule.st=rule.st)
#config2.ci.length_31.100 <- list(ci.tbm=ci.tbm, ci.bm=ci.bm, ci.st=ci.st)

u.range.tbm <- list()
u.range.bm <- list()
u.range.st <- list()

for(s in 1:length(SS)){ #SS from 31 to 100
  u.range.tbm[[s]] <- list()
  u.range.bm[[s]] <- list()
  u.range.st[[s]] <- list()
  for (p in 1:length(fr)) {
    u.range.tbm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,1])
    u.range.bm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,1])
    u.range.st[[s]][[p]] <- tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,1]
  }
}
#config2.urange.length_31.100 <- list(u.range.tbm=u.range.tbm, u.range.bm=u.range.bm, u.range.st=u.range.st)

#results
sum(unlist(rule.tbm))
sum(unlist(rule.bm))
sum(unlist(rule.st))

sum(unlist(ci.tbm))
sum(unlist(ci.bm))
sum(unlist(ci.st))

cbind((unlist(ci.tbm)),(unlist(ci.bm)),(unlist(ci.st)))

cbind(unlist(u.range.tbm), unlist(u.range.bm), unlist(u.range.st))

#remove heavy object
rm(tbm.frac2)
rm(tbm.frac.config2_31.45)
rm(tbm.frac.config2_46.60)
rm(tbm.frac.config2_61.75)
rm(tbm.frac.config2_76.90)
rm(tbm.frac.config2_91.100)

ls()
gc()

#conclusion the difference is on local mean because overall they all work and the range between tbm and bm are almost the same, the stationary mean is not far but it is larger

#################################################################################################################################
#CONFIG2 FROM 1 TO 20
load("plots/coverage/config2_1.20.RData")
tbm.frac2 <- tbm.frac.config2

rule.tbm <- list()
rule.bm <- list()
rule.st <- list()
ci.tbm <- list() #length of confidence interval
ci.bm <- list()
ci.st <- list()

for(s in 1:length(tbm.frac2)){ #SS from 31 to 100
  rule.tbm[[s]] <- list() # rule to chek if mean is inside ci
  rule.bm[[s]] <- list()
  rule.st[[s]] <- list()
  
  ci.tbm[[s]] <- list() #length of confidence interval
  ci.bm[[s]] <- list()
  ci.st[[s]] <- list()
  
  for (p in 1:length(tbm.frac2[[1]])) {
    #rule
    if(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,3] < tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,1] &&
       tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,1] < tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,5]) {
      rule.tbm[[s]][[p]] <- 1
    } else {
      rule.tbm[[s]][[p]] <- 0
    }
    
    if(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,3] < tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,1] &&
       tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,1] < tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,5]) {
      rule.bm[[s]][[p]] <- 1
    } else {
      rule.bm[[s]][[p]] <- 0
    }
    
    if(tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,3] < tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,1] &&
       tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,1] < tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,5]) {
      rule.st[[s]][[p]] <- 1
    } else {
      rule.st[[s]][[p]] <- 0
    }
    
    #length 
    ci.tbm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,5]) - exp(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,3])
    ci.bm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,5]) - exp(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,3])
    ci.st[[s]][[p]] <- tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,5] - tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,3]
  }
}

config2.rule_1.20 <- list(rule.tbm=rule.tbm, rule.bm=rule.bm, rule.st=rule.st)
config2.ci.length_1.20 <- list(ci.tbm=ci.tbm, ci.bm=ci.bm, ci.st=ci.st)

u.range.tbm <- list()
u.range.bm <- list()
u.range.st <- list()

for(s in 1:length(tbm.frac2)){ #SS from 31 to 100
  u.range.tbm[[s]] <- list()
  u.range.bm[[s]] <- list()
  u.range.st[[s]] <- list()
  for (p in 1:length(tbm.frac2[[1]])) {
    u.range.tbm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,1])
    u.range.bm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,1])
    u.range.st[[s]][[p]] <- tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,1]
  }
}
config2.urange.length_1.20 <- list(u.range.tbm=u.range.tbm, u.range.bm=u.range.bm, u.range.st=u.range.st)
summary.config2_1.20 <- list(rule = config2.rule_1.20, ci.length = config2.ci.length_1.20, u.range = config2.urange.length_1.20)
#saveRDS(summary.config2_1.20, "plots/coverage/summary.config2_1.20.rds")

#RESULTS
sum(unlist(rule.tbm))
sum(unlist(rule.bm))
sum(unlist(rule.st))

sum(unlist(ci.tbm))
sum(unlist(ci.bm))
sum(unlist(ci.st))

cbind((unlist(ci.tbm)),(unlist(ci.bm)),(unlist(ci.st)))

cbind(unlist(u.range.tbm), unlist(u.range.bm), unlist(u.range.st))


#################################################################################################################################
tbm.frac2 <- readRDS("plots/coverage/config2_21.30.rds")

rule.tbm <- list()
rule.bm <- list()
rule.st <- list()
ci.tbm <- list() #length of confidence interval
ci.bm <- list()
ci.st <- list()

for(s in 1:length(tbm.frac2)){ #SS from 31 to 100
  rule.tbm[[s]] <- list() # rule to chek if mean is inside ci
  rule.bm[[s]] <- list()
  rule.st[[s]] <- list()
  
  ci.tbm[[s]] <- list() #length of confidence interval
  ci.bm[[s]] <- list()
  ci.st[[s]] <- list()
  
  for (p in 1:length(tbm.frac2[[1]])) {
    #rule
    if(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,3] < tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,1] &&
       tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,1] < tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,5]) {
      rule.tbm[[s]][[p]] <- 1
    } else {
      rule.tbm[[s]][[p]] <- 0
    }
    
    if(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,3] < tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,1] &&
       tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,1] < tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,5]) {
      rule.bm[[s]][[p]] <- 1
    } else {
      rule.bm[[s]][[p]] <- 0
    }
    
    if(tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,3] < tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,1] &&
       tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,1] < tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,5]) {
      rule.st[[s]][[p]] <- 1
    } else {
      rule.st[[s]][[p]] <- 0
    }
    
    #length 
    ci.tbm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,5]) - exp(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,3])
    ci.bm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,5]) - exp(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,3])
    ci.st[[s]][[p]] <- tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,5] - tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,3]
  }
}

config2.rule_1.30 <- list(rule.tbm=rule.tbm, rule.bm=rule.bm, rule.st=rule.st)
config2.ci.length_1.30 <- list(ci.tbm=ci.tbm, ci.bm=ci.bm, ci.st=ci.st)

u.range.tbm <- list()
u.range.bm <- list()
u.range.st <- list()

for(s in 1:length(tbm.frac2)){ #SS from 31 to 100
  u.range.tbm[[s]] <- list()
  u.range.bm[[s]] <- list()
  u.range.st[[s]] <- list()
  for (p in 1:length(tbm.frac2[[1]])) {
    u.range.tbm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.tbm$res$summary.hyperpar[3,1])
    u.range.bm[[s]][[p]] <- exp(tbm.frac2[[s]][[p]]$pos.bm$res$summary.hyperpar[3,1])
    u.range.st[[s]][[p]] <- tbm.frac2[[s]][[p]]$pos.st$res$summary.hyperpar[2,1]
  }
}
config2.urange.length_1.30 <- list(u.range.tbm=u.range.tbm, u.range.bm=u.range.bm, u.range.st=u.range.st)
summary.config2_1.30 <- list(rule = config2.rule_1.30, ci.length = config2.ci.length_1.30, u.range = config2.urange.length_1.30)
saveRDS(summary.config2_1.30, "plots/coverage/summary.config2_1.30.rds")

#################################################################################################################################

#RESULTS for ALL CONFIG2

sum.config2.30 <- readRDS("plots/coverage/summary.config2_1.30.rds")
sum.config2.31.100 <- readRDS("plots/coverage/summary.config2_31.100.rds")

ones.tbm <- list()
ones.bm <- list()
ones.st <- list()
for (s in 1:30) {
  ones.tbm[[s]] <- sum.config2.30$rule$rule.tbm[[s]][[1]]
  ones.bm[[s]] <- sum.config2.30$rule$rule.bm[[s]][[1]]
  ones.st[[s]] <- sum.config2.30$rule$rule.st[[s]][[1]]
}

sum(unlist(ones.tbm)) 
sum(unlist(ones.bm))
sum(unlist(ones.st))

#RESULTS:
#ALL SUM UP TO 30 from 1 to 30

ones.tbm <- list()
ones.bm <- list()
ones.st <- list()
for (s in 1:70) {
  ones.tbm[[s]] <- sum.config2.31.100$rule$rule.tbm[[s]][[4]]
  ones.bm[[s]] <- sum.config2.31.100$rule$rule.bm[[s]][[4]]
  ones.st[[s]] <- sum.config2.31.100$rule$rule.st[[s]][[4]]
}

sum(unlist(ones.tbm)) 
sum(unlist(ones.bm))
sum(unlist(ones.st))

#RESULTS:
#ALL SUM UP TO 30 from 31 to 70
# 0.01
length.tbm <- list()
length.bm <- list()
length.st <- list()
for (s in 1:70) {
  length.tbm[[s]] <- sum.config2.31.100$ci.length$ci.tbm[[s]][[1]]
  length.bm[[s]] <- sum.config2.31.100$ci.length$ci.bm[[s]][[1]]
  length.st[[s]] <- sum.config2.31.100$ci.length$ci.st[[s]][[1]]
}

length.tbm_30 <- list()
length.bm_30 <- list()
length.st_30 <- list()
for (s in 1:30) {
  length.tbm_30[[s]] <- sum.config2.30$ci.length$ci.tbm[[s]][[1]]
  length.bm_30[[s]] <- sum.config2.30$ci.length$ci.bm[[s]][[1]]
  length.st_30[[s]] <- sum.config2.30$ci.length$ci.st[[s]][[1]]
}

u.length.tbm1 <- mean(c(unlist(length.tbm), unlist(length.tbm_30)))
u.length.bm1 <- mean(c(unlist(length.bm), unlist(length.bm_30)))
u.length.st1 <- mean(c(unlist(length.st), unlist(length.st_30)))

# 0.5
length.tbm <- list()
length.bm <- list()
length.st <- list()
for (s in 1:70) {
  length.tbm[[s]] <- sum.config2.31.100$ci.length$ci.tbm[[s]][[2]]
  length.bm[[s]] <- sum.config2.31.100$ci.length$ci.bm[[s]][[2]]
  length.st[[s]] <- sum.config2.31.100$ci.length$ci.st[[s]][[2]]
}

length.tbm_30 <- list()
length.bm_30 <- list()
length.st_30 <- list()
for (s in 1:30) {
  length.tbm_30[[s]] <- sum.config2.30$ci.length$ci.tbm[[s]][[2]]
  length.bm_30[[s]] <- sum.config2.30$ci.length$ci.bm[[s]][[2]]
  length.st_30[[s]] <- sum.config2.30$ci.length$ci.st[[s]][[2]]
}

u.length.tbm2 <- mean(c(unlist(length.tbm), unlist(length.tbm_30)))
u.length.bm2 <- mean(c(unlist(length.bm), unlist(length.bm_30)))
u.length.st2 <- mean(c(unlist(length.st), unlist(length.st_30)))

#0.8
length.tbm <- list()
length.bm <- list()
length.st <- list()
for (s in 1:70) {
  length.tbm[[s]] <- sum.config2.31.100$ci.length$ci.tbm[[s]][[3]]
  length.bm[[s]] <- sum.config2.31.100$ci.length$ci.bm[[s]][[3]]
  length.st[[s]] <- sum.config2.31.100$ci.length$ci.st[[s]][[3]]
}

length.tbm_30 <- list()
length.bm_30 <- list()
length.st_30 <- list()
for (s in 1:30) {
  length.tbm_30[[s]] <- sum.config2.30$ci.length$ci.tbm[[s]][[3]]
  length.bm_30[[s]] <- sum.config2.30$ci.length$ci.bm[[s]][[3]]
  length.st_30[[s]] <- sum.config2.30$ci.length$ci.st[[s]][[3]]
}

u.length.tbm3 <- mean(c(unlist(length.tbm), unlist(length.tbm_30)))
u.length.bm3 <- mean(c(unlist(length.bm), unlist(length.bm_30)))
u.length.st3 <- mean(c(unlist(length.st), unlist(length.st_30)))

#1
length.tbm <- list()
length.bm <- list()
length.st <- list()
for (s in 1:70) {
  length.tbm[[s]] <- sum.config2.31.100$ci.length$ci.tbm[[s]][[4]]
  length.bm[[s]] <- sum.config2.31.100$ci.length$ci.bm[[s]][[4]]
  length.st[[s]] <- sum.config2.31.100$ci.length$ci.st[[s]][[4]]
}

length.tbm_30 <- list()
length.bm_30 <- list()
length.st_30 <- list()
for (s in 1:30) {
  length.tbm_30[[s]] <- sum.config2.30$ci.length$ci.tbm[[s]][[4]]
  length.bm_30[[s]] <- sum.config2.30$ci.length$ci.bm[[s]][[4]]
  length.st_30[[s]] <- sum.config2.30$ci.length$ci.st[[s]][[4]]
}

u.length.tbm4 <- mean(c(unlist(length.tbm), unlist(length.tbm_30)))
u.length.bm4 <- mean(c(unlist(length.bm), unlist(length.bm_30)))
u.length.st4 <- mean(c(unlist(length.st), unlist(length.st_30)))

data.frame(length.tbm = c(u.length.tbm1, u.length.tbm2, u.length.tbm3, u.length.tbm4),
           length.bm = c(u.length.bm1, u.length.bm2, u.length.bm3, u.length.bm4),
           length.st = c(u.length.st1, u.length.st2, u.length.st3, u.length.st4))
#################################################################################################################################


















#################################################################################################################################
