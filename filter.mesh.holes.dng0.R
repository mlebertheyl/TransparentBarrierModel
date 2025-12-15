# https://www.flaticon.com/
# dng_icon <- "dng.data/dugong.png"

bathym_ <- raster("dng.data/RSPbathym100m.tif")
crs_longlat <- CRS("+proj=longlat +datum=WGS84")
crs_km <- CRS("+proj=utm +zone=39 ellps=WGS84 +units=km")
bathym_longlat <- projectRaster(bathym_, crs = crs_longlat)
bathym_km <- projectRaster(bathym_, crs = crs_km)

# Data for all species from helicopter surveys.
Data_allspp <- readRDS("dng.data/RSP_survey_sightings.rds")
# Data Incidental Sightings Dugongs

# Dugong Example
## Data

Data_ISD <- readRDS("dng.data/incidental_sightings_dugongs.rds")

# Data for the 10 different species. 
# The list number corresponds to sppID, e.g. species 7 is spp[[7]]
spp_ <- 10
spp <- list()
for (s in 1:spp_) {
  data <- Data_allspp %>% filter(sp_grp_nam == paste0("spp", as.character(s)))
  spp[[s]] <- data
}

# Dugong = spp7
Data_7spp <- spp[[7]]

# bind ISD with Survey helicopter data
# Dugong is the only one that has data coming from 2 sources.
Data_ISD_ <- dplyr::select(Data_ISD, c(Lat, Long, geometry)) %>% 
  rename(lat = Lat, long = Long) %>% 
  mutate(method = "ISD")

spp7 <- dplyr::select(Data_7spp, c(lat, long, geometry)) %>% 
  mutate(method = "Survey")

all_dugong <- rbind(Data_ISD_, spp7)
st_crs(Data_ISD_)
Data_ISD <- st_transform(Data_ISD_, crs_km)
cords_Dng.xy <- st_coordinates(Data_ISD)
cords_Dng.sp <- SpatialPoints(cords_Dng.xy, proj4string = crs_km)

# line 15 "mybath"
# creating polygons with terra
bathym <- terra::rast(bathym_longlat)
plot(bathym)
bathym_poly <- terra::as.polygons(bathym > -Inf)
plot(bathym_poly)
# getting coordinates for all the polygons
cords <- terra::crds(bathym)
# getting bathymetry values for all polygons
vals <- terra::extract(bathym, as.data.frame(cords))
# create df with cords and bathym values.
df_bath <- cbind(as.data.frame(cords), vals) %>% rename(bathymetry = layer)
# filter polygons corresponding to islands and sand patches, select x, y ordinates.
sea_cords <- df_bath %>% dplyr::filter(bathymetry <= -2) %>% dplyr::select(x,y)
island_cords <- df_bath %>% dplyr::filter(bathymetry > -2) %>% dplyr::select(x,y)

# boundary of the study area:
# raster w/o the boundary
raster_study.area <- bathym_longlat
raster_study.area[raster_study.area < -500] <- NA
# opposite
raster_deep.sea<- bathym_longlat
raster_deep.sea[raster_deep.sea>= -500] <- NA
# plot
plot(raster_study.area); plot(raster_deep.sea)

# raster w/o islands
# islands have r = 0.01
raster_no.islands <- raster_study.area
raster_no.islands[raster_no.islands >= 0] <- NA # create raster with NA values for islands.
# opposite
# boundary and islands have range = 0.01
raster_islands <- bathym_longlat
raster_islands[raster_islands < 0 & raster_islands >= -500] <- NA 
# plot
plot(raster_no.islands); plot(raster_islands)
raster_b1 <- raster_islands

# raster w/o islands (b1) and b2 
# 90% between -10 and 0 depth, https://www.nature.com/articles/s41598-021-04412-3
# ratio = 0.1/0.9 ~ 0.11
raster_no.b2 <- raster_no.islands
raster_no.b2[raster_no.b2 < -50] <- NA
# opposite
raster_b2 <- raster_no.islands
raster_b2[raster_b2 >= -50] <- NA
# plot
plot(raster_no.b2); plot(raster_b2)
#-50 effectively filters the barriers and leaves the area in the middle as normal area

# i could, in theory, have a third barrier
raster_no.b3 <- raster_no.b2
raster_no.b3[raster_no.b3 < -30] <- NA
# opposite
raster_b3 <- raster_no.b2
raster_b3[raster_b3 >= -30] <- NA
# plot
plot(raster_no.b3); plot(raster_b3)


# from terra lib
# SpatRaster
spatr_no.island <- terra::rast(raster_no.islands)
# SpatVector
spatvect_no.island <- terra::as.polygons(spatr_no.island > -Inf) # bathym_poly
#SpatialPolygonsDataFrame 
spdf.water <- as(spatvect_no.island, "Spatial") # bathym_sp_df
# for the mesh I need b2 without NA

# There's one observation that is out of the bathymetry map
spdf.water@bbox
Data_ISD <-
  Data_ISD %>% filter(long > spdf.water@bbox[1,1] &long < spdf.water@bbox[1,2] 
                      &lat > spdf.water@bbox[2,1] &lat < spdf.water@bbox[2,2])
# data inside the mesh box                  
cords_Dng.xy <- st_coordinates(Data_ISD)
cords_Dng.sp <- SpatialPoints(cords_Dng.xy, proj4string = crs_km)

sp4msh <- spTransform(spdf.water, crs_km)
spdf.water <- sp4msh

pl.sel <- SpatialPolygons(list(Polygons(list(Polygon(
  cbind(c(-945, -920, -880, -880, -940, -965, -945), 
        c(2880, 2856, 2856, 2900, 2987, 2985, 2880)),
  FALSE)), '0')), proj4string = crs_km)
plot(pl.sel)

spatr_bathym <- terra::rast(bathym_longlat)
# SpatVector
spatvect_bathym <- terra::as.polygons(spatr_bathym > -Inf) # bathym_poly
#SpatialPolygonsDataFrame 
spdf.bathym <- as(spatvect_bathym, "Spatial") # bathym_sp_df

sp4perimeter <- spTransform(spdf.bathym, crs_km)
poly.water.book <- gDifference(pl.sel, sp4perimeter)

plot(pl.sel, col = alpha("skyblue", 0.5), asp = 1)
plot(sp4perimeter, add = TRUE, col = alpha(gray(0.9), 0.5))

plot(pl.sel, col = alpha("skyblue", 0.5), asp = 1)
plot(sp4perimeter, add = TRUE, col = alpha(gray(0.9), 0.5))
plot(poly.water.book, col = alpha("skyblue", 0.5))

plot(pl.sel, col = alpha("skyblue", 0.5), asp = 1)
plot(sp4msh, add = TRUE, col = alpha(gray(0.9), 0.5))

spatr.b1 <- terra::rast(raster_b1)
spatv.b1 <- terra::as.polygons(spatr.b1 > -Inf) # bathym_poly
spdf.b1_ <- as(spatv.b1, "Spatial") 
spdf.b1 <- spTransform(spdf.b1_, crs_km)

sp.b1_ <- geometry(spdf.b1)
n.spb1 <- length(sp.b1_@polygons[[1]]@Polygons)
idx.spb1 <- seq(1:n.spb1)
sea_list = lapply(idx.spb1, function(n) Polygon(sp.b1_@polygons[[1]]@Polygons[[n]]@coords, hole = F)) 

sp.b1 <- SpatialPolygons(list(Polygons(sea_list, ID = runif(1)))) #poly1

sp.b1.ex <- SpatialPolygons(c(sp.b1@polygons, poly.water.book@polygons))

max.edge.length <- 1 #1
bound.outer <- diff(range(cords_Dng.xy[,1][-c(5, 19)]))/3
bound.outer2 <- diff(range(cords_Dng.xy[,1][-c(5, 19)]))
mesh.dng <- inla.mesh.2d(boundary = sp4msh,
                         max.edge = c(1,5)*max.edge.length,
                         cutoff = 0.05,
                         offset = c(max.edge.length, bound.outer2))
#plot(mesh.dng)
mesh.dng$crs <- crs_km

gg_mesh <- ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='purple',size=1.7,alpha=0.5)


# Spatial polygons for barrier 2, between -10 and -500, poly2
spatr.b2 <- terra::rast(raster_b2)
spatv.b2 <- terra::as.polygons(spatr.b2 > -Inf) # bathym_poly
spdf.b2_ <- as(spatv.b2, "Spatial") 
spdf.b2 <- spTransform(spdf.b2_, crs_km)

sp.b2_ <- geometry(spdf.b2)
n.spb2 <- length(sp.b2_@polygons[[1]]@Polygons)
idx.spb2 <- seq(1:n.spb2)
sea_list2 = lapply(idx.spb2, function(n) Polygon(sp.b2_@polygons[[1]]@Polygons[[n]]@coords, hole = F)) 
sp.b2 <- SpatialPolygons(list(Polygons(sea_list2, ID = runif(1)))) #poly2

# WITH spb1 <- spb1.ex
sp.b1 <- sp.b1.ex
sp.bars <- SpatialPolygons(c(sp.b1@polygons, sp.b2@polygons)) #poly.original

# BARRIER TRIANGLES
tl <- length(mesh.dng$graph$tv[,1])
# - the number of triangles in the mesh.pp
posTri <- matrix(0, tl, 2)

for (t in 1:tl){
  temp = mesh.dng$loc[mesh.dng$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}

posTri <- SpatialPoints(posTri)

# - the positions of the triangle centers
bars.centers <- over(sp.bars, SpatialPoints(posTri), returnList=T)
# - checking which mesh.dng triangles are inside the barrier area
bars.centers <- unlist(bars.centers) #bar.original
poly.bars <- inla.barrier.polygon(mesh.dng, barrier.triangles = bars.centers)
# - the Barrier model's polygon
# - in most cases this should be the same as poly.original

# BARRIER 1
bar1 <- over(sp.b1, SpatialPoints(posTri), returnList=T)
# - checking which mesh.dng triangles are inside the barrier area
bar1 <- unlist(bar1)
poly.bar1 <- inla.barrier.polygon(mesh.dng, barrier.triangles = bar1)

# BARRIER 2
bar2 <- over(sp.b2, SpatialPoints(posTri), returnList=T)
# - checking which mesh.dng triangles are inside the barrier area
bar2 <- unlist(bar2)
poly.bar2 <- inla.barrier.polygon(mesh.dng, barrier.triangles = bar2)

fem <-  inla.barrier.fem.plus(mesh.dng, list(bar1, bar2))

barrier.triangles <- list(bar1, bar2)

plot(mesh.dng, main="Mesh and Omega")
plot(poly.bars, add=T, col='lightblue')


######## getting rid of the holes

poly.water <- spdf.water
zlim = c(0.1, 1)
xlim = poly.water@bbox[1, ] 
ylim = poly.water@bbox[2, ]

proj = inla.mesh.projector(mesh.dng, xlim = xlim, 
                           ylim = ylim, dims=c(300, 300))

# Plots for the mesh
poly.water_sf <- st_as_sf(poly.water)
poly.bars_sf <- st_as_sf(poly.bars)
poly.bar1_sf <- st_as_sf(poly.bar1)
poly.bar2_sf <- st_as_sf(poly.bar2)

st_crs(poly.water_sf) <- crs_km
st_crs(poly.bars_sf) <- crs_km
st_crs(poly.bar1_sf) <- crs_km
st_crs(poly.bar2_sf) <- crs_km

gg_mesh + 
  geom_sf(data = poly.water_sf,
          col='red', alpha=0.5) 
# covering barrier 1, islands and open deep sea
gg_mesh + 
  geom_sf(data = poly.bar1_sf,
          col='red', alpha=0.5) 
# covering barrier 2, water < -20 (and > -500)
gg_mesh + 
  geom_sf(data = poly.bar2_sf,
          col='red', alpha=0.5) 
# covering barrier 1 and barrier 2
gg_mesh + 
  geom_sf(data = poly.bars_sf,
          col='red', alpha=0.5) 

# covering barrier 1, islands and open deep sea
gg_mesh.nob1 <- 
  gg_mesh + 
  geom_sf(data = poly.bar1_sf, fill = "white")

# covering barrier 1 + partially covering bar 2
gg_mesh.bars <-
  gg_mesh.nob1 +
  geom_sf(data = poly.bar2_sf,
          col='deeppink3', fill = "pink", alpha=0.2)

location <- matrix(c(c(-901.5532), 
                     c(2862.0542)), ncol = 2)
location.sp <- SpatialPoints(location, proj4string = crs_km)

gg_mesh.bars +
  geom_sf(data = st_as_sf(location.sp), col = "red") +
  xlim(c(-910, -890)) +
  ylim(c(2855,2870))

gg_mesh.bars +
  xlim(c(-950, -900)) +
  ylim(c(2870,2890))

poly.bars@polygons
poly.water@polygons
sp4msh

#######
#######
#trying with -50

dmesh <- book.mesh.dual(mesh.dng)

domainSP <- sp4msh
domainSPsf <- st_as_sf(domainSP)
dmesh_sf <- st_as_sf(dmesh) 
st_crs(dmesh_sf) <- crs_km

# with sapply
w.bath50 <- sapply(1:length(dmesh), function(i) {
  if (length(st_intersects(dmesh_sf[i, ], domainSPsf)[[1]]) == 1)
    return(st_area(st_intersection(dmesh_sf[i, ], domainSPsf)))
  else {
    return(0)
  }
})

# Saving on object in RData format
#save(w.bath50, file = "w's/w.bath50.dng.RData")
# Save multiple objects
#save(data1, data2, file = "data.RData")
# To load the data again
#load("data.RData")
save.image()

sum(w.bath50)
table(w.bath50>0); table(w.bath50==0)
pal <- wes_palette("Zissou1")
colr = rep(c(pal[2]), length = as.numeric(length(w.bath50)))
colr[w.bath50>0] = pal[5]
plot(dmesh, col = colr)

#####
#####
#projection matrices so I can checl plot.field

w <- w.bath50
n <- nrow(cords_Dng.xy)
nv <- mesh.dng$n
xy <- cords_Dng.xy
xy.sp <- cords_Dng.sp

y.pp <- rep(0:1, c(nv, n))
#The exposure vector can be defined as:
e.pp <- c(w, rep(0, n)) 
length(y.pp); length(e.pp)
# The projection matrix is defined in two steps. For the integration points this is just a diagonal matrix because these locations are just the mesh vertices:
imat <- Diagonal(nv, rep(1, nv))
# For the observed points, another projection matrix is defined:
lmat <- inla.spde.make.A(mesh.dng, xy.sp)
# The entire projection matrix is:
A.pp <- rbind(imat, lmat)

# We set up the data stack as follows:
stk.pp <- inla.stack(
  data = list(y = y.pp, e = e.pp), 
  A = list(1, A.pp),
  effects = list(list(b0 = rep(1, nv + n)), list(i = 1:nv)),
  tag = 'pp')


A.pred <- inla.spde.make.A(mesh.dng, mesh.dng$loc[, 1:2])

#Stack for prediction at mesh nodes
stk.pred <- inla.stack(
  data = list(y = NA, e = 0),
  A = list(A.pred, 1),
  effects =list(list(i = 1:mesh.dng$n), list(b0 = rep(1, nrow(A.pred)))), #data = list(y = rep(NA, nrow(coop)), e = rep(0, nrow(coop))),
  tag = 'pred')

joint.stk <- inla.stack(stk.pp, stk.pred)
idx.pred <- inla.stack.index(joint.stk, 'pred')$data
#####
#####
#test plot.field
#from [[24]]

x <- matrix(c(1,1, 0.01,0.01, 0.01,1, 0.01,0.2, 0.01,0.3, 0.01,0.5, 0.01,0.7, 0.01,0.8), ncol = 2, byrow = TRUE)

for (i in 1:nrow(x)) {
  tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                    fem = fem, 
                                    barrier.triangles = barrier.triangles, 
                                    prior.range = prior.range[1,], #21.0  0.5
                                    prior.sigma = prior.sigma, #1.0 0.1
                                    range.fraction = c(x[i,1],x[i,2]))
  model[[i+23]] <- tbm
  formula[[i+23]] <- y ~ 0 + b0 + f(i, model =  tbm)
  
  res.dng[[i+23]] <- inla(formula[[i+23]],
                       data = inla.stack.data(joint.stk),
                       family = 'poisson', 
                       control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                       E = inla.stack.data(joint.stk)$e, 
                       inla.mode = "experimental")
  
  res.dng.pred[[i+23]] <- inla(formula[[i+23]],
                            data = inla.stack.data(joint.stk),
                            family = 'poisson', 
                            control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                            E = inla.stack.data(joint.stk)$e, 
                            inla.mode = "experimental")
}

local.plot.field.book <- function(field, xlim, ylim, ...){
  if (missing(xlim)) xlim = poly.water@bbox[1, ] 
  if (missing(ylim)) ylim = poly.water@bbox[2, ]
  proj = inla.mesh.projector(mesh, xlim = xlim,
                             ylim = ylim, dims=c(300, 300))
  field.proj = inla.mesh.project(proj, field)
  image.plot(list(x = proj$x, y = proj$y, z = field.proj),
             xlim = xlim, ylim = ylim, ...)
}

# [24] 1,1 [25] 0.01,0.01 [26] 0.01,1 [27] 0.01,0.2 [28] 0.01,0.3 [29] 0.01,0.5 [30] 0.01,0.7 [31] 0.01,0.8
# so take 24 as the refernce cuz the mean shouldn't be (a lot) higher than without barriers
mesh <- mesh.dng

local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean) #, zlim = c(0, 10))
plot(poly.bar1, add = TRUE, col = alpha("snow3", 0.5), lty = 0)

plot(poly.bar2, add = TRUE, col = alpha("pink", 0.5), lty = 0) #alpha("pink", 0.8)
#plot(poly.water.book, add = T, col = "white", lty = 0)

local.plot.field.book(
  res.dng[[31]]$summary.random$i$mean)
plot(poly.bar1, add = TRUE, col = alpha("snow3", 0.5), lty = 0)

plot(poly.bar2, add = TRUE, col = alpha("pink", 0.5), lty = 0) #alpha("pink", 0.8)

local.plot.field.book(
  res.dng[[31]]$summary.random$i$sd)#, zlim = c(5, 15))
plot(poly.bar1, add = TRUE, col = alpha("snow3", 0.5), lty = 0)

plot(poly.bar2, add = TRUE, col = alpha("pink", 0.5), lty = 0) #alpha("pink", 0.8)

#close the polygons with really high values

local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = c(-950, -910), ylim = c(2860, 2895))

######
######
#identifying and coloring the normal triangles green 

# SpatRaster
spatr_no.barriers <- terra::rast(raster_no.b2)
# SpatVector
spatvect_no.barriers <- terra::as.polygons(spatr_no.barriers > -Inf) # bathym_poly
#SpatialPolygonsDataFrame 
spdf.no.bars <- as(spatvect_no.barriers, "Spatial") # bathym_sp_df
# for the mesh I need b2 without NA
spdf.no.bars <- spTransform(spdf.no.bars, crs_km)

spdf.no.bars_ <- geometry(spdf.no.bars)
n.spdf.no.bars <- length(spdf.no.bars_@polygons[[1]]@Polygons)
idx.spdf.no.bars <- seq(1:n.spdf.no.bars)
no.bars_list = lapply(idx.spdf.no.bars, function(n) Polygon(spdf.no.bars_@polygons[[1]]@Polygons[[n]]@coords, hole = F)) 
sp.no.bars <- SpatialPolygons(list(Polygons(no.bars_list, ID = runif(1)))) 
#####

sea.centers <- over(sp.no.bars, SpatialPoints(posTri), returnList=T)
# - checking which mesh.dng triangles are inside the normal area
sea.centers <- unlist(sea.centers) 
#poly.seas <- inla.barrier.polygon(mesh.dng, barrier.triangles = sea.centers)

# - the number of triangles in the mesh.pp
posTri_ <- matrix(0, tl, 2)

for (t in 1:tl){
  temp = mesh.dng$loc[mesh.dng$graph$tv[t, ], ]
  posTri_[t,] = colMeans(temp)[c(1,2)] 
}

id.tri <- intersect(1:tl, sea.centers) #all triangles in normal area
coord.sea.center <- posTri_[id.tri,]
loc.sea.center.sp <- SpatialPoints(coord.sea.center, proj4string = crs_km)

par(mfrow = c(2,2))

gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(c(-950, -910)) +
  ylim(c(2860,2895))

local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = c(-950, -910), ylim = c(2860, 2895))

#using coord.sea.center identify the triangles that are in the middle of barrier and are inflating results
#-942, ~2280, there's and island, then normal sea, and then deep sea outside. 
xlim_ <- c(-943, -941)
ylim_ <- c(2879, 2881)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)

local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

#coord sea centers I need to move from normal to bars, these I have to move for bar2
csc <- coord.sea.center[which(coord.sea.center[,1] < -941 & coord.sea.center[,1] > -943),]
csc <- csc[which(csc[,2] < 2881 & csc[,2] > 2879),]

# I don't need to be super precise with the coordinates 
# as long as I am getting the triangles of normal area I want ot get rid off
# I also need to find these islands in the middle, 
# but let's see if they are causing problems or not doing the model

id.csc2.list <- list()
bar2.new.list <- list()
bar2.new.id.list <- list()

id.csc <- match(csc[,1], posTri_[,1]) #if I had more than 26 i would have to filter by [,2] too
# id of triangles I need to move
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
intersect(bar2.new$'0.376775514567271', id.csc)
bar2.new$'0.376775514567271' <- c(bar2.new$'0.376775514567271', id.csc)
bar2.new.id <- bar2.new$'0.376775514567271'
bar2.new <- unlist(bar2.new)

#agregar [[1]]
id.csc2.list[[1]] <- id.csc
bar2.new.list[[1]] <- bar2.new
bar2.new.id.list[[1]] <- bar2.new.id

#check it works
coord.bar2.new <- posTri_[bar2.new.id,]
loc.bar2.new.sp <- SpatialPoints(coord.bar2.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar2.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#do the model with the new bar2

fem <-  inla.barrier.fem.plus(mesh.dng, list(bar1, bar2.new))
barrier.triangles <- list(bar1, bar2.new)

#x <- matrix(c(1,1, 0.01,0.01, 0.01,1, 0.01,0.2, 0.01,0.3, 0.01,0.5, 0.01,0.7, 0.01,0.8), ncol = 2, byrow = TRUE)

for (i in 2:3) {
  tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                    fem = fem, 
                                    barrier.triangles = barrier.triangles, 
                                    prior.range = prior.range[1,], #21.0  0.5
                                    prior.sigma = prior.sigma, #1.0 0.1
                                    range.fraction = c(x[i,1],x[i,2]))
  model[[i+31]] <- tbm
  formula[[i+31]] <- y ~ 0 + b0 + f(i, model =  tbm)
  
  res.dng[[i+31]] <- inla(formula[[i+31]],
                          data = inla.stack.data(joint.stk),
                          family = 'poisson', 
                          control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                          E = inla.stack.data(joint.stk)$e, 
                          inla.mode = "experimental")
  
  res.dng.pred[[i+31]] <- inla(formula[[i+31]],
                               data = inla.stack.data(joint.stk),
                               family = 'poisson', 
                               control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                               E = inla.stack.data(joint.stk)$e, 
                               inla.mode = "experimental")
}

local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = c(-950, -910), ylim = c(2860, 2895))

local.plot.field.book(
  res.dng[[34]]$summary.random$i$mean, 
  xlim = c(-950, -910), ylim = c(2860, 2895))

local.plot.field.book(
  res.dng[[34]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

##it looks fine so I have to do it for the others
local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean)

xlim_ <- c(-933, -931)
ylim_ <- c(2888.5, 2890)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)

local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

#coord sea centers I need to move from normal to bars, these I have to move for bar1
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[1], posTri_[,1]) #csc[1] only cuz there's one tri, otherwise csc[,1]
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)
intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list <- list()
bar1.new.list <- list()
bar1.new.id.list <- list()
#agregar [[1]]
id.csc1.list[[1]] <- id.csc
bar1.new.list[[1]] <- bar1.new
bar1.new.id.list[[1]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#do the model with the new bar1

fem <-  inla.barrier.fem.plus(mesh.dng, list(bar1.new, bar2.new))
barrier.triangles <- list(bar1.new, bar2.new)

#x <- matrix(c(1,1, 0.01,0.01, 0.01,1, 0.01,0.2, 0.01,0.3, 0.01,0.5, 0.01,0.7, 0.01,0.8), ncol = 2, byrow = TRUE)

  tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                    fem = fem, 
                                    barrier.triangles = barrier.triangles, 
                                    prior.range = prior.range[1,], #21.0  0.5
                                    prior.sigma = prior.sigma, #1.0 0.1
                                    range.fraction = c(x[3,1],x[3,2]))
  model[[35]] <- tbm
  formula[[35]] <- y ~ 0 + b0 + f(i, model =  tbm)
  
  res.dng[[35]] <- inla(formula[[35]],
                          data = inla.stack.data(joint.stk),
                          family = 'poisson', 
                          control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                          E = inla.stack.data(joint.stk)$e, 
                          inla.mode = "experimental")
  
  res.dng.pred[[35]] <- inla(formula[[35]],
                               data = inla.stack.data(joint.stk),
                               family = 'poisson', 
                               control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                               E = inla.stack.data(joint.stk)$e, 
                               inla.mode = "experimental")

local.plot.field.book(
    res.dng[[35]]$summary.random$i$sd, 
    xlim = c(-950, -910), ylim = c(2860, 2895))
local.plot.field.book(
  res.dng[[25]]$summary.random$i$sd, 
  xlim = c(-950, -910), ylim = c(2860, 2895))

#####
##### try to do all triangles
#####

par(mfrow = c(2,2))
local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean)
xlim_ <- c(-975, -962)
ylim_ <- c(2950, 2962)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

#coord sea centers I need to move from normal to bars, these I have to move for bar1
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
intersect(bar2.new$'0.376775514567271', id.csc)
bar2.new$'0.376775514567271' <- c(bar2.new$'0.376775514567271', id.csc)
bar2.new.id <- bar2.new$'0.376775514567271'
bar2.new <- unlist(bar2.new)

id.csc2.list[[2]] <- id.csc
bar2.new.list[[2]] <- bar2.new
bar2.new.id.list[[2]] <- bar2.new.id

#check it works
coord.bar2.new <- posTri_[bar2.new.id,]
loc.bar2.new.sp <- SpatialPoints(coord.bar2.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar2.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#######
#######

local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean)
xlim_ <- c(-964.6, -963.8)
ylim_ <- c(2919, 2920.05)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)

local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

#coord sea centers I need to move from normal to bars, these I have to move for bar1
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[2]] <- id.csc
bar1.new.list[[2]] <- bar1.new
bar1.new.id.list[[2]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#######
#######
xlim_ <- c(-962.6, -962)
ylim_ <- c(c(2921, 2921.6))
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[3]] <- id.csc
bar1.new.list[[3]] <- bar1.new
bar1.new.id.list[[3]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#######
#######
xlim_ <- c(-957, -956.2)
ylim_ <- c(2925, 2927)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[4]] <- id.csc
bar1.new.list[[4]] <- bar1.new
bar1.new.id.list[[4]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#######
#######
xlim_ <- c(-959.75, -959.3)
ylim_ <- c(2933.8, 2934.25)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[5]] <- id.csc
bar1.new.list[[5]] <- bar1.new
bar1.new.id.list[[5]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)


#######
#######
xlim_ <- c(-960.1, -959.5)
ylim_ <- c(2932.4, 2933.05)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[6]] <- id.csc
bar1.new.list[[6]] <- bar1.new
bar1.new.id.list[[6]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#######
#######
xlim_ <- c(-961.5, -959.2)
ylim_ <- c(2926.6, 2930.8)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)

#coord sea centers I need to move from normal to bars, these I have to move for bar1
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
intersect(bar2.new$'0.376775514567271', id.csc)
bar2.new$'0.376775514567271' <- c(bar2.new$'0.376775514567271', id.csc)
bar2.new.id <- bar2.new$'0.376775514567271'
bar2.new <- unlist(bar2.new)

id.csc2.list[[3]] <- id.csc
bar2.new.list[[3]] <- bar2.new
bar2.new.id.list[[3]] <- bar2.new.id

#check it works
coord.bar2.new <- posTri_[bar2.new.id,]
loc.bar2.new.sp <- SpatialPoints(coord.bar2.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar2.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#######
#######
xlim_ <- c(-942.485, -941)
ylim_ <- c(2929.95, 2932)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
#I need to divided in 2
#first half
xlim_ <- c(-942.485, -941)
ylim_ <- c(2930.45, 2932)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[7]] <- id.csc
bar1.new.list[[7]] <- bar1.new
bar1.new.id.list[[7]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#second half
xlim_ <- c(-942.2, -941)
ylim_ <- c(2929.93, 2930.365)

gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[8]] <- id.csc
bar1.new.list[[8]] <- bar1.new
bar1.new.id.list[[8]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#######
#######
xlim_ <- c(-941.6, -940.5)
ylim_ <- c(2946, 2947.6)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[25]]$summary.random$i$mean, 
  xlim = xlim_, ylim = ylim_)
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[9]] <- id.csc
bar1.new.list[[9]] <- bar1.new
bar1.new.id.list[[9]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_) #ylim(c(2946, 2948))

#######
####### do the model
#######

#bar 1 total
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)
id.csc1 <- unlist(id.csc1.list)
intersect(bar1.new$'0.504909630864859', id.csc1)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc1)
bar1.new.id <- bar1.new$'0.504909630864859'
bar1.new <- unlist(bar1.new)

#bar 2 total
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
id.csc2 <- unlist(id.csc2.list)
intersect(bar2.new$'0.376775514567271', id.csc2)
bar2.new$'0.376775514567271' <- c(bar2.new$'0.376775514567271', id.csc2)
bar2.new.id <- bar2.new$'0.376775514567271'
bar2.new <- unlist(bar2.new)

fem <-  inla.barrier.fem.plus(mesh.dng, list(bar1.new, bar2.new))
barrier.triangles <- list(bar1.new, bar2.new)

#x <- matrix(c(1,1, 0.01,0.01, 0.01,1, 0.01,0.2, 0.01,0.3, 0.01,0.5, 0.01,0.7, 0.01,0.8), ncol = 2, byrow = TRUE)
tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                  fem = fem, 
                                  barrier.triangles = barrier.triangles, 
                                  prior.range = prior.range[1,], #21.0  0.5
                                  prior.sigma = prior.sigma, #1.0 0.1
                                  range.fraction = c(x[3,1],x[3,2]))
model[[36]] <- tbm
formula[[36]] <- y ~ 0 + b0 + f(i, model =  tbm)

res.dng[[36]] <- inla(formula[[36]],
                      data = inla.stack.data(joint.stk),
                      family = 'poisson', 
                      control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                      E = inla.stack.data(joint.stk)$e, 
                      inla.mode = "experimental")

res.dng.pred[[36]] <- inla(formula[[36]],
                           data = inla.stack.data(joint.stk),
                           family = 'poisson', 
                           control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                           E = inla.stack.data(joint.stk)$e, 
                           inla.mode = "experimental")

par(mfrow=c(1,1))
local.plot.field.book(
  res.dng[[36]]$summary.random$i$sd)
local.plot.field.book(
  res.dng[[36]]$summary.random$i$sd, zlim = c(50, 200), 
  xlim = c(-965, -935), ylim = c(2915, 2935))

#######
#######
xlim_ <- c(-964.5, -963.845)
ylim_ <- c(2919.5, 2920.3)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[36]]$summary.random$i$sd,
  xlim = xlim_, ylim = ylim_)

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[10]] <- id.csc
bar1.new.list[[10]] <- bar1.new
bar1.new.id.list[[10]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_) 

#######
#######
xlim_ <- c(-957, -956.15)
ylim_ <- c(2925, 2927)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[36]]$summary.random$i$sd,
  xlim = xlim_, ylim = ylim_)

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[11]] <- id.csc
bar1.new.list[[11]] <- bar1.new
bar1.new.id.list[[11]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_) 

#######
#######
par(mfrow=c(2,2))
xlim_ <- c(-942.8, -941)
ylim_ <- c(2929.8, 2932)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
local.plot.field.book(
  res.dng[[36]]$summary.random$i$sd,
  xlim = xlim_, ylim = ylim_)
#I need to divided in 2
#first half
xlim_ <- c(-942.8, -941)
ylim_ <- c(2930.4, 2932)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[12]] <- id.csc
bar1.new.list[[12]] <- bar1.new
bar1.new.id.list[[12]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)
#second half
xlim_ <- c(-942.2, -941)
ylim_ <- c(2929.85, 2930.4)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.504909630864859', id.csc)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc)
bar1.new.id <- bar1.new$'0.504909630864859'

id.csc1.list[[13]] <- id.csc
bar1.new.list[[13]] <- bar1.new
bar1.new.id.list[[13]] <- bar1.new.id

bar1.new <- unlist(bar1.new)
#check it works
coord.bar1.new <- posTri_[bar1.new.id,]
loc.bar1.new.sp <- SpatialPoints(coord.bar1.new, proj4string = crs_km)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.bar1.new.sp), col = "orange") +
  xlim(xlim_) +
  ylim(ylim_)

#######
####### do the model again
#######

#bar 1 total
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)
id.csc1 <- unlist(id.csc1.list)
intersect(bar1.new$'0.504909630864859', id.csc1)
bar1.new$'0.504909630864859' <- c(bar1.new$'0.504909630864859', id.csc1)
bar1.new.id <- bar1.new$'0.504909630864859'
bar1.new <- unlist(bar1.new)

#bar 2 total
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
id.csc2 <- unlist(id.csc2.list)
intersect(bar2.new$'0.376775514567271', id.csc2)
bar2.new$'0.376775514567271' <- c(bar2.new$'0.376775514567271', id.csc2)
bar2.new.id <- bar2.new$'0.376775514567271'
bar2.new <- unlist(bar2.new)

fem <-  inla.barrier.fem.plus(mesh.dng, list(bar1.new, bar2.new))
barrier.triangles <- list(bar1.new, bar2.new)

#x <- matrix(c(1,1, 0.01,0.01, 0.01,1, 0.01,0.2, 0.01,0.3, 0.01,0.5, 0.01,0.7, 0.01,0.8), ncol = 2, byrow = TRUE)
tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                  fem = fem, 
                                  barrier.triangles = barrier.triangles, 
                                  prior.range = prior.range[1,], #21.0  0.5
                                  prior.sigma = prior.sigma, #1.0 0.1
                                  range.fraction = c(x[3,1],x[3,2]))
model[[37]] <- tbm
formula[[37]] <- y ~ 0 + b0 + f(i, model =  tbm)

res.dng[[37]] <- inla(formula[[37]],
                      data = inla.stack.data(joint.stk),
                      family = 'poisson', 
                      control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                      E = inla.stack.data(joint.stk)$e, 
                      inla.mode = "experimental")

res.dng.pred[[37]] <- inla(formula[[37]],
                           data = inla.stack.data(joint.stk),
                           family = 'poisson', 
                           control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                           E = inla.stack.data(joint.stk)$e, 
                           inla.mode = "experimental")

par(mfrow=c(1,1))
local.plot.field.book(
  res.dng[[37]]$summary.random$i$sd)
local.plot.field.book(
  res.dng[[37]]$summary.random$i$mean)

#######
####### for more range fractions
#######

x <- matrix(c(1,1, 0.01,0.01, 0.01,1, 0.01,0.2, 0.01,0.3, 0.01,0.5, 0.01,0.7, 0.01,0.8), ncol = 2, byrow = TRUE)

for (k in 1:nrow(x)) {
  tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                    fem = fem, 
                                    barrier.triangles = barrier.triangles, 
                                    prior.range = prior.range[1,], #21.0  0.5
                                    prior.sigma = prior.sigma, #1.0 0.1
                                    range.fraction = c(x[k,1],x[k,2]))
  model[[k+37]] <- tbm
  formula[[k+37]] <- y ~ 0 + b0 + f(i, model =  tbm)
  
  res.dng[[k+37]] <- inla(formula[[k+37]],
                          data = inla.stack.data(joint.stk),
                          family = 'poisson', 
                          control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                          E = inla.stack.data(joint.stk)$e, 
                          inla.mode = "experimental")
  
  res.dng.pred[[k+37]] <- inla(formula[[k+37]],
                               data = inla.stack.data(joint.stk),
                               family = 'poisson', 
                               control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                               E = inla.stack.data(joint.stk)$e, 
                               inla.mode = "experimental")
}

# [24] 1,1 [25] 0.01,0.01 [26] 0.01,1 [27] 0.01,0.2 [28] 0.01,0.3 [29] 0.01,0.5 [30] 0.01,0.7 [31] 0.01,0.8

#without any fixes
par(mfrow=c(2,2))
local.plot.field.book(
  res.dng[[24]]$summary.random$i$sd)
local.plot.field.book(
  res.dng[[25]]$summary.random$i$sd)
local.plot.field.book(
  res.dng[[26]]$summary.random$i$sd)
local.plot.field.book(
  res.dng[[29]]$summary.random$i$sd)

#first fixes for 0.01,1
par(mfrow = c(2,2))
local.plot.field.book(
  res.dng[[36]]$summary.random$i$sd)
#second fixes based on results from [[36]]
#for c(1,1)
local.plot.field.book( 
  res.dng[[37]]$summary.random$i$sd)
#for c(0.01,1)
local.plot.field.book(
  res.dng[[38]]$summary.random$i$sd)
#for c(0.01,1)
local.plot.field.book( 
  res.dng[[39]]$summary.random$i$sd)





# to get rid of the normal tri, take the out of posTri

# let's get rid of these triangles and see what happens with the model.
# maybe the reason. why it worked with both barriers being the same it's because I have 3 ranges too close?




# I was trying to identigy the triangles, but it might be easier to just change the bath values to NA




##########

# this identifies nodes not centers
# make rule to get rid of identify triangle jumps lower than 2 for example?
id.coord39 <- c(mesh.dng$loc[39, 1], mesh.dng$loc[39, 2])
loc39 <- matrix(id.coord39, ncol = 2, byrow = F)
loc39.sp <- SpatialPoints(loc39, proj4string = crs_km)

id.tri <- c(39, 40)
id.coord.mesh <- c(mesh.dng$loc[id.tri, 1], mesh.dng$loc[id.tri, 2])
locmesh <- matrix(id.coord.mesh, ncol = 2, byrow = F)
locmesh.sp <- SpatialPoints(locmesh, proj4string = crs_km)


gg_mesh.bars +
  geom_sf(data = st_as_sf(loc39.sp), col = "red") +
  xlim(c(-962, -956)) +
  ylim(c(2930,2935))

gg_mesh.bars +
  geom_sf(data = st_as_sf(locmesh.sp), col = "red") +
  xlim(c(-962, -956)) +
  ylim(c(2930,2935))

id.tri <- list()
id.coord.mesh <- list()
locmesh <- list()
locmesh.sp <- list()

# these are id nodes for triangles in water (normal area) and barriers
id.tri[[1]] <- c(39,43,43,44,45,47,48,49,50,51,53,57) #in water
id.tri[[2]] <- c(36,37,38,40,42,46,53,54,55,56) #in barrier

id.coord.mesh[[1]] <- c(mesh.dng$loc[id.tri[[1]], 1], mesh.dng$loc[id.tri[[1]], 2])
locmesh[[1]] <- matrix(id.coord.mesh[[1]], ncol = 2, byrow = F)
locmesh.sp[[1]] <- SpatialPoints(locmesh[[1]], proj4string = crs_km)

id.coord.mesh[[2]] <- c(mesh.dng$loc[id.tri[[2]], 1], mesh.dng$loc[id.tri[[2]], 2])
locmesh[[2]] <- matrix(id.coord.mesh[[2]], ncol = 2, byrow = F)
locmesh.sp[[2]] <- SpatialPoints(locmesh[[2]], proj4string = crs_km)

gg_mesh.bars +
  geom_sf(data = st_as_sf(locmesh.sp[[1]]), col = "green") +
  geom_sf(data = st_as_sf(locmesh.sp[[2]]), col = "red") +
  xlim(c(-962, -956)) +
  ylim(c(2930,2935))



##########
# After figuring out the holes...
###########






#