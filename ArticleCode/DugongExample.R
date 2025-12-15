#fp <- file.path("ArticleCode")
#if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)

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
set.seed(17)
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
set.seed(602)
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

###########################################################################################################
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
#load("w's/w.bath50.dng.RData")

sum(w.bath50)
table(w.bath50>0); table(w.bath50==0)
pal <- wes_palette("Zissou1")
colr = rep(c(pal[2]), length = as.numeric(length(w.bath50)))
colr[w.bath50>0] = pal[5]
plot(dmesh, col = colr)

#####
#####
#projection matrices so I can check plot.field

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

###########################################################################################################
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
set.seed(880)
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

gg_mesh.nob1 <- 
  gg_mesh + 
  geom_sf(data = poly.bar1_sf, fill = "white")

# covering barrier 1 + partially covering bar 2
gg_mesh.bars <-
  gg_mesh.nob1 +
  geom_sf(data = poly.bar2_sf,
          col='deeppink3', fill = "pink", alpha=0.2)

gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(c(-950, -910)) +
  ylim(c(2860,2895))

#using coord.sea.center identify the triangles that are in the middle of barrier and are inflating results
#-942, ~2280, there's and island, then normal sea, and then deep sea outside. 
xlim_ <- c(-943, -941)
ylim_ <- c(2879, 2881)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)

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
intersect(bar2.new$'0.0930902801919729', id.csc) #put bar2.new$'0.0930902801919729' the number you get after bar2.new$...
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc)
bar2.new.id <- bar2.new$'0.0930902801919729'
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

xlim_ <- c(-933, -931)
ylim_ <- c(2888.5, 2890)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)

#coord sea centers I need to move from normal to bars, these I have to move for bar1
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[1], posTri_[,1]) #csc[1] only cuz there's one tri, otherwise csc[,1]
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)
intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

xlim_ <- c(-975, -962)
ylim_ <- c(2950, 2962)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)

#coord sea centers I need to move from normal to bars, these I have to move for bar1
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
intersect(bar2.new$'0.0930902801919729', id.csc)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc)
bar2.new.id <- bar2.new$'0.0930902801919729'
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

xlim_ <- c(-964.6, -963.8)
ylim_ <- c(2919, 2920.1)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)

#coord sea centers I need to move from normal to bars, these I have to move for bar1
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)
intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

#coord sea centers I need to move from normal to bars, these I have to move for bar1
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
intersect(bar2.new$'0.0930902801919729', id.csc)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc)
bar2.new.id <- bar2.new$'0.0930902801919729'
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

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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
csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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


###########################################################################################################
####### do the model
#######

#bar 1 total
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)
id.csc1 <- unlist(id.csc1.list)
intersect(bar1.new$'0.155050833243877', id.csc1)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc1)
bar1.new.id <- bar1.new$'0.155050833243877'
bar1.new <- unlist(bar1.new)

#bar 2 total
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
id.csc2 <- unlist(id.csc2.list)
intersect(bar2.new$'0.0930902801919729', id.csc2)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc2)
bar2.new.id <- bar2.new$'0.0930902801919729'
bar2.new <- unlist(bar2.new)

#######
xlim_ <- c(-964.5, -963.845)
ylim_ <- c(2919.5, 2920.3)
gg_mesh.bars +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

csc <- coord.sea.center[which(coord.sea.center[,1] < xlim_[2] & coord.sea.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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

intersect(bar1.new$'0.155050833243877', id.csc)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc)
bar1.new.id <- bar1.new$'0.155050833243877'

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


###########################################################################################################
####### do the model again

#bar 1 total
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)
id.csc1 <- unlist(id.csc1.list)
intersect(bar1.new$'0.155050833243877', id.csc1)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc1)
bar1.new.id <- bar1.new$'0.155050833243877'
bar1.new <- unlist(bar1.new)

#bar 2 total
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
id.csc2 <- unlist(id.csc2.list)
intersect(bar2.new$'0.0930902801919729', id.csc2)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc2)
bar2.new.id <- bar2.new$'0.0930902801919729'
bar2.new <- unlist(bar2.new)

##add white
pl.sel2 <- SpatialPolygons(list(Polygons(list(Polygon(
  cbind(c(-990, -945, -880, -880, -880, -940, -990, -990), 
        c(2880, 2840, 2840, 2880, 2900, 2987, 2987, 2880)),
  FALSE)), '0')), proj4string = crs_km)
plot(pl.sel2)
poly.water.book2 <- gDifference(pl.sel2, sp4perimeter)

plot(poly.water.book2, col = alpha("skyblue", 0.5))


###########################################################################################################
#plots for leaflet########
###########################################################################################################

#poly.bar1.new@proj4string <- crs_km #not necessary
###########################################################################################################
#add a third barrier
###########################################################################################################
#run again to make sure i'm using the riht barnew
#bar 1 total
bar1.new <- over(sp.b1, SpatialPoints(posTri), returnList=T)
id.csc1 <- unlist(id.csc1.list)
bar1.new$'0.155050833243877' <- c(bar1.new$'0.155050833243877', id.csc1)
#this '..' is not the same for all
bar1.new.id <- c(bar1.new$'0.155050833243877', bar1.new$'1', bar1.new$'0.155050833243877')
bar1.new <- unlist(bar1.new)

#bar 2 total
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
id.csc2 <- unlist(id.csc2.list)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc2)
#this '..' is not the same for all
bar2.new.id <- c(bar2.new$'0.0930902801919729', bar2.new$'0.0930902801919729')
bar2.new <- unlist(bar2.new)

#in principle the third barrier will be between -25 and -50 just to filter according bathymetry, 
#however I will filter later using the holes, etc
raster_no.b3 <- raster_no.b2
raster_no.b3[raster_no.b3 < -25] <- NA
# opposite
raster_b3 <- raster_no.b2
raster_b3[raster_b3 >= -25] <- NA
# plot
plot(raster_no.b3); plot(raster_b3)

# Spatial polygons for barrier 3, between -25 and -50, poly3
spatr.b3 <- terra::rast(raster_b3)
spatv.b3 <- terra::as.polygons(spatr.b3 > -Inf) # bathym_poly
spdf.b3_ <- as(spatv.b3, "Spatial") 
spdf.b3 <- spTransform(spdf.b3_, crs_km)

sp.b3_ <- geometry(spdf.b3)
n.spb3 <- length(sp.b3_@polygons[[1]]@Polygons)
idx.spb3 <- seq(1:n.spb3)
sea_list3 = lapply(idx.spb3, function(n) Polygon(sp.b3_@polygons[[1]]@Polygons[[n]]@coords, hole = F)) 
set.seed(947)
sp.b3 <- SpatialPolygons(list(Polygons(sea_list3, ID = runif(1)))) #poly2

sp.bars3 <- SpatialPolygons(c(sp.b1@polygons, sp.b2@polygons, sp.b3@polygons)) #poly.original

bars.centers3 <- over(sp.bars3, SpatialPoints(posTri), returnList=T)
# - checking which mesh.dng triangles are inside the barrier area
bars.centers3 <- unlist(bars.centers3) #bar.original
poly.bars3 <- inla.barrier.polygon(mesh.dng, barrier.triangles = bars.centers3)

# BARRIER 3
bar3 <- over(sp.b3, SpatialPoints(posTri), returnList=T)
bar3.id <- bar3
bar3.id <- bar3.id$'0.872436952777207'

#i need to use bar1.new so take out the new from bar3 in case they intersect
intersect(bar2.new.id, bar1.new.id)
intersect(bar1.new.id, bar3.id)
intersect(bar2.new.id, bar3.id)
#take out the triangles that are in both barnew and bar3 because bars cannot intersect
sum(setdiff(bar3$'0.872436952777207', bar2.new.id)) + sum(intersect(bar2.new.id, bar3.id)) == sum(bar3$'0.872436952777207')
bar3$'0.872436952777207' <- setdiff(bar3$'0.872436952777207', bar2.new.id)
#bar3 again without the tri in bar2.new
bar3.id <- bar3$'0.872436952777207'
bar3 <- unlist(bar3)
poly.bar3 <- inla.barrier.polygon(mesh.dng, barrier.triangles = bar3)

#color green the normal triangles again
# SpatRaster
spatr_no.barriers <- terra::rast(raster_no.b3)
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
set.seed(880)
sp.no.bars <- SpatialPolygons(list(Polygons(no.bars_list, ID = runif(1)))) 
#####
sea.centers <- over(sp.no.bars, SpatialPoints(posTri), returnList=T)
# - checking which mesh.dng triangles are inside the normal area
sea.centers <- unlist(sea.centers) 
#poly.seas <- inla.barrier.polygon(mesh.dng, barrier.triangles = sea.centers)

id.tri <- intersect(1:tl, sea.centers) #all triangles in normal area
#+ the new ones from barnew
id.tri3 <- setdiff(id.tri, bar1.new.id)
id.tri3 <- setdiff(id.tri3, bar2.new.id)
#normal triangle ids, 3 means it includes bar3
coord.sea.center3 <- posTri_[id.tri3,]
loc.sea.center.sp3 <- SpatialPoints(coord.sea.center3, proj4string = crs_km)

#to plot the barriers
# BARRIER 1
poly.bar1.new <- inla.barrier.polygon(mesh.dng, barrier.triangles = bar1.new)
# BARRIER 2
poly.bar2.new <- inla.barrier.polygon(mesh.dng, barrier.triangles = bar2.new)


##gg_mesh with new bar and bar3
poly.bar1_sf.new <- st_as_sf(poly.bar1.new)
poly.bar2_sf.new <- st_as_sf(poly.bar2.new)
poly.bar3_sf <- st_as_sf(poly.bar3)

st_crs(poly.bar1_sf.new) <- crs_km
st_crs(poly.bar2_sf.new) <- crs_km
st_crs(poly.bar3_sf) <- crs_km

# covering barrier 1 + partially covering bar 2
gg_mesh.bars3 <- 
  gg_mesh + 
  geom_sf(data = poly.bar1_sf.new, fill = "wheat") +
  geom_sf(data = poly.bar2_sf.new,
          col='dodgerblue4', fill = "blue", alpha=0.2) +
  geom_sf(data = poly.bar3_sf,
          col='turquoise4', fill = "turquoise1", alpha=0.2)
  

############################################################################################################
#merge barriers get rid of weird triangles
#use common sense to merge some areas
#identify bar2 centers so I can merge bar3 with bar2
coord.bar2.center <- posTri_[bar2.new.id,]
loc.bar2.center.sp <- SpatialPoints(coord.bar2.center, proj4string = crs_km)
#identify bar3 centers so I can merge bar3 with bar2
coord.bar3.center <- posTri_[bar3.id,]
loc.bar3.center.sp <- SpatialPoints(coord.bar3.center, proj4string = crs_km)

gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar2.center.sp), col = "blue") +
  xlim(c(-950, -910)) +
  ylim(c(2860,2895))

gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar3.center.sp), col = "turquoise") +
  xlim(c(-950, -910)) +
  ylim(c(2860,2895))

gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(c(-960,-900)) +
  ylim(c(2860,2960))

#merge barriers
par(mfrow=c(2,3))

xlim_ <- c(-936, -924)
ylim_ <- c(2863, 2880)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar2.center.sp), col = "blue") +
  xlim(xlim_) +
  ylim(ylim_)

#normal area and bar3 should now be bar2
#normal area first
csc <- coord.sea.center3[which(coord.sea.center3[,1] < xlim_[2] & coord.sea.center3[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc)
bar2.new.id <- bar2.new$'0.0930902801919729'

id.csc2.list3 <- list()
bar2.new.list3 <- list()
bar2.new.id.list3 <- list()

id.csc2.list3[[1]] <- id.csc
bar2.new.list3[[1]] <- bar2.new
bar2.new.id.list3[[1]] <- bar2.new.id

coord.bar2.new <- posTri_[bar2.new.id,]
loc.bar2.new.sp <- SpatialPoints(coord.bar2.new, proj4string = crs_km)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar2.new.sp), col = "blue") +
  xlim(xlim_) +
  ylim(ylim_)

#bar3 to bar2
csc <- coord.bar3.center[which(coord.bar3.center[,1] < xlim_[2] & coord.bar3.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc)
bar2.new.id <- bar2.new$'0.0930902801919729'
bar2.new <- unlist(bar2.new)

id.csc2.list3[[2]] <- id.csc
bar2.new.list3[[2]] <- bar2.new
bar2.new.id.list3[[2]] <- bar2.new.id

coord.bar2.new <- posTri_[bar2.new.id,]
loc.bar2.new.sp <- SpatialPoints(coord.bar2.new, proj4string = crs_km)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar2.new.sp), col = "blue") +
  xlim(xlim_) +
  ylim(ylim_)


########## [[3]] y [[4]]
xlim_ <- c(-923, -914)
ylim_ <- c(2864, 2872)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  xlim(xlim_) +
  ylim(ylim_)
#normal area first
csc <- coord.sea.center3[which(coord.sea.center3[,1] < xlim_[2] & coord.sea.center3[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc)
bar2.new.id <- bar2.new$'0.0930902801919729'

id.csc2.list3[[3]] <- id.csc
bar2.new.list3[[3]] <- bar2.new
bar2.new.id.list3[[3]] <- bar2.new.id

coord.bar2.new <- posTri_[bar2.new.id,]
loc.bar2.new.sp <- SpatialPoints(coord.bar2.new, proj4string = crs_km)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar2.new.sp), col = "blue") +
  xlim(xlim_) +
  ylim(ylim_)

#bar3 to bar2
csc <- coord.bar3.center[which(coord.bar3.center[,1] < xlim_[2] & coord.bar3.center[,1] > xlim_[1]),]
csc <- csc[which(csc[,2] < ylim_[2] & csc[,2] > ylim_[1]),]

id.csc <- match(csc[,1], posTri_[,1]) 
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc)
bar2.new.id <- bar2.new$'0.0930902801919729'
bar2.new <- unlist(bar2.new)

id.csc2.list3[[4]] <- id.csc
bar2.new.list3[[4]] <- bar2.new
bar2.new.id.list3[[4]] <- bar2.new.id

coord.bar2.new <- posTri_[bar2.new.id,]
loc.bar2.new.sp <- SpatialPoints(coord.bar2.new, proj4string = crs_km)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar2.new.sp), col = "blue") +
  xlim(xlim_) +
  ylim(ylim_)

####################
#done up
xlim_ <- c(-960, -900)
ylim_ <- c(2840, 2880)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar3.center.sp), col = "turquoise") +
  xlim(xlim_) +
  ylim(ylim_)

########
xlim_ <- c(-980, -900)
ylim_ <- c(2880, 2980)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar3.center.sp), col = "turquoise") +
  xlim(xlim_) +
  ylim(ylim_)

################
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.sea.center.sp), col = "green") +
  ylim(ylim_)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar2.center.sp), col = "blue") +
  ylim(ylim_)
gg_mesh.bars3 +
  geom_sf(data = st_as_sf(loc.bar3.center.sp), col = "turquoise") +
  xlim(xlim_) +
  ylim(ylim_)

###########################################################################################################
####### do the model again

#bar2 total
bar2.new <- over(sp.b2, SpatialPoints(posTri), returnList=T)
id.csc2 <- unlist(id.csc2.list3)
intersect(bar2.new$'0.0930902801919729', id.csc2)
bar2.new$'0.0930902801919729' <- c(bar2.new$'0.0930902801919729', id.csc2)
bar2.new.id <- bar2.new$'0.0930902801919729'
bar2.new <- unlist(bar2.new)

#new bar3 without bar2
bar3.new <- over(sp.b3, SpatialPoints(posTri), returnList=T)
bar3.id <- bar3.new
bar3.id <- bar3.id$'0.872436952777207'
intersect(bar3.id, bar2.new.id)
bar3.new.id <- bar3.id[!bar3.id %in% bar2.new.id] #same as:
bar3.new$'0.872436952777207' <- bar3.new$'0.872436952777207'[!bar3.id %in% bar2.new.id]
intersect(bar3.new.id, bar2.new.id)
bar3.new <- unlist(bar3.new)

#save(id.csc2.list3, file = "plots/dng/lists/id.csc2.list3.RData")
#save(id.csc1.list, file = "plots/dng/lists/id.csc1.list.RData")
#save(id.csc2.list, file = "plots/dng/lists/id.csc2.list.RData")
#save(id.csc3.list, file = "plots/dng/lists/id.csc3.list.RData")

fem <-  inla.barrier.fem.plus(mesh.dng, list(bar1.new, bar2.new, bar3.new))
barrier.triangles <- list(bar1.new, bar2.new, bar3.new)
triangle.list <- list(bar1=bar1.new, bar2=bar2.new, bar3=bar3.new, barrier.triangles=barrier.triangles)
#saveRDS(triangle.list, "plots/dng/lists/barrier.triangles.list.rds")
x3 <- matrix(c(0.01, 0.01, 0.01,
               0.01, 1, 0.01,
               0.01, 0.01, 0.5,
               0.01, 0.1, 0.5,
               0.01, 0.1, 0.8,
               0.01, 0.2, 0.5,
               0.01, 0.2, 0.8,
               0.01, 0.3, 0.5,
               0.01, 0.3, 0.8,
               0.01, 0.4, 0.5,
               0.01, 0.4, 0.8), ncol = 3, byrow = TRUE)

for (k in 1:nrow(x3)) {
  tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                    fem = fem, 
                                    barrier.triangles = barrier.triangles, 
                                    prior.range = prior.range[1,], #21.0  0.5
                                    prior.sigma = prior.sigma, #1.0 0.1
                                    range.fraction = c(x3[k,1],x3[k,2],x3[k,3]))
  model[[k+48]] <- tbm
  formula[[k+48]] <- y ~ 0 + b0 + f(i, model =  tbm)
  
  res.dng[[k+48]] <- inla(formula[[k+48]],
                          data = inla.stack.data(joint.stk),
                          family = 'poisson', 
                          control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                          E = inla.stack.data(joint.stk)$e, 
                          inla.mode = "experimental")
}

#results
#[[49]] 0.01 0.01 0.01 
#[[50]] 0.01 1.00 0.01 
#[[51]] 0.01 0.01 0.50
#[[52]] 0.01 0.10 0.50
#[[53]] 0.01 0.10 0.80
#[[54]] 0.01 0.20 0.50
#[[55]] 0.01 0.20 0.80
#[[56]] 0.01 0.30 0.50
#[[57]] 0.01 0.30 0.80
#[[58]] 0.01 0.40 0.50
#[[59]] 0.01 0.40 0.80

#more range fractions
x3 <- matrix(c(0.01, 0.01, 0.01,
               0.01, 1, 0.01,
               0.01, 0.01, 0.5,
               0.01, 0.1, 0.5,
               0.01, 0.1, 0.8,
               0.01, 0.2, 0.5,
               0.01, 0.2, 0.8,
               0.01, 0.3, 0.5,
               0.01, 0.3, 0.8,
               0.01, 0.4, 0.5,
               0.01, 0.4, 0.8), ncol = 3, byrow = TRUE)

x3.2 <- matrix(c(0.01, 0.1, 0.2,
                 0.01, 0.1, 0.7,
                 0.01, 0.1, 0.9,
                 0.01, 0.1, 1,
                 
                 0.01, 0.2, 0.2,
                 0.01, 0.2, 0.7,
                 0.01, 0.2, 0.9,
                 0.01, 0.2, 1,
                 
                 0.01, 0.3, 0.2,
                 0.01, 0.3, 0.7,
                 0.01, 0.3, 0.9,
                 0.01, 0.3, 1,
                 
                 0.01, 0.4, 0.2,
                 0.01, 0.4, 0.7,
                 0.01, 0.4, 0.9,
                 0.01, 0.4, 1,
                 
                 0.01, 0.5, 0.2,
                 0.01, 0.5, 0.7,
                 0.01, 0.5, 0.9,
                 0.01, 0.5, 1,
                 
                 0.01, 0.6, 0.2,
                 0.01, 0.6, 0.7,
                 0.01, 0.6, 0.9,
                 0.01, 0.6, 1), ncol = 3, byrow = TRUE)


for (k in 1:nrow(x3.2)) {
  tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                    fem = fem, 
                                    barrier.triangles = barrier.triangles, 
                                    prior.range = prior.range[1,], #21.0  0.5
                                    prior.sigma = prior.sigma, #1.0 0.1
                                    range.fraction = c(x3.2[k,1],x3.2[k,2],x3.2[k,3]))
  model[[k+59]] <- tbm
  formula[[k+59]] <- y ~ 0 + b0 + f(i, model =  tbm)
  
  res.dng[[k+59]] <- inla(formula[[k+59]],
                          data = inla.stack.data(joint.stk),
                          family = 'poisson', 
                          control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                          E = inla.stack.data(joint.stk)$e, 
                          inla.mode = "experimental")
}

#prev results
#[[49]] 0.01 0.01 0.01 si
#[[50]] 0.01 1.00 0.01 
#[[51]] 0.01 0.01 0.50 
#[[52]] 0.01 0.10 0.50 
#[[53]] 0.01 0.10 0.80 si


#[[54]] 0.01 0.20 0.50
#[[55]] 0.01 0.20 0.80 

#[[56]] 0.01 0.30 0.50 maybe
#[[57]] 0.01 0.30 0.80
#[[58]] 0.01 0.40 0.50 maybe
#[[59]] 0.01 0.40 0.80 maybe

#results
#[,1] [,2] [,3]
#[1,] 0.01  0.1  0.2 maybe
#[2,] 0.01  0.1  0.7
#[3,] 0.01  0.1  0.9
#[4,] 0.01  0.1  1.0

#[5,] 0.01  0.2  0.2
#[6,] 0.01  0.2  0.7
#[7,] 0.01  0.2  0.9
#[8,] 0.01  0.2  1.0

#[9,] 0.01  0.3  0.2
#[10,] 0.01  0.3  0.7
#[11,] 0.01  0.3  0.9
#[12,] 0.01  0.3  1.0

#[13,] 0.01  0.4  0.2
#[14,] 0.01  0.4  0.7 maybe
#[15,] 0.01  0.4  0.9 maybe
#[16,] 0.01  0.4  1.0

#[17,] 0.01  0.5  0.2
#[18,] 0.01  0.5  0.7
#[19,] 0.01  0.5  0.9 maybe
#[20,] 0.01  0.5  1.0

#[21,] 0.01  0.6  0.2
#[22,] 0.01  0.6  0.7 maybe
#[23,] 0.01  0.6  0.9 maybe
#[24,] 0.01  0.6  1.0

#sd
#sd


#plots with all the results and good colors
###########################################################################################################
#no close up, zlim low

zoom <- 0
zlim1 <- -1
##mean
statist <- "mean"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

x3.total <- rbind(x3, x3.2)

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  #png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}



##sd
statist <- "sd"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  #png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$sd,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.025
zlim1 <- c(-12)
statist <- "0.025"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
#fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist, "wo.bar"))

if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp
for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$'0.025quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}
#change the values lower than -8 to -8 so  they are colored blue for 0.01, [[49]]
new0.025quant <- list()
for (j in 1:nrow(x3.total)) {
  new0.025quant[[j]] <- 
    ifelse(res.dng[[48+j]]$summary.random$i$'0.025quant' < -12, -12, res.dng[[49]]$summary.random$i$'0.025quant')
  
}
#do the second 0.025 with the new lower. values = to 7.99 so they appear blue
zlim1 <- c(-12)
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
#fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist, "wo.bar"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new0.025quant[[j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 1), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$'0.975quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.5
zlim1 <- c(-12)
statist <- "0.5"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$'0.5quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##mode
zlim1 <- c(-1)
statist <- "mode"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,5.5))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#meand and mode zlim same

zlim1 <- 6
##mean
statist <- "mean"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

x3.total <- rbind(x3, x3.2)

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  #png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

statist <- "mode"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#wo colored barriers 2 and 3

zoom <- 0
zlim1 <- -1
##mean
statist <- "mean"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

x3.total <- rbind(x3, x3.2)

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  #points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##sd
statist <- "sd"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$sd,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.025
zlim1 <- c(-12)
statist <- "0.025"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
#fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist, "wo.bar"))

if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp
for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$'0.025quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  ##plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  ##plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}
#change the values lower than -8 to -8 so  they are colored blue for 0.01, [[49]]
new0.025quant <- list()
for (j in 1:nrow(x3.total)) {
  new0.025quant[[j]] <- 
    ifelse(res.dng[[48+j]]$summary.random$i$'0.025quant' < -12, -12, res.dng[[49]]$summary.random$i$'0.025quant')
  
}
#do the second 0.025 with the new lower. values = to 7.99 so they appear blue
zlim1 <- c(-12)
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
#fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist, "wo.bar"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new0.025quant[[j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 1), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$'0.975quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.5
zlim1 <- c(-12)
statist <- "0.5"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$'0.5quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##mode
zlim1 <- c(-1)
statist <- "mode"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,5.5))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#meand and mode zlim same

zlim1 <- 6
##mean
statist <- "mean"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

x3.total <- rbind(x3, x3.2)

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

statist <- "mode"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[48+j]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#sd and 0.975 correction

##sd
zlim1 <- -1
statist <- "sd"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist, "corrected"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.sd <- list()
for (j in 1:nrow(x3.total)) {
  new.sd[[48+j]] <- ifelse(res.dng[[48+j]]$summary.random$i$sd > 3.7, 3.7, res.dng[[48+j]]$summary.random$i$sd)
}


for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.sd[[48+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist, "corrected"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.0.975 <- list()
for (j in 1:nrow(x3.total)) {
  new.0.975[[48+j]] <- ifelse(res.dng[[48+j]]$summary.random$i$'0.975quant' > 12, 12, res.dng[[48+j]]$summary.random$i$'0.975quant')
}


for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.0.975[[48+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##sd
zlim1 <- -1
statist <- "sd"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist, "corrected"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.sd <- list()
for (j in 1:nrow(x3.total)) {
  new.sd[[48+j]] <- ifelse(res.dng[[48+j]]$summary.random$i$sd > 3.7, 3.7, res.dng[[48+j]]$summary.random$i$sd)
}


for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.sd[[48+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist, "corrected"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.0.975 <- list()
for (j in 1:nrow(x3.total)) {
  new.0.975[[48+j]] <- ifelse(res.dng[[48+j]]$summary.random$i$'0.975quant' > 12, 12, res.dng[[48+j]]$summary.random$i$'0.975quant')
}


for (j in 1:nrow(x3.total)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.total[j,1], "b2_", x3.total[j,2], "b3_", x3.total[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.0.975[[48+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

####################################################################################################################
# c(1,1,1); c(0.01, 1, 1)
prior.range <- c(21, 0.5)
prior.sigma <- c(1, 0.1)
x3.0 <- matrix(c(1,1,1, 0.01,1,1), ncol = 3, byrow = T) #control
for (k in 1:nrow(x3.0)) {
  tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                    fem = fem, 
                                    barrier.triangles = barrier.triangles, 
                                    prior.range = prior.range, #21.0  0.5
                                    prior.sigma = prior.sigma, #1.0 0.1
                                    range.fraction = c(x3.0[k,1],x3.0[k,2],x3.0[k,3]))
  model[[k+83]] <- tbm
  formula[[k+83]] <- y ~ 0 + b0 + f(i, model =  tbm)
  
  res.dng[[k+83]] <- inla(formula[[k+83]],
                          data = inla.stack.data(joint.stk),
                          family = 'poisson', 
                          control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                          E = inla.stack.data(joint.stk)$e, 
                          inla.mode = "experimental")
}


####################################################################################################################

zoom <- 0
zlim1 <- -1
##mean
statist <- "mean"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp


for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##sd
statist <- "sd"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$sd,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.025
zlim1 <- c(-12)
statist <- "0.025"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
#fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist, "wo.bar"))

if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp
for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$'0.025quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}
#change the values lower than -8 to -8 so  they are colored blue for 0.01, [[49]]
new0.025quant <- list()
for (j in 1:nrow(x3.0)) {
  new0.025quant[[j]] <- 
    ifelse(res.dng[[83+j]]$summary.random$i$'0.025quant' < -12, -12, res.dng[[49]]$summary.random$i$'0.025quant')
  
}
#do the second 0.025 with the new lower. values = to 7.99 so they appear blue
zlim1 <- c(-12)
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
#fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist, "wo.bar"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new0.025quant[[j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 1), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-1)
statist <- "0.975"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$'0.975quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.5
zlim1 <- c(-12)
statist <- "0.5"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$'0.5quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##mode
zlim1 <- c(-1)
statist <- "mode"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,5.5))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#meand and mode zlim same

zlim1 <- 6
##mean
statist <- "mean"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp


for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

statist <- "mode"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#wo colored barriers 2 and 3

zoom <- 0
zlim1 <- -1
##mean
statist <- "mean"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp



for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##sd
statist <- "sd"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$sd,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.025
zlim1 <- c(-12)
statist <- "0.025"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
#fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist, "wo.bar"))

if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp
for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$'0.025quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  ##plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  ##plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}
#change the values lower than -8 to -8 so  they are colored blue for 0.01, [[49]]
new0.025quant <- list()
for (j in 1:nrow(x3.0)) {
  new0.025quant[[j]] <- 
    ifelse(res.dng[[83+j]]$summary.random$i$'0.025quant' < -12, -12, res.dng[[49]]$summary.random$i$'0.025quant')
  
}
#do the second 0.025 with the new lower. values = to 7.99 so they appear blue
zlim1 <- c(-12)
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
#fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist, "wo.bar"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new0.025quant[[j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 1), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$'0.975quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.5
zlim1 <- c(-12)
statist <- "0.5"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$'0.5quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##mode
zlim1 <- c(-1)
statist <- "mode"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,5.5))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#mean and mode zlim same

zlim1 <- 6
##mean
statist <- "mean"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp



for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

statist <- "mode"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[83+j]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#sd and 0.975 correction
#color 2, w barriers
##sd
zlim1 <- -1
statist <- "sd"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist, "corrected"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  new.sd[[83+j]] <- ifelse(res.dng[[83+j]]$summary.random$i$sd > 3.7, 3.7, res.dng[[83+j]]$summary.random$i$sd)
}


for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.sd[[83+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng", paste0("zoom", zoom, "color2/zlim", zlim1, "/", statist, "corrected"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x3.0)) {
  new.0.975[[83+j]] <- ifelse(res.dng[[83+j]]$summary.random$i$'0.975quant' > 12, 12, res.dng[[83+j]]$summary.random$i$'0.975quant')
}


for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.0.975[[83+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

# wo b2 and b3
##sd
zlim1 <- -1
statist <- "sd"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist, "corrected"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.sd <- list()
for (j in 1:nrow(x3.0)) {
  new.sd[[83+j]] <- ifelse(res.dng[[83+j]]$summary.random$i$sd > 3.7, 3.7, res.dng[[83+j]]$summary.random$i$sd)
}


for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.sd[[83+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng", paste0("zoom", zoom, "wob2b3/zlim", zlim1, "/", statist, "corrected"))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.0.975 <- list()
for (j in 1:nrow(x3.0)) {
  new.0.975[[83+j]] <- ifelse(res.dng[[83+j]]$summary.random$i$'0.975quant' > 12, 12, res.dng[[83+j]]$summary.random$i$'0.975quant')
}


for (j in 1:nrow(x3.0)) {
  fpp <- file.path(paste0(path, "/", "b1_", x3.0[j,1], "b2_", x3.0[j,2], "b3_", x3.0[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.0.975[[83+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

####################################################################################################################
#Sensitivity analysis
#use different priors
0.25*diff(range(cords_Dng.xy[,2]))
0.5*diff(range(cords_Dng.xy[,2][-c(5, 19)])) #= 21.09291
0.25*diff(range(cords_Dng.xy[,2][-c(5, 19)])) #= 10.54646

prior.range <- matrix(c(21,0.5, 21,0.9, 
                        41,0.5, 41,0.9,
                        10,0.5, 10,0.9,
                        6,0.5, 6,0.9), ncol = 2, byrow = TRUE)

prior.sigma <- c(1, 0.1)
x.sensi <- rbind(x3.total[c(1,5,12,15),], x3.0)

for (j in 1:nrow(x.sensi)) { #
  model[[j+85]] <- list()
  formula[[j+85]] <- list()
  res.dng[[j+85]] <- list()
  
  for (k in 1:nrow(prior.range)) { #
    tbm <- inla.barrier.pcmatern.plus(mesh =  mesh.dng, 
                                      fem = fem, 
                                      barrier.triangles = barrier.triangles, 
                                      prior.range = prior.range[k,], #21.0  0.5
                                      prior.sigma = prior.sigma, #1.0 0.1
                                      range.fraction = c(x.sensi[j,1],x.sensi[j,2],x.sensi[j,3]))
    model[[j+85]][[k]] <- tbm
    formula[[j+85]][[k]] <- y ~ 0 + b0 + f(i, model =  tbm)
    
    res.dng[[j+85]][[k]] <- inla(formula[[j+85]][[k]],
                                 data = inla.stack.data(joint.stk),
                                 family = 'poisson', 
                                 control.predictor = list(A = inla.stack.A(joint.stk), link = 1), 
                                 E = inla.stack.data(joint.stk)$e, 
                                 inla.mode = "experimental")
  }
}

par(mfrow=c(2,4))
for (j in 1:length(model[[87]])) {
  local.plot.field.book(
    res.dng[[87]][[j]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F)#, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  
}

##mode
prior.range <- matrix(c(21,0.5, 21,0.9, 
                        41,0.5, 41,0.9,
                        10,0.5, 10,0.9,
                        6,0.5, 6,0.9), ncol = 2, byrow = TRUE)

prior.sigma <- c(1, 0.1)
x.sensi <- rbind(x3.total[c(1,5,12,15),], x3.0)

statist <- "mode"
for (k in 1:nrow(prior.range)) {
  fp <- file.path("plots/dng", paste0("sensitivity/", prior.range[k,1], "_", prior.range[k,2], "/", statist))
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  path = fp
  
  for (j in 1:nrow(x.sensi)) {
    fpp <- file.path(paste0(path, "/", "b1_", x.sensi[j,1], "b2_", x.sensi[j,2], "b3_", x.sensi[j,3], ".png"))
    
    png(fpp)
    local.plot.field.book(
      res.dng[[85+j]][[k]]$summary.random$i$mode,
      #  asp = 1,
      xlim = c(-960, -905),
      ylim = c(2870, 2919),
      axes = F, zlim = c(-1.6,6))
    plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
    #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
    #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
    plot(poly.water.book2, add = T, col = "white", lty = 0)
    plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
    points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
    
    dev.off()
  }
}

statist <- "mean"
for (k in 1:nrow(prior.range)) {
  fp <- file.path("plots/dng", paste0("sensitivity/", prior.range[k,1], "_", prior.range[k,2], "/", statist))
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  path = fp
  
  for (j in 1:nrow(x.sensi)) {
    fpp <- file.path(paste0(path, "/", "b1_", x.sensi[j,1], "b2_", x.sensi[j,2], "b3_", x.sensi[j,3], ".png"))
    
    png(fpp)
    local.plot.field.book(
      res.dng[[85+j]][[k]]$summary.random$i$mean,
      #  asp = 1,
      xlim = c(-960, -905),
      ylim = c(2870, 2919),
      axes = F, zlim = c(-1.6,6))
    plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
    #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
    #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
    plot(poly.water.book2, add = T, col = "white", lty = 0)
    plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
    points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
    
    dev.off()
  }
}

##tables
prior.table <- list()
for (j in 1:nrow(x.sensi)) {
  prior.table[[j]] <- list()
  for (k in 1:nrow(prior.range)) {
    prior.table[[j]][[k]] <- res.dng[[85+j]][[k]]$summary.hyperpar
  }
}

#sensitivity table for article/overleaf
prior.sensi <- matrix(c(21,0.5, 21,0.9, 
                        41,0.5, 41,0.9,
                        10,0.5, 10,0.9,
                        6,0.5, 6,0.9), ncol = 2, byrow = TRUE)
#prior.sigma c(1,0.1)

sel.j <- x.sensi[c(5,6,4,3),]
sensi.table <- list()
for (j in 1:nrow(sel.j)) {
  sensi.table[[j]] <- list()
  for (k in 1:nrow(prior.sensi)) {
    sensi.table[[j]][[k]]  <- as.data.frame(exp(res.dng[[85+j]][[k]]$summary.hyperpar))
  }
}

##########################################################################################################################
#mesh plots
spatr_deep.sea <- terra::rast(raster_deep.sea)
# SpatVector
spatvect_deep.sea <- terra::as.polygons(spatr_deep.sea > -Inf) # bathym_poly
#SpatialPolygonsDataFrame 
spdf.deep.sea <- as(spatvect_deep.sea, "Spatial") # bathym_sp_df
sp4deep.sea <- spTransform(spdf.deep.sea, crs_km)
poly.b1.wo.deepsea <- gDifference(poly.bar1.new, sp4deep.sea)
poly.b1.wo.deepsea_sf <- st_as_sf(poly.b1.wo.deepsea)
poly.deepsea_sf <- st_as_sf(sp4deep.sea)
st_crs(poly.deepsea_sf) <- crs_km

poly.bar1.new_sf <- st_as_sf(poly.bar1.new) 
poly.bar2.new_sf <- st_as_sf(poly.bar2.new)
poly.bar3_sf <- st_as_sf(poly.bar3)
poly.water.book2_sf <- st_as_sf(poly.water.book2)
sp4deep.sea_sf <- st_as_sf(sp4deep.sea)

st_crs(poly.bar1.new_sf) <- crs_km
st_crs(poly.bar2.new_sf) <- crs_km
st_crs(poly.bar3_sf) <- crs_km
st_crs(poly.water.book2_sf) <- crs_km
st_crs(sp4deep.sea_sf) <- crs_km

#fpp <- file.path("plots/dng/mesh/mesh.void/mesh.dng.col1.png")
png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.6) +
  xlab("long") + ylab("lat") +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='red',size=1,alpha=0.6) +
  theme_void()
dev.off()

fpp <- file.path("plots/dng/mesh/mesh.void/mesh.dng.col3.png")
#png(fpp)

ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='red',size=1,alpha=0.6)  +
  xlab("long") + ylab("lat") +
  theme_void()
dev.off()

fpp <- file.path("plots/dng/mesh/mesh.void/mesh.dng.col2.png")
#png(fpp)

ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='red',size=1,alpha=0.6)  +
  xlab("long") + ylab("lat") +
  theme_void()
dev.off()

fpp <- file.path("plots/dng/mesh/mesh.dng.col0.png")
#png(fpp)

ggplot() +
  inlabru::gg(mesh.dng) +
  #geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  #geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='red',size=1,alpha=0.6)  +
  xlab("long") + ylab("lat")
dev.off()

###########
fpp <- file.path("plots/dng/mesh/mesh.dng.zoom1col2.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  xlim(c(-978, -907)) +
  ylim(c(2875,2942)) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=1.5,alpha=1) +
  xlab("long") + ylab("lat") 
dev.off()

fpp <- file.path("plots/dng/mesh/mesh.dng.zoom1col2.2.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  xlim(c(-978, -907)) +
  ylim(c(2875,2942)) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=1.5,alpha=1) +
  xlab("long") + ylab("lat") 
dev.off()


fpp <- file.path("plots/dng/mesh/mesh.dng.zoom2col2.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=1.5,alpha=1) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  xlim(c(-967, -907)) +
  ylim(c(2875,2942)) +
  xlab("long") + ylab("lat") 
dev.off()

fpp <- file.path("plots/dng/mesh/mesh.dng.zoom2col2.2.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=1.5,alpha=1) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  xlim(c(-967, -907)) +
  ylim(c(2875,2942)) +
  xlab("long") + ylab("lat") 
dev.off()

fpp <- file.path("plots/dng/mesh/mesh.dng.zoom1col3.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=1.5,alpha=1) +
  xlim(c(-978, -907)) +
  ylim(c(2875,2942)) +
  
  xlab("long") + ylab("lat")
dev.off()


fpp <- file.path("plots/dng/mesh/mesh.dng.zoom1col4.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=1.5,alpha=1) +
  xlim(c(-978, -907)) +
  ylim(c(2875,2942)) +
  
  xlab("long") + ylab("lat")
dev.off()


fpp <- file.path("plots/dng/mesh/mesh.dng.zoom2col3.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=1.5,alpha=1) +
  xlim(c(-967, -907)) +
  ylim(c(2875,2942)) +
  
  xlab("long") + ylab("lat")
dev.off()


fpp <- file.path("plots/dng/mesh/mesh.dng.zoom2col4.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=1.5,alpha=1) +
  xlim(c(-967, -907)) +
  ylim(c(2875,2942)) +
  xlab("long") + ylab("lat")
dev.off()

#############################################################################################################
#corr plots for points in dugong mesh
#from Dugong_article.Rmd line 2000
#locations chosen
loc.corr <- matrix(c(-919, 2880,
                     -909, 2875,
                     -911, 2880,
                     -914.0942, 2882.206,
                     
                     -923.0821, 2888.942,
                     -924, 2892,
                     -922.7535, 2896.290,
                     
                     -930.4101, 2900.680,
                     -935.7784, 2899.498,
                     -932.2049, 2900.386), ncol = 2, byrow = TRUE)

colnames(loc.corr) <- c("x", "y")
loc.corr  <- as.matrix(loc.corr)
poly_loc.corr <- Polygon(loc.corr)

sp_loc.corr  <- SpatialPoints(poly_loc.corr@coords)
sf_loc.corr  <- st_as_sf(sp_loc.corr)
st_crs(sf_loc.corr) <- crs_km

sf_group1 <- st_as_sf(SpatialPoints(Polygon(loc.corr[1:4,])@coords))
st_crs(sf_group1) <- crs_km

sf_group2 <- st_as_sf(SpatialPoints(Polygon(loc.corr[5:7,])@coords))
st_crs(sf_group2) <- crs_km

sf_group3 <- st_as_sf(SpatialPoints(Polygon(loc.corr[8:10,])@coords))
st_crs(sf_group3) <- crs_km

#find the mesh node closer to the chosen locations.
#Dugon_article: id.node.dng[[19]], change it to res.dng[[49]]
id.node.dng[[49]] <- id.node.tbm(mesh.dng, location = loc.corr, npoint = nrow(loc.corr))

df_group1_ <- matrix(c(id.node.dng[[49]]$id.coord[[1]][1], id.node.dng[[49]]$id.coord[[1]][2]), 
                     ncol = 2)
for (g in 2:4) {
  df_group1_ <- rbind(df_group1_,
                      matrix(c(id.node.dng[[49]]$id.coord[[g]][1], id.node.dng[[49]]$id.coord[[g]][2]), 
                             ncol = 2))
}
df_group1 <- st_as_sf(SpatialPoints(df_group1_))
st_crs(df_group1) <- crs_km

### group2
df_group2_ <- matrix(c(id.node.dng[[49]]$id.coord[[5]][1], id.node.dng[[49]]$id.coord[[5]][2]), 
                     ncol = 2)
for (g in 6:7) {
  df_group2_ <- rbind(df_group2_,
                      matrix(c(id.node.dng[[49]]$id.coord[[g]][1], id.node.dng[[49]]$id.coord[[g]][2]), 
                             ncol = 2))
}
df_group2 <- st_as_sf(SpatialPoints(df_group2_))
st_crs(df_group2) <- crs_km
###group 3
df_group3_ <- matrix(c(id.node.dng[[49]]$id.coord[[8]][1], id.node.dng[[49]]$id.coord[[8]][2]), 
                     ncol = 2)
for (g in 9:10) {
  df_group3_ <- rbind(df_group3_,
                      matrix(c(id.node.dng[[49]]$id.coord[[g]][1], id.node.dng[[49]]$id.coord[[g]][2]), 
                             ncol = 2))
}
df_group3 <- st_as_sf(SpatialPoints(df_group3_))
st_crs(df_group3) <- crs_km


##############################################################################################################

fpp <- file.path("plots/dng/corr.points/mesh.col3.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-950, -900)) +
  ylim(c(2860,2920)) +
  geom_sf(data =df_group1, shape = 17, col='red1',size=2) +
  geom_sf(data =df_group2, shape = 15, col='red1',size=2) +
  geom_sf(data =df_group3, shape = 19, col='red1',size=2) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.col4.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-950, -900)) +
  ylim(c(2860,2920)) +
  geom_sf(data =df_group1, shape = 17, col='red1',size=2) +
  geom_sf(data =df_group2, shape = 15, col='red1',size=2) +
  geom_sf(data =df_group3, shape = 19, col='red1',size=2) +
  xlab("long") + ylab("lat")
dev.off()

#zoom 2
fpp <- file.path("plots/dng/corr.points/mesh.col3.zoom2.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-955, -905)) +
  ylim(c(2862,2919)) +
  geom_sf(data =df_group1, shape = 17, col='red1',size=2) +
  geom_sf(data =df_group2, shape = 15, col='red1',size=2) +
  geom_sf(data =df_group3, shape = 19, col='red1',size=2) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.col4.zoom2.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-955, -905)) +
  ylim(c(2862,2919)) +
  geom_sf(data =df_group1, shape = 17, col='red1',size=2) +
  geom_sf(data =df_group2, shape = 15, col='red1',size=2) +
  geom_sf(data =df_group3, shape = 19, col='red1',size=2) +
  xlab("long") + ylab("lat")
dev.off()


##############################################################################################################

x.sel <- rbind(x3.total[c(1,5,12,15),], x3.0)
#prior.range <- matrix(c(21,0.5, 21,0.9, 41,0.5, 41,0.9,  10,0.5, 10,0.9,  6,0.5, 6,0.9), ncol = 2, byrow = TRUE)
#prior.sigma <- c(1, 0.1)

# all the points with the same range fraction x[1,]

corr.tbm <- function(barrier.model,
                     range.fraction = range.fraction,
                     id.node = id.node,
                     prior.range = prior.range,
                     prior.sigma = prior.sigma) {
  
  Q <- inla.rgeneric.q(
    barrier.model, 
    "Q",
    theta = c(log(prior.sigma[1]),
              log(prior.range[1])))
  
  sd <- sqrt(diag(inla.qinv(Q)))
  Inode <- rep(0, dim(Q)[1])
  Inode[id.node] <- 1
  
  covar.column <- solve(Q, Inode)
  corr = drop(matrix(covar.column))/(sd*sd[id.node])
  #cov = drop(matrix(covar.column))
  return(corr)
}


#IM DOING THIS LOOP

for (j in 1:nrow(x.sel)) {
  corr.dng[[j+85]] <- list()
  
  for (k in 1:nrow(prior.range)) {
    corr.dng[[j+85]][[k]] <- list()
    
    for (l in 1:nrow(loc.corr)) {
      corr.dng[[j+85]][[k]][[l]] <- 
        corr.tbm(barrier.model = model[[j+85]][[k]],
                 range.fraction = x.sel[j,], 
                 id.node = id.node.dng[[49]]$id.node[[l]],
                 prior.range = prior.range[k,],
                 prior.sigma = prior.sigma)
    }
  }
}

#save(model, file = "plots/dng/lists/model.tbm.RData")
#save(formula, file = "plots/dng/lists/formula.tbm.RData")
#save(res.dng, file = "plots/dng/lists/res.dng.RData")
#save(res.dng.pred, file = "plots/dng/lists/res.dng.pred.RData")
#save(corr.dng, file = "plots/dng/lists/corr.dng.points.RData")
#save(mesh.dng, file = "plots/dng/lists/mesh.dng.RData")

##############################################################################################################
for (j in 1:nrow(x.sel)) {
  corr.dng[[j+85]] <- list()
  
  for (k in 1:nrow(prior.range)) {
    corr.dng[[j+85]][[k]] <- list()
    
    for (l in 1:nrow(loc.corr)) {
      corr.dng[[j+85]][[k]][[l]] <- 
        corr.tbm(barrier.model = model[[j+85]][[k]],
                 range.fraction = x.sel[j,], 
                 id.node = id.node.dng[[49]]$id.node[[l]],
                 prior.range = prior.range[k,],
                 prior.sigma = prior.sigma)
    }
  }
}


for (l in 1:nrow(loc.corr)) {
  fp <- file.path("plots/dng/points_dng", paste0("point", point[l]))
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  path = fp
  
  for (j in 1:nrow(x.sel)){
    fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
    
    png(fpp)
    local.plot.field.book(
      corr.dng[[j+85]][[1]][[l]],
      #  asp = 1,
      xlim = c(-950, -900),
      ylim = c(2860,2920),
      axes = F, zlim = c(0.2, 1))
    plot(poly.bar1.new, add = TRUE, col = alpha("#996633", 0.5), lty = 0)
    plot(poly.bar2.new, add = TRUE, col = alpha("blue", 0.1), lty = 0)
    plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.5), lty = 0)
    #plot(poly.bar3, add = T, col = alpha("white", 0.3), lty = 0)
    plot(poly.bar3, add = T, col = alpha("skyblue", 0.5), lty = 0)
    #plot(poly.bar3, add = T, col = alpha("yellow", 0.3), lty = 0)
    plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
    points(id.node.dng[[49]]$id.coord[[l]][1], id.node.dng[[49]]$id.coord[[l]][2], pch = 20)
    
    dev.off()
    
  }
}

#zoom2

for (l in 1:nrow(loc.corr)) {
  fp <- file.path("plots/dng/points_dng.zoom0", paste0("point", point[l]))
  if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
  path = fp
  
  for (j in 1:nrow(x.sel)){
    fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
    
    png(fpp)
    local.plot.field.book(
      corr.dng[[j+85]][[1]][[l]],
      #  asp = 1,
      xlim = c(-960, -900),
      ylim = c(2860,2920),
      axes = F, zlim = c(0.2, 1))
    plot(poly.bar1.new, add = TRUE, col = alpha("#996633", 0.5), lty = 0)
    plot(poly.bar2.new, add = TRUE, col = alpha("blue", 0.1), lty = 0)
    plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.5), lty = 0)
    #plot(poly.bar3, add = T, col = alpha("white", 0.3), lty = 0)
    plot(poly.bar3, add = T, col = alpha("skyblue", 0.5), lty = 0)
    #plot(poly.bar3, add = T, col = alpha("yellow", 0.3), lty = 0)
    plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
    points(id.node.dng[[49]]$id.coord[[l]][1], id.node.dng[[49]]$id.coord[[l]][2], pch = 20)
    
    dev.off()
    
  }
}

###############################################################################################################
#some corrected plots for overleaf
loc.corr <- matrix(c(-909, 2875,
                     -930.4101, 2900.680,
                     -935.7784, 2899.498), ncol = 2, byrow = TRUE)

colnames(loc.corr) <- c("x", "y")
loc.corr  <- as.matrix(loc.corr)
poly_loc.corr <- Polygon(loc.corr)

sp_loc.corr  <- SpatialPoints(poly_loc.corr@coords)
sf_loc.corr  <- st_as_sf(sp_loc.corr)
st_crs(sf_loc.corr) <- crs_km

#find the mesh node closer to the chosen locations.
#Dugon_article: id.node.dng[[19]], change it to res.dng[[49]]
id.node.dng[[1]] <- id.node.tbm(mesh.dng, location = loc.corr, npoint = nrow(loc.corr))

df_group.sel <- matrix(c(id.node.dng[[1]]$id.coord[[1]][1], id.node.dng[[1]]$id.coord[[1]][2]), 
                       ncol = 2)
for (g in 2:3) {
  df_group.sel <- rbind(df_group.sel,
                        matrix(c(id.node.dng[[1]]$id.coord[[g]][1], id.node.dng[[1]]$id.coord[[g]][2]), 
                               ncol = 2))
}
df_group.sel <- st_as_sf(SpatialPoints(df_group.sel))
st_crs(df_group.sel) <- crs_km


#zoom = to the field plots
fpp <- file.path("plots/dng/mesh/zoom3.col0.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  #geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  #geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  #geom_sf(data =df_group.sel, shape = 17, col='red1',size=3) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=2.5,alpha=0.6) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()


fpp <- file.path("plots/dng/mesh/zoom3.col1.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=2.5,alpha=0.6) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/mesh/zoom3.col2.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=2.5,alpha=0.6) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/mesh/zoom3.col3.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=2.5,alpha=0.6) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/mesh/zoom3.col4.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "#000333", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "#000333", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data=st_as_sf(cords_Dng.sp),
          col='black',size=2.5,alpha=0.6) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()

#zoom = to the field plots
fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col0.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  #geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  #geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=2,alpha=0.5) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col0.void.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  #geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  #geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=2,alpha=0.5) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat") +
  theme_void()
dev.off()


fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col1.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col1.void.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat") +
  theme_void()
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col2.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col2.void.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  #geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat") +
  theme_void()
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col3.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col3.void.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "azure4", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "blue", alpha = 0.1) +
  geom_sf(data = poly.bar3_sf, fill = "white", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat") +
  theme_void()
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col4.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "#000333", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "#000333", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat")
dev.off()

fpp <- file.path("plots/dng/corr.points/mesh.zoom3.col4.void.png")
#png(fpp)
ggplot() +
  inlabru::gg(mesh.dng) +
  geom_sf(data = poly.bar1.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar2.new_sf, fill = "#000333", alpha = 0.3) +
  #geom_sf(data = poly.bar2.new_sf, fill = "#996633", alpha = 0.3) +
  geom_sf(data = poly.bar3_sf, fill = "#000333", alpha = 0.3) +
  #geom_sf(data = poly.bar3_sf, fill = "skyblue", alpha = 0.2) +
  #geom_sf(data = poly.bar3_sf, fill = "yellow", alpha = 0.1) +
  geom_sf(data = sp4deep.sea_sf, fill = "#000333", alpha = 0.7) +
  geom_sf(data =df_group.sel, shape = 17, col='red1',size=5) +
  #geom_sf(data=st_as_sf(cords_Dng.sp),
  #       col='black',size=1.5,alpha=1) +
  xlim(c(-960, -905)) +
  ylim(c(2870,2919)) +
  xlab("long") + ylab("lat") +
  theme_void()
dev.off()

###############################################################################
#corrected field plots for overleaf
#colored b0

zoom <- 0
zlim1 <- -1
b.col <- "b0"
##mean
statist <- "mean"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1,6))
  #plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  #plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##0.025
statist <- "0.025"

#change the values lower than -8 to -8 so  they are colored blue for 0.01, [[49]]
new0.025quant <- list()
for (j in 1:nrow(x.sel)) {
  new0.025quant[[j]] <- 
    ifelse(res.dng[[85+j]][[1]]$summary.random$i$'0.025quant' < -12, -12, res.dng[[49]]$summary.random$i$'0.025quant')
  
}
#do the second 0.025 with the new lower. values = to 7.99 so they appear blue
zlim1 <- c(-12)
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new0.025quant[[j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  #plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  #plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.5
zlim1 <- c(-12)
statist <- "0.5"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$'0.5quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  #plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  #plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##mode
zlim1 <- c(-1)
statist <- "mode"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,5.5))
  #plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  #plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#sd and 0.975 correction
##sd
zlim1 <- -1
statist <- "sd"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.sd <- list()
for (j in 1:nrow(x.sel)) {
  new.sd[[85+j]] <- ifelse(res.dng[[85+j]][[1]]$summary.random$i$sd > 3.7, 3.7, res.dng[[85+j]][[1]]$summary.random$i$sd)
}


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.sd[[85+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  #plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  #plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.0.975 <- list()
for (j in 1:nrow(x.sel)) {
  new.0.975[[85+j]] <- ifelse(res.dng[[85+j]][[1]]$summary.random$i$'0.975quant' > 12, 12, res.dng[[85+j]][[1]]$summary.random$i$'0.975quant')
}


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.0.975[[85+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  #plot(poly.bar1.new, add = TRUE, col = alpha("#996633",1), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  #plot(sp4deep.sea, add = T, col = alpha("#000333", 1), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

################################################################################################
#colored b1

zoom <- 0
zlim1 <- -1
b.col <- "b1"
##mean
statist <- "mean"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##0.025
statist <- "0.025"

#change the values lower than -8 to -8 so  they are colored blue for 0.01, [[49]]
new0.025quant <- list()
for (j in 1:nrow(x.sel)) {
  new0.025quant[[j]] <- 
    ifelse(res.dng[[85+j]][[1]]$summary.random$i$'0.025quant' < -12, -12, res.dng[[49]]$summary.random$i$'0.025quant')
  
}
#do the second 0.025 with the new lower. values = to 7.99 so they appear blue
zlim1 <- c(-12)
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new0.025quant[[j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.5
zlim1 <- c(-12)
statist <- "0.5"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$'0.5quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##mode
zlim1 <- c(-1)
statist <- "mode"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,5.5))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#sd and 0.975 correction
##sd
zlim1 <- -1
statist <- "sd"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.sd <- list()
for (j in 1:nrow(x.sel)) {
  new.sd[[85+j]] <- ifelse(res.dng[[85+j]][[1]]$summary.random$i$sd > 3.7, 3.7, res.dng[[85+j]][[1]]$summary.random$i$sd)
}


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.sd[[85+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.0.975 <- list()
for (j in 1:nrow(x.sel)) {
  new.0.975[[85+j]] <- ifelse(res.dng[[85+j]][[1]]$summary.random$i$'0.975quant' > 12, 12, res.dng[[85+j]][[1]]$summary.random$i$'0.975quant')
}


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.0.975[[85+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure4", 0.7), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

################################################################################################
#colored b1 and b2, so name b2

zoom <- 0
zlim1 <- -1
b.col <- "b2"
##mean
statist <- "mean"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##0.025
statist <- "0.025"

#change the values lower than -8 to -8 so  they are colored blue for 0.01, [[49]]
new0.025quant <- list()
for (j in 1:nrow(x.sel)) {
  new0.025quant[[j]] <- 
    ifelse(res.dng[[85+j]][[1]]$summary.random$i$'0.025quant' < -12, -12, res.dng[[49]]$summary.random$i$'0.025quant')
  
}
#do the second 0.025 with the new lower. values = to 7.99 so they appear blue
zlim1 <- c(-12)
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new0.025quant[[j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.5
zlim1 <- c(-12)
statist <- "0.5"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$'0.5quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##mode
zlim1 <- c(-1)
statist <- "mode"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,5.5))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#sd and 0.975 correction
##sd
zlim1 <- -1
statist <- "sd"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.sd <- list()
for (j in 1:nrow(x.sel)) {
  new.sd[[85+j]] <- ifelse(res.dng[[85+j]][[1]]$summary.random$i$sd > 3.7, 3.7, res.dng[[85+j]][[1]]$summary.random$i$sd)
}


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.sd[[85+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.0.975 <- list()
for (j in 1:nrow(x.sel)) {
  new.0.975[[85+j]] <- ifelse(res.dng[[85+j]][[1]]$summary.random$i$'0.975quant' > 12, 12, res.dng[[85+j]][[1]]$summary.random$i$'0.975quant')
}


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.0.975[[85+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

################################################################################################
#colored b1, b2 and b3, so name b3

zoom <- 0
zlim1 <- -1
b.col <- "b3"
##mean
statist <- "mean"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$mean,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1,6))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##0.025
statist <- "0.025"

#change the values lower than -8 to -8 so  they are colored blue for 0.01, [[49]]
new0.025quant <- list()
for (j in 1:nrow(x.sel)) {
  new0.025quant[[j]] <- 
    ifelse(res.dng[[85+j]][[1]]$summary.random$i$'0.025quant' < -12, -12, res.dng[[49]]$summary.random$i$'0.025quant')
  
}
#do the second 0.025 with the new lower. values = to 7.99 so they appear blue
zlim1 <- c(-12)
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new0.025quant[[j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.5
zlim1 <- c(-12)
statist <- "0.5"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$'0.5quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}


##mode
zlim1 <- c(-1)
statist <- "mode"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    res.dng[[85+j]][[1]]$summary.random$i$mode,
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-1.6,5.5))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

#sd and 0.975 correction
##sd
zlim1 <- -1
statist <- "sd"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.sd <- list()
for (j in 1:nrow(x.sel)) {
  new.sd[[85+j]] <- ifelse(res.dng[[85+j]][[1]]$summary.random$i$sd > 3.7, 3.7, res.dng[[85+j]][[1]]$summary.random$i$sd)
}


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.sd[[85+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,3.7))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

##0.975
zlim1 <- c(-12)
statist <- "0.975"
fp <- file.path("plots/dng/sel", paste0("col_", b.col, "/", statist))
if (!dir.exists(fp)) dir.create(fp, recursive = TRUE)
path = fp

new.0.975 <- list()
for (j in 1:nrow(x.sel)) {
  new.0.975[[85+j]] <- ifelse(res.dng[[85+j]][[1]]$summary.random$i$'0.975quant' > 12, 12, res.dng[[85+j]][[1]]$summary.random$i$'0.975quant')
}


for (j in 1:nrow(x.sel)) {
  fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  png(fpp)
  local.plot.field.book(
    new.0.975[[85+j]],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(-12,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  dev.off()
}

################################################################################################

local.plot.field.book(
  res.dng[[90]][[1]]$summary.fitted.values[idx.pred, "mean"],
  #  zlim = c(-30.22, 30),
  main = "", asp = 1, col = book.color.c(100),
  axes = FALSE)

par(mfrow=c(2,2))
for (j in 1:(nrow(x.sel)-2)) {
  #fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  #png(fpp)
  local.plot.field.book(
    res.dng[[87+j]][[1]]$summary.fitted.values[idx.pred, "0.5quant"],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0.01, 0.2))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  #points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  #dev.off()
}

for (j in 1:(nrow(x.sel)-2)) {
  #fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  #png(fpp)
  local.plot.field.book(
    res.dng[[87+j]][[1]]$summary.fitted.values[idx.pred, "0.975quant"],
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0.00, 2.5))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  #points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  #dev.off()
}

for (j in 1:(nrow(x.sel)-2)) {
  #fpp <- file.path(paste0(path, "/", "b1_", x.sel[j,1], "b2_", x.sel[j,2], "b3_", x.sel[j,3], ".png"))
  
  #png(fpp)
  local.plot.field.book(
    res.dng[[87+j]][[1]]$summary.random$i$'0.5quant',
    #  asp = 1,
    xlim = c(-960, -905),
    ylim = c(2870, 2919),
    axes = F, zlim = c(0,12))
  plot(poly.bar1.new, add = TRUE, col = alpha("white",1), lty = 0)
  plot(poly.bar1.new, add = TRUE, col = alpha("#996633",0.5), lty = 0)
  #plot(poly.bar2.new, add = TRUE, col = alpha("azure3", 0.6), lty = 0)
  #plot(poly.bar3, add = T, col = alpha("white", 0.8), lty = 0)
  plot(poly.water.book2, add = T, col = "white", lty = 0)
  plot(sp4deep.sea, add = T, col = alpha("#000333", 0.7), lty = 0)
  #points(xy[,1], xy[,2], pch = 20, col = alpha("black", 0.6))
  
  #dev.off()
}



