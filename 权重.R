
library(tmap)
library(sf)
library(spdep)
library(ggplot2)
library(deldir)
library(sp)
library(purrr)
library(here)


library(tmaptools)
library(tidyverse)
library(car)
library(janitor)
library(tidymodels)
library(spatialreg)
library(spgwr)
library(RColorBrewer)



# 加载数据------
clev.points <- st_read(here::here("regression2.5.json"))%>%
  st_transform(.,27700)

a <- read.csv("regression2.5.csv")

vehicle.points = st_as_sf(a, coords = c("lon", "lat"), crs = 4326, agr = "constant")
class(vehicle.points)
plot(vehicle.points)


summary(clev.points)
plot(clev.points)

# 可视化点数据------
LondonBoroughs <- st_read(here::here("shenzhen",
                                     "Shenzhen.shp"))%>%
  st_transform(.,27700)
qtm(LondonBoroughs) 

tmap_mode("view")
tm_shape(LondonBoroughs) +
  tm_polygons(col = NA, alpha = 0.5) +
  tm_shape(vehicle.points) +
  tm_dots(col = "blue")



# 泰森多边形------------
vtess <- deldir(a$lon, a$lat)
plot(vtess, wlines="tess", wpoints="none",
     lty=1)

voronoipolygons = function(thiess) {
  w = tile.list(thiess)
  polys = vector(mode='list', length=length(w))
  for (i in seq(along=polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1,])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  SP = SpatialPolygons(polys)
  voronoi = SpatialPolygonsDataFrame(SP, data=data.frame(dummy = seq(length(SP)), row.names=sapply(slot(SP, 'polygons'), 
                                                                                                   function(x) slot(x, 'ID'))))
}

v <- voronoipolygons(vtess)
plot(v)

vtess.sf <- st_as_sf(v)
plot(vtess.sf$geometry)

# 泰森多边形邻接权重
st_queen <- function(a, b = a) st_relate(a, b, pattern = "F***T****")

queen.sgbp <- st_queen(vtess.sf)

as.nb.sgbp <- function(x, ...) {
  attrs <- attributes(x)
  x <- lapply(x, function(i) { if(length(i) == 0L) 0L else i } )
  attributes(x) <- attrs
  class(x) <- "nb"
  x
}

queen.nb <- as.nb.sgbp(queen.sgbp)
queen.weights <- nb2listw(queen.nb)
queen.weights

queen.nb.card <- card(queen.nb)
ggplot() +
  geom_histogram(aes(x=queen.nb.card)) +
  xlab("Number of Neighbors")

summary(queen.nb)


plot(queen.nb,coords, lwd=.2, col="blue", cex = .5)


# 莫兰-----------------------
# 全局莫兰
moran <- moran(a$price, queen.weights, length(queen.nb), Szero(queen.weights))
moran$I

## Moran scatter plot
moran.plot(a$price,queen.weights)

## traffic accident Moran's I---------------------------------------------------

Trafficaccidentnumbermoran <- a %>%
  dplyr::select(price)%>%
  pull()%>%
  moran.test(., queen.weights)
Trafficaccidentnumbermoran


## traffic accident Geary's C---------------------------------------------------
TrafficaccidentnumberGeary <- a %>%
  dplyr::select(price)%>%
  pull()%>%
  geary.test(., queen.weights)
TrafficaccidentnumberGeary




## local Moran's I-------------------------------------------------------
Trafficaccidentnumberlocalmoran <- a %>%
  dplyr::select(price)%>%
  pull()%>%
  localmoran(., queen.weights)%>%
  as_tibble()
Trafficaccidentnumberlocalmoran



# dependent variable distribution plot-----------------------------------------


ggplot(clev.points, aes(x=price^0.75)) + 
  geom_histogram(aes(y = ..density..),
                 binwidth = 1000) + 
  geom_density(colour="red", 
               size=1, 
               adjust=1)
# Tukey------------------------------------------------------------------------
symbox(~price, 
       Trafficaccidentnumber, 
       na.rm=T,
       powers=seq(-3,3,by=.25)) 



## Multiple Linear regression---------------------------------------------------
Regressiondata<- clev.points%>%
  clean_names()%>%
  dplyr::select(price, 
                lon,
                lat,
                floor_area,
                number_floors,
                room_number,
                number_bathroom,
                school_district,
                rate,
                distance_metro,
                distance_hospital,
                distance_sea)

# model1
model1 <- Regressiondata %>%
  lm(price^0.75~
       lon+
       lat+
       floor_area+
       number_floors+
       room_number +
       number_bathroom+
       school_district+
       rate+
       distance_metro+
       distance_hospital+
       distance_sea,
     data=.)
summary(model1)

model2 <- Regressiondata %>%
  lm(price^0.75~
       lon+
       lat+
       floor_area+
       number_floors+
       number_bathroom+
       school_district+
       rate+
       distance_hospital,
     data=.)
summary(model2)

AIC(model2)



# vif test
vif(model2)

#print some model diagnositcs.
par(mfrow=c(2,2))
plot(model2)

# Durbin-Watson test
DW <- durbinWatsonTest(model2)
tidy(DW)

#and for future use, write the residuals out
model_data2 <- model2 %>%
  augment(., Regressiondata)





# lagrange Multiplier Diagnostics
lm.LMtests(model2,queen.weights, zero.policy=NULL, test="all", spChk=NULL, naSubset=TRUE)

# also add them to the shapelayer
Trafficaccidentnumber <- clev.points %>%
  mutate(model1resids = residuals(model2))



## Spatial autocorrelation------------------------------------------------------

# residual moran's I
Nearest_neighbour <- Trafficaccidentnumber %>%
  st_drop_geometry()%>%
  dplyr::select(model1resids)%>%
  pull()%>%
  moran.test(., queen.weights)%>%
  tidy()

Nearest_neighbour


## spatially lagged model-------------------------------------------------------

# slm
slag_dv_model3_knn4 <- lagsarlm(price^0.75~
                                  lon+
                                  lat+
                                  floor_area+
                                  number_floors+
                                  number_bathroom+
                                  school_district+
                                  rate+
                                  distance_hospital,
                                data=Trafficaccidentnumber, 
                                nb2listw(queen.nb, 
                                         style="C"), 
                                method = "eigen")

summary(slag_dv_model3_knn4)


# write out the residuals

Trafficaccidentnumber <- Trafficaccidentnumber %>%
  mutate(slag_dv_model3_knn_resids = residuals(slag_dv_model3_knn4))

KNN4Moran <- Trafficaccidentnumber %>%
  st_drop_geometry()%>%
  dplyr::select(slag_dv_model3_knn_resids)%>%
  pull()%>%
  moran.test(., queen.weights)%>%
  tidy()
KNN4Moran


## spatial error model----------------------------------------------------------
sem_model1 <- errorsarlm(price^0.75~
                           lon+
                           lat+
                           floor_area+
                           number_floors+
                           number_bathroom+
                           school_district+
                           rate+
                           distance_hospital,
                         data=Trafficaccidentnumber,
                         nb2listw(queen.nb, style="C"), 
                         method = "eigen")

tidy(sem_model1)
summary(sem_model1)

Trafficaccidentnumber <- Trafficaccidentnumber %>%
  mutate(sem_model1 = residuals(sem_model1))

KNN4Moranerr <- Trafficaccidentnumber %>%
  st_drop_geometry()%>%
  dplyr::select(sem_model1)%>%
  pull()%>%
  moran.test(., queen.weights)%>%
  tidy()
KNN4Moranerr









## GWR--------------------------------------------------------------------------


st_crs(Trafficaccidentnumber) = 27700

TrafficaccidentnumberSP <- Trafficaccidentnumber %>%
  as(., "Spatial")



#calculate kernel bandwidth
GWRbandwidth <- gwr.sel(price^0.75~
                          lon+
                          lat+
                          floor_area+
                          number_floors+
                          number_bathroom+
                          school_district+
                          rate+
                          distance_hospital,
                        data = TrafficaccidentnumberSP, 
                        coords=cbind(a$lon,a$lat),
                        adapt=T)
GWRbandwidth




#run the gwr model
gwr.model = gwr(price^0.75~
                  lon+
                  lat+
                  floor_area+
                  number_floors+
                  number_bathroom+
                  school_district+
                  rate+
                  distance_hospital,
                data = TrafficaccidentnumberSP, 
                coords=cbind(a$lon,a$lat), 
                adapt=GWRbandwidth, 
                hatmatrix=TRUE, 
                se.fit=TRUE)
gwr.model









# print regression results to html
library(stargazer)
stargazer(model2, slag_dv_model3_knn4, sem_model1, type = "html",
          title="Title: Regression Results")









