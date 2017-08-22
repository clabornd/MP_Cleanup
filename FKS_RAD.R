library(ggplot2)
library(ggmap)
library(dplyr)
library(lubridate)
library(scales)
library(geoR)
library(spacetime)
library(sp)
library(zoo)
library(readr)
library(gstat)
library(raster)
library(purrr)
library(gtools)

############### THIS FILE IS ONLY TO BUILD THE DATA, VARIOUS CHUNKS CAN BE  ##################
############### RUN TO GET AN APPROPRIATE DATAFRAME TO USE IN KRIGING FILES ##################
############### SOME PLOTTING CODE INCLUDED AS WELL  ---------------------  ##################

Rad <- read_csv("D:/Data/Rad.ExZone.csv")
#Rad.grouped <- read.csv("Rad.grouped.csv")

time.origin <- 15044

######CONSTRUCT DATA FRAME##############
#Begin with big data, aggregate as necessary:

#Remove future data, NA's.
Rad <- Rad[year(Rad$`Captured Time`) <= 2017,]
Rad <- Rad[is.na(Rad$`Captured Time`) == FALSE & is.na(Rad$Value)==FALSE,]
Rad <- Rad[Rad$Longitude < 141.03,]


#Set variables, define spatial resolution.
Rad.rounded <-  Rad %>% mutate(Lat.dec = round(Rad$Latitude, 2),
                      Long.dec = round(Rad$Longitude, 2),
                      year = as.integer(year(Rad$`Captured Time`))-min(as.integer(year(Rad$`Captured Time`)))+1,
                      month = as.integer(month(Rad$`Captured Time`))) %>%
                      dplyr::select(year, month, Lat.dec, Long.dec, dist, Value) %>%
                      mutate(time.int = 12*year+month)

Rad.rounded <- Rad.rounded %>% group_by(year, month, Lat.dec, Long.dec) %>% summarise_each(funs(mean), Value, time.int, year, Lat.dec, Long.dec, dist)

#Work with only saturated years
Rad.rounded <- Rad.rounded[Rad.rounded$year == 2 | Rad.rounded$year == 3 | Rad.rounded$year == 4 | Rad.rounded$year == 5,] #Get only certain years

#####SMALLER EXCLUSION ZONE#####

P2 <- Polygon(cbind(bbox_small$X2, bbox_small$X1))
Ps2 = SpatialPolygons(list(Polygons(list(P2), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(Ps2, axes = TRUE)

inner <- SpatialPointsDataFrame(cbind(Rad.rounded$Long.dec, Rad.rounded$Lat.dec), data = Rad.rounded[,-c(3,4)],
                        proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

inner <- inner[Ps2,]
inner <- inner[inner@coords[,2] > 37.50 | inner$dist < 18.642,]

inner.df <- as.data.frame(inner)

inner.df <- inner.df %>% rename(Long.dec = coords.x1, Lat.dec = coords.x2)
#various data storing

save(inner, file = "SPDF_small_601.RData")
write.csv(inner.df, "D:/Data/Rad_smallzone_601.csv")
write.csv(Rad.rounded, "D:/Data/Rad_ymon_531.csv")
write.csv(Rad.rounded, "D:/Data/Rad_rounded_coarse_529.csv")


################ USE NLS TO GET RESIDUALS ####################

estimate <- nls(data = Rad.clean, control = nls.control(maxiter = 1000),
                log(Value+1) ~ z + 1/dist^a+time.int*b+Long.dec*c+Lat.dec*d, 
                start = list(a = 2, b = 1, c = 1, d = 1, z =1))

#log scale
estimate <- nls(data = Rad.rounded, control = nls.control(maxiter = 1000),
                log(Value+1) ~ z + 1/dist^a+time.int*b+Long.dec*c+Lat.dec*d, 
                start = list(a = 2, b = 1, c = 1, d = 1, z = 1))

summary(estimate)
qplot(fitted(estimate), residuals(estimate))

#Add Residuals to dataframe
Rad.rounded["resids"] <- residuals(estimate)
Rad.rounded["fitted"] <- fitted(estimate)
Rad.rounded["time.ymon"] <- rep(as.Date("2011-01-01"), length(Rad.rounded$year))+
  years(Rad.rounded$year)+
  months(Rad.rounded$month-1)

#Average by time and location, add POSITX dates
Rad.clean <- Rad.rounded %>% 
              group_by(year,Lat.dec, Long.dec) %>% 
              summarise_each(funs(mean), Value, dist, resids, fitted)

Rad.clean["time"] <- rep(as.Date("2011-01-01"), length(Rad.clean$year))+years(Rad.clean$year)+months(Rad.clean$month)

Rad.clean <- Rad.clean[Rad.clean$Value < 1000,]

write.csv(Rad.clean, "D:/Data/Rad.full524_Resids.csv")


############# GENERIC PLOTTING CODEBLOCK #####################

longlims <- c(min(Rad.rounded$Long.dec), max(Rad.rounded$Long.dec))
latlims <- c(min(Rad.rounded$Lat.dec), max(Rad.rounded$Lat.dec))

map_image <- get_map(location = c(lon = mean(longlims), lat = mean(latlims)), zoom = 9, source = "google")

ggmap(map_image)+
  geom_tile(data = inner.df, aes(x = coords.x1, y = coords.x2, fill = Value), shape = 21)+
  geom_point(aes(x = 141.03, y = 37.42, size = 6), color = "red")

ggmap(map_image)+
  geom_point(data = inner_grid, aes(x = Var1, y = Var2, fill = estimates), shape = 21)+
  geom_point(aes(x = 141.03, y = 37.42, size = 6), color = "red")

###############PLOTTING FOR ALL DATA###################
#Get Latitude and Longitude Limits for Background Map
longlims <- c(min(inner.df$coords.x1), max(inner.df$coords.x1))
latlims <- c(min(inner.df$coords.x2), max(inner.df$coords.x2))

map_image <- get_map(location = c(lon = mean(longlims), lat = mean(latlims)), zoom = 10, source = "google")

ggmap(map_image)+
  geom_point(data = subset(inner.df, Value >= 629), aes(x = coords.x1, y = coords.x2), color = "red")+
  geom_point(data = subset(inner.df, Value < 629 & Value >= 50), aes(x = coords.x1, y = coords.x2, fill = Value), shape = 21)+
  scale_fill_gradient("CPM", low = "blue", high = "orange")+
  geom_point(data = subset(inner.df, Value < 50 & Value > 0), aes(x = coords.x1, y = coords.x2, color = Value))+
  scale_color_gradient("CPM")+
  geom_point(data = subset(inner.df, Value == 0), aes(x = coords.x1, y = coords.x2), color = "black")+
  scale_size_continuous(name = "")+
  geom_point(aes(x = 141.03, y = 37.42, size = 6), color = "red")


#####By Year#####
ggmap(map_image)+
  geom_point(data = subset(Rad.rounded, Value >= 629 & year == 4), aes(x = Long.dec, y = Lat.dec), color = "red")+
  geom_point(data = subset(Rad.rounded, Value < 629 & Value >= 50 & year == 4), aes(x = Long.dec, y = Lat.dec, fill = Value), shape = 21)+
  scale_fill_gradient("CPM", low = "blue", high = "orange")+
  geom_point(data = subset(Rad.rounded, Value < 50 & Value > 0 & year == 4), aes(x = Long.dec, y = Lat.dec, color = Value))+
  scale_color_gradient("CPM")+
  geom_point(data = subset(Rad.rounded, Value == 0 & year == 4), aes(x = Long.dec, y = Lat.dec), color = "black")+
  scale_size_continuous(name = "")+
  geom_point(aes(x = 141.03, y = 37.42, size = 6), color = "red")+
  theme(legend.position = "none")

ggmap(map_image)+
  geom_point(data = subset(Rad.grouped, Value >= 629 & year == 1), aes(x = Longitude, y = Latitude), color = "red")+
  geom_point(data = subset(Rad.grouped, Value < 629 & Value >= 50 & year == 1), aes(x = Longitude, y = Latitude, fill = Value), shape = 21)+
  scale_fill_gradient("CPM", low = "blue", high = "orange")+
  geom_point(data = subset(Rad.grouped, Value < 50 & Value > 0 & year == 1), aes(x = Longitude, y = Latitude, color = Value))+
  scale_color_gradient("CPM")+
  geom_point(data = subset(Rad.grouped, Value == 0 & year == 1), aes(x = Longitude, y = Latitude), color = "black")+
  scale_size_continuous(name = "")+
  geom_point(aes(x = 141.03, y = 37.42, size = 6), color = "red")

ggmap(map_image)+
  geom_point(data = subset(Rad.grouped, Value >= 629 & year == 2), aes(x = Longitude, y = Latitude), color = "red")+
  geom_point(data = subset(Rad.grouped, Value < 629 & Value >= 50 & year == 2), aes(x = Longitude, y = Latitude, fill = Value), shape = 21)+
  scale_fill_gradient("CPM", low = "blue", high = "orange")+
  geom_point(data = subset(Rad.grouped, Value < 50 & Value > 0 & year == 2), aes(x = Longitude, y = Latitude, color = Value))+
  scale_color_gradient("CPM")+
  geom_point(data = subset(Rad.grouped, Value == 0 & year == 2), aes(x = Longitude, y = Latitude), color = "black")+
  scale_size_continuous(name = "")+
  geom_point(aes(x = 141.03, y = 37.42, size = 6), color = "red")

ggmap(map_image)+
  geom_point(data = subset(Rad.grouped, Value >= 629 & year == 3), aes(x = Longitude, y = Latitude), color = "red")+
  geom_point(data = subset(Rad.grouped, Value < 629 & Value >= 50 & year == 3), aes(x = Longitude, y = Latitude, fill = Value), shape = 21)+
  scale_fill_gradient("CPM", low = "blue", high = "orange")+
  geom_point(data = subset(Rad.grouped, Value < 50 & Value > 0 & year == 3), aes(x = Longitude, y = Latitude, color = Value))+
  scale_color_gradient("CPM")+
  geom_point(data = subset(Rad.grouped, Value == 0 & year == 3), aes(x = Longitude, y = Latitude), color = "black")+
  scale_size_continuous(name = "")+
  geom_point(aes(x = 141.03, y = 37.42, size = 6), color = "red")

ggmap(map_image)+
  geom_point(data = subset(Rad.grouped, Value >= 629 & year == 4), aes(x = Longitude, y = Latitude), color = "red")+
  geom_point(data = subset(Rad.grouped, Value < 629 & Value >= 50 & year == 4), aes(x = Longitude, y = Latitude, fill = Value), shape = 21)+
  scale_fill_gradient("CPM", low = "blue", high = "orange")+
  geom_point(data = subset(Rad.grouped, Value < 50 & Value > 0 & year == 4), aes(x = Longitude, y = Latitude, color = Value))+
  scale_color_gradient("CPM")+
  geom_point(data = subset(Rad.grouped, Value == 0 & year == 4), aes(x = Longitude, y = Latitude), color = "black")+
  scale_size_continuous(name = "")+
  geom_point(aes(x = 141.03, y = 37.42, size = 6), color = "red")

ggmap(map_image)+
  geom_point(data = subset(Rad.full, Value >= 629 & year == 5), aes(x = Long.dec, y = Lat.dec), color = "red")+
  geom_point(data = subset(Rad.full, Value < 629 & Value >= 50 & year == 5), aes(x = Long.dec, y = Lat.dec, fill = Value), shape = 21)+
  scale_fill_gradient("CPM", low = "blue", high = "orange")+
  geom_point(data = subset(Rad.full, Value < 50 & Value > 0 & year == 5), aes(x = Long.dec, y = Lat.dec, color = Value))+
  scale_color_gradient("CPM")+
  geom_point(data = subset(Rad.full, Value == 0 & year == 5), aes(x = Long.dec, y = Lat.dec), color = "black")+
  scale_size_continuous(name = "")+
  geom_point(aes(x = 141.03, y = 37.42, size = 6), color = "red")



#### DUMP ####

#year = as.integer(year(Rad$`Captured Time`))-min(as.integer(year(Rad$`Captured Time`))

#Radiation %>% mutate(time = strsplit(as.character(Radiation$`Captured Time`)," ")[[1]][1]
#time = rep(as.Date("2011-01-01"), length(Radiation$year))+years(Radiation$year)+months(Radiation$month)    #Old Time Statements
#jdate = as.integer(julian(Radiation$`Captured Time`))-min(as.integer(julian(Radiation$`Captured Time`)))

#USE NONLINEAR LEAST SQUARES TO GET ESTIMATES FOR ALPHA IN REGRESSION EQUATION:#
##y~1/d^alpha+Vt+epsilon???? (something of this form)

#Create Distance Vector
#dist <- 0
#for(i in 1:dim(Radiation.clean)[1]){
#  dist[i]<-sqrt((Radiation.clean$Long.dec[i]-141)^2+(Radiation.grouped$Lat.dec[i]-37.4)^2)
#}

#dist <- 0
#for(i in 1:dim(Radiation.flat.time)[1]){
#  dist[i]<-sqrt((Radiation.flat.time$Longitude[i]-141)^2+(Radiation.flat.time$Latitude[i]-37.4)^2)
#}

#Radiation.clean["dist"] <- dist
#Radiation.flat.time["dist"] <- dist
#write.csv(Rad.clean, "Rad.full510.csv")
