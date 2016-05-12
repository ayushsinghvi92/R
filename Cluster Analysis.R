##Preparing Data

#data for 2014 Well Reports
wellReports = read.csv('C:/Users/Ayush/Desktop/R/BA Project/Cluster Analysis/OK_Aggregated_UIC_Well_Reports.csv')
wellReports = subset(wellReports, wellReports$ReportYear =='2014')

#data for 2014 Wells
wells = read.csv('C:/Users/Ayush/Desktop/R/BA Project/Cluster Analysis/OK_ALL_UIC_Wells_v3.csv')
wells = subset(wells, wells$ReportYear=='2014')
#removing duplicates
wells = unique(wells[,1:3])

#writing CSV to remove NULL in Excel
write.csv(wells,file="C:/Users/Ayush/Desktop/R/BA Project/unique.csv")
wells = read.csv("C:/Users/Ayush/Desktop/R/BA Project/unique.csv")

#removing NAs
wells = na.omit(wells)
#Removing Trash data (some data with 6000 coordinates)
wells = subset(wells, wells$Latitude<100 & wells$Latitude>0 & wells$Longitude<100)


#data for 2014 earthquakes
earthquakes = read.csv('C:/Users/Ayush/Desktop/R/BA Project/Cluster Analysis/OK_EQ_v3.csv')
#Lubridate package required
require(lubridate)
earthquakes$EQ_Timestamp = mdy(earthquakes$EQ_Timestamp)
earthquakes = earthquakes[earthquakes$EQ_Timestamp > 
                            mdy(01012014) & earthquakes$EQ_Timestamp < 
                            mdy(01012015),]

#distances for earthquakes happening in 2014
require(lubridate)
wellToEarthquake = read.csv('C:/Users/Ayush/Desktop/R/BA Project/Cluster Analysis/OK_SWD_w2EQ_distance.csv')

wellToEarthquake$eq_timestamp = ymd_hms(wellToEarthquake$eq_timestamp)
wellToEarthquake = wellToEarthquake[wellToEarthquake$eq_timestamp > 
                                      ymd_hms(20140101000000) & 
                                      wellToEarthquake$eq_timestamp<
                                      ymd_hms(20150101000000),]
#removing abnormal rows
wellToEarthquake = subset (wellToEarthquake, wellToEarthquake$distance<15)

##Now, all wells are unique and from 2014, all earthquakes are from 2014, and well to earthquake distances
##are calculated for earthquakes that happened in 2014 without anomalies

#SP package required
require(sp)

#converting earthquake and wells into a spatial dataframes
coordinates(earthquakes)=~LONGITUDE+LATITUDE
coordinates(wells)=~Longitude+Latitude

#loading map of Oklahoma
#Requires package rgdal and raster
require(rgdal)
require(raster)

polygons <- shapefile("C:/Users/Ayush/Desktop/R/BA Project/Cluster Analysis/GeoData/oklahoma_administrative.shp")

##Creating a picture of the earthquakes and the wells

jpeg("C:/Users/Ayush/Desktop/R/BA Project/Cluster Analysis/Oklahoma.jpg",4000,2000,res=300)
plot(polygons)
title("Earthquakes and Wells",cex.main=3)
points(wells,cex=0.7,col="dark red")
points(earthquakes,col="blue",cex=0.5)
dev.off()

#setting the coordinate reference system for everything we're plotting

projection(earthquakes)=projection(polygons)
projection(wells)=projection(polygons)

##Creating a Distance Matrix

#assinging the standard mercator
#plot(wells)

wellsUTM <- spTransform(wells,CRS("+init=epsg:3395"))
earthquakesUTM <- spTransform(earthquakes,CRS("+init=epsg:3395"))


#Requires package rgeos for gdistance
require(rgeos)

distance.matrix <- matrix(0,nrow(earthquakes),7,dimnames=
    list(c(),c("EQ_ID","Lat","Long","Mag","EQ_Timestamp","Depth","DistWell")))
#loop to calculate distances (~47M calculations!)
for(i in 1:nrow(earthquakesUTM)){
  sub <- earthquakesUTM[i,]
  DistWell <- gDistance(sub,wellsUTM)
  distance.matrix[i,] <- matrix(c(sub$EQ_ID,sub@coords,sub$MAGNITUDE,sub$EQ_Timestamp,sub$DEPTH,gDistance(sub,wellsUTM)),ncol=7)
}

distanceMatrix = as.data.frame(distance.matrix)

#write.csv(distanceMatrix, "C:/Users/Ayush/Desktop/R/BA Project/Cluster Analysis/clusterData.csv")

distanceMatrix=read.csv("C:/Users/Ayush/Desktop/R/BA Project/Cluster Analysis/distance_matrix_2.0.csv")

#Finding optimum number of clusters
mydata <-  scale(distanceMatrix[,7])
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:30) wss[i] <- sum(kmeans(mydata,centers=i)$withinss)
plot(1:30, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups SSE")

##The graph indicates the relationship of SSE and number of clusters
#After 6 clusters, the SSE value doesnt decrease significantly so we choose to make 6 clusters

help("kmeans")
clust <- kmeans(mydata,6)


distanceMatrix$Clusters = clust$cluster

#Converting distance from meters to KM
distanceMatrix$DistWell=distanceMatrix$DistWell/1000

#Seeing what this actually means through 3D graph
require(scatterplot3d)

scatterplot3d(distanceMatrix$DistWell,xlab="Distance to Well",distanceMatrix$Depth,
              ylab="Depth",distanceMatrix$Mag,zlab="Magnitude", color = clust$cluster,
     pch=16,scale=1,grid=T,box=F, angle=100)

scatterplot3d(distanceMatrix$DistWell,xlab="Distance to Well",distanceMatrix$Depth,
              ylab="Depth",distanceMatrix$Mag,zlab="Magnitude", color = clust$cluster,
              pch=16,angle=320, scale=2,grid=T,box=F)

scatterplot3d(distanceMatrix$DistWell,xlab="Distance to Well",distanceMatrix$Depth,
              ylab="Depth",distanceMatrix$Mag,zlab="Magnitude", color = clust$cluster,
              pch=16,angle=220, scale=2,grid=T,box=F)

scatterplot3d(distanceMatrix$DistWell,xlab="Distance to Well",distanceMatrix$Depth,
              ylab="Depth",distanceMatrix$Mag,zlab="Magnitude", color = clust$cluster,
              pch=16,angle=120, scale=2,grid=T,box=F)

##Plotting Clusters on a Map
clustSP <- SpatialPointsDataFrame(coords=earthquakes@coords,data=data.frame(Clusters=clust$cluster))

jpeg("C:/Users/Ayush/Desktop/R/BA Project/Cluster Analysis/Eq.jpg",4000,2000,res=300)
plot(polygons)
title("Earthquakes",cex.main=3)
points(clustSP,col=clustSP$Clusters,cex=1, pch="*")
points(wells,cex=1,col="red",pch="*")

dev.off()