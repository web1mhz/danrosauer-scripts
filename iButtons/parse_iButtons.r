# script to read the iButton data

setwd("~/Dropbox/ARC Laureate/iButtons/lawnhillibuttons")
# library(plyr)
library(stringr)

rm(list=ls())

file_pattern <- "LawnHill site"
files <- list.files(path="./",pattern = file_pattern, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

combined <- data.frame("site"=character(),"timecode"=numeric(),"temp"=numeric(),"date"=numeric(), "time_in_hours"=numeric(), "first_date"=numeric(),"last_date"=numeric, "serial_no"=character(),"internalName"=character(),stringsAsFactors=F)

for (file in files) {

  button_data <- read.csv(file,stringsAsFactors=F)
  
  # extract header data
  serialNo <- button_data[2,2]
  internalName <- button_data[3,2]
  measurementNum <- nrow(button_data) - 6
  
  button_data <- button_data[-(1:6),]
  names(button_data) <- c("timecode","temp")
  button_data$timecode <- as.numeric(button_data$timecode)
  button_data$temp <- as.numeric(button_data$temp)
  button_data <- na.omit(button_data)

  
  site <- str_replace(file,".csv","")
  site_rep <- rep(site,nrow(button_data))
  button_data <- cbind(site=site_rep,button_data)
  
  # convert dates
  date_origin <- as.Date("1900-01-01")
  button_data$date <- as.Date(button_data$timecode,origin=date_origin)
  button_data$time_in_hours <- round((button_data$timecode %% 1) * 24,2)  
  
  # remove data before placement
  button_data <- button_data[button_data$date >= as.Date("30/4/13", "%d/%m/%y"),]
  # and after collection
  button_data <- button_data[button_data$date < as.Date("16/4/14", "%d/%m/%y"),]
  
  first <- as.Date(min(button_data$timecode),date_origin)
  last  <- as.Date(max(button_data$timecode),date_origin)
  button_data$first_date <- first
  button_data$last_date <- last
  button_data$serial_no <- serialNo
  button_data$internalName <- internalName
  
  combined <- rbind(combined,button_data)
}

# write a combined iButton data table to file
write.csv(combined,"Lawn_Hill_iButtons.csv",row.names=F)


### some data summary and plots ###

# summarise data by site
temp_median_by_site <- tapply(combined$temp,INDEX=as.factor(combined$site),FUN=median)
temp_mean_by_site <- tapply(combined$temp,INDEX=as.factor(combined$site),FUN=mean)
temp_95_by_site     <- tapply(combined$temp,INDEX=as.factor(combined$site),FUN=quantile,0.95)
temp_05_by_site     <- tapply(combined$temp,INDEX=as.factor(combined$site),FUN=quantile,0.05)
temp_range_90       <- temp_95_by_site - temp_05_by_site


par(mfrow=c(2,2))
boxplot(temp ~ site,data=combined,cex.axis=0.8,notch=T,range=1)

# summarise data by time of day
median_temp_by_hour <- tapply(button_data$temp,INDEX=as.factor(button_data$time_in_hours),FUN=median)
median_temp_by_hour <- data.frame(time=row.names(median_temp_by_hour),median_temp=median_temp_by_hour,row.names=NULL,stringsAsFactors=F)
mean_temp_by_hour <- tapply(button_data$temp,INDEX=as.factor(button_data$time_in_hours),FUN=mean)
mean_temp_by_hour <- data.frame(time=row.names(mean_temp_by_hour),mean_temp=mean_temp_by_hour,row.names=NULL,stringsAsFactors=F)

windows()
plot(button_data[,c("time_in_hours","temp")],pch=20)
points(median_temp_by_hour$time,median_temp_by_hour$median_temp,col="red",pch=20,cex=2)
points(mean_temp_by_hour$time,mean_temp_by_hour$mean_temp,col="blue",pch=20,cex=2)

