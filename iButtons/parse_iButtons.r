# script to read the iButton data

#setwd("~/Dropbox/ARC Laureate/iButtons/lawnhillibuttons")
setwd("//wallace.uds.anu.edu.au/shared data/Research/EEG/Moritz Lab/Dan/TempLogging/Kimberley")
library(plyr)
library(stringr)
library(reshape2)

rm(list=ls())

toMonthNumber <- function(textDate) {
  return(as.numeric(format(as.Date(daily_stats$date), format="%m")))
}

#file_pattern <- "LawnHill site"
file_pattern <- "Kimberley site "

start_date <- as.Date("1/7/13", "%d/%m/%y")  # remove data before this date
end_date <- as.Date("5/7/14", "%d/%m/%y")  # remove data before this date

files <- list.files(path="./",pattern = file_pattern, recursive = FALSE,ignore.case = TRUE, include.dirs = FALSE)

combined <- data.frame("site"=character(),"timecode"=numeric(),"temp"=numeric(),"date"=numeric(), "time_in_hours"=numeric(), "first_date"=numeric(),"last_date"=numeric(), "serial_no"=character(),"internalName"=character(),stringsAsFactors=F)

for (file in files) {

  button_data <- read.csv(file,stringsAsFactors=F)
  
  # extract header data
  serialNo <- button_data[2,2]
  internalName <- button_data[3,2]
  measurementNum <- nrow(button_data) - 6
  
  button_data <- button_data[-(1:6),]  #remove header rows
  button_data <- button_data[- which(button_data[,1]=="-end-"),]
  names(button_data) <- c("timecode","temp")
  button_data$timecode <- as.numeric(button_data$timecode)
  button_data$temp <- as.numeric(button_data$temp)
  button_data <- na.omit(button_data)

  
  site <- str_replace(file,".csv","")
  site <- str_replace(site,"LawnHill ","")
  site <- str_replace_all(site," ","_")
  site_rep <- rep(site,nrow(button_data))
  button_data <- cbind(site=site_rep,button_data)
  
  # convert dates
  date_origin <- as.Date("1900-01-01")
  button_data$date <- as.Date(button_data$timecode,origin=date_origin)
  button_data$time_in_hours <- round((button_data$timecode %% 1) * 24,2)  
  
  # remove data before placement
  button_data <- button_data[button_data$date >= start_date,]
  # and after collection
  button_data <- button_data[button_data$date < end_date,]
  
  first <- as.Date(min(button_data$timecode),date_origin)
  last  <- as.Date(max(button_data$timecode),date_origin)
  button_data$first_date <- first
  button_data$last_date <- last
  button_data$serial_no <- serialNo
  button_data$internalName <- internalName
  
  combined <- rbind(combined,button_data)
}

# write a combined iButton data table to file
# write.csv(combined,"Lawn_Hill_iButtons.csv",row.names=F)

### some data summary and plots ###


# summarise data by time of day
median_temp_by_hour <- tapply(combined$temp,INDEX=as.factor(round(combined$time_in_hours,1)),FUN=median)
median_temp_by_hour <- data.frame(time=row.names(median_temp_by_hour),median_temp=median_temp_by_hour,row.names=NULL,stringsAsFactors=F)
mean_temp_by_hour <- tapply(combined$temp,INDEX=as.factor(round(combined$time_in_hours,1)),FUN=mean)
mean_temp_by_hour <- data.frame(time=row.names(mean_temp_by_hour),mean_temp=mean_temp_by_hour,row.names=NULL,stringsAsFactors=F)

# summarise data by site
site_stats <- ddply(
                    .data = combined, 
                    .variables = c("site"),
                    .fun=summarize,
                      site_median = median(temp),
                      site_mean = mean(temp),
                      perc_95   = quantile(temp,0.95),
                      perc_05    = quantile(temp,0.05),
                      temp_range_90   = perc_95 - perc_05
                  )

# get daily extremes
site_date_temp <- combined[,c("site","date","temp")]
site_date_temp$date <- as.factor(site_date_temp$date)

daily_stats <- ddply(
  .data = site_date_temp, 
  .variables = c("site","date"),
  .fun=summarize,
  max_temp = max(temp),
  min_temp = min(temp),
  hours_above_32 = 1.5 * length(which(temp>32)),
  hours_below_20 = 1.5 * length(which(temp<20)),
  .parallel = F
)

daily_stats$month <- toMonthNumber(daily_stats$date)

stats_from_daily <- ddply(
  .data = daily_stats, 
  .variables = c("site"),
  .fun=summarize,
  N = length(date),
  daily_max_98 = quantile(max_temp,0.98),
  daily_min_02 = quantile(min_temp,0.02),
  .parallel = F
)

monthly_stats <- ddply(
  .data = daily_stats, 
  .variables = c("site","month"),
  .fun=summarize,
  N = length(date),
  mean_daily_max = mean(max_temp),
  mean_daily_min = mean(min_temp),
  mean_hours_above_32 = mean(hours_above_32),
  mean_hours_below_20 = mean(hours_below_20),
  .parallel = F
)


daily_max_by_site <- dcast(daily_stats, date ~ site, value.var="max_temp")
daily_min_by_site <- dcast(daily_stats, date ~ site, value.var="min_temp")

rm(site_date_temp, file, file_pattern, files, date_origin, internalName, serialNo, site, button_data, site_rep)


############### PLOTS ###############

windows()
par(mfrow=c(2,2))

boxplot(temp ~ site,data=combined,cex.axis=0.8,notch=T,range=1,main="Lawn Hill temperature range by site",xlab="Site",ylab="Temp")

# plot temperatures by hour
plot(combined[,c("time_in_hours","temp")],pch=20,xlab="Time of day",ylab="Temp",main="Lawn Hill temperature range by time of day")
points(median_temp_by_hour$time,median_temp_by_hour$median_temp,col="red",pch=20,cex=2)
#points(mean_temp_by_hour$time,mean_temp_by_hour$mean_temp,col="blue",pch=20,cex=2)

# plot daily maxima by site
boxplot(max_temp ~ site,data=daily_stats,cex.axis=0.8,notch=T,range=0.01,xlab="Time of day",ylab="Daily max temp",main="Lawn Hill daily maximum temperature")
plot(stats_from_daily$site,stats_from_daily$daily_max_98,col="red",add=T,axes=F)

# plot daily minima by site
boxplot(min_temp ~ site,data=daily_stats,cex.axis=0.8,notch=T,range=0.01,xlab="Time of day",ylab="Daily min temp",main="Lawn Hill daily minimum temperature")
plot(stats_from_daily$site,stats_from_daily$daily_min_02,col="blue",add=T,axes=F)

# a new plot for temperatures above a threshold
windows()
monthly_stats_above <- monthly_stats[,c("site","month","mean_hours_above_32")]
monthly_matrix   <- acast(monthly_stats_above, site ~ month, value.var="mean_hours_above_32")

monthly_stats_below <- monthly_stats[,c("site","month","mean_hours_below_20")]
monthly_matrix      <- acast(monthly_stats_above, site ~ month, value.var="mean_hours_above_32")

sites <- c(12,9,3,4) # row numbers of  the sites in the barplot

barplot(height=monthly_matrix[sites,],beside=T,col=rainbow(length(sites)), 
        xlab="Month", cex.lab=1.5, ylab="Daily hours above 32 degrees")
ymax <- max(monthly_matrix[sites,])
legend(15,ymax*.98,fill=rainbow(length(sites)),legend=row.names(monthly_matrix)[sites],cex=1.5)

# a new plot for temperatures below a threshold
windows()
monthly_stats_below <- monthly_stats[,c("site","month","mean_hours_below_20")]
monthly_matrix   <- acast(monthly_stats_below, site ~ month, value.var="mean_hours_belwo_20")

monthly_stats_below <- monthly_stats[,c("site","month","mean_hours_below_20")]
monthly_matrix      <- acast(monthly_stats_below, site ~ month, value.var="mean_hours_below_20")

sites <- c(12,9,3) # row numbers of  the sites in the barplot

barplot(height=monthly_matrix[sites,],beside=T,col=rainbow(length(sites)), 
        xlab="Month", cex.lab=1.5, ylab="Daily hours below 20 degrees")
ymax <- max(monthly_matrix[sites,])
legend(2,ymax*.98,fill=rainbow(length(sites)),legend=row.names(monthly_matrix)[sites],cex=1.5)
