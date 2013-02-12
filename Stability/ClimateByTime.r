# Dan Rosauer - February 2013

point.time.series = function (point_x,point_y,grid_name,directory) {
  
  library(SDMTools)
 
  ################################################################################
  #list the projections, cycle thorugh them and project the models onto them
  times <- data.frame(year=list.files(directory)) #list the projections
  times$value <- rep(NA,nrow(times))

  model_count = 0
  
  #cycle through the projections
  for (time in times$year) {
 
    model_count = model_count+1
  
    point <- data.frame(x=point_x,y=point_y)
    cat("Reading data for year", time,"\n")
    grid_path <- paste(directory,time,grid_name,sep="/")
    envgrid.asc <- read.asc.gz(grid_path)
    times$value[model_count] <- extract.data(point,envgrid.asc)

  }
  
  return(times)
  
}


model.time.series = function (point_x,point_y,grid_name,directory) {
  
  library(SDMTools)
  
  ################################################################################
  #list the projections, cycle thorugh them and project the models onto them
  times <- data.frame(year=list.files(directory,pattern='\\.asc')) #list the projections
  times$year = gsub('\\.asc','',times$year)
  
  times$value <- rep(NA,nrow(times))
  
  model_count = 0
  
  #cycle through the projections
  for (time in times$year) {
    
    model_count = model_count+1
    
    point <- data.frame(x=point_x,y=point_y)
    cat("Reading data for year", time,"\n")
    grid_path <- paste(directory,"/",time,".asc",sep="")
  
    modelgrid.asc <- read.asc(grid_path)
    times$value[model_count] <- extract.data(point,modelgrid.asc)
    
  }
  
  return(times)
  
}
