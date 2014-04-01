# combine Marxan outputs
library(raster)

combo.shp <- shapefile("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/marxan_result_combined3.shp")
sp_25.shp <- shapefile("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/sp_step/run_many_step/Marxan_result_sp_25.shp")
ph_25.shp <- shapefile("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/ph_step_batch/run__1_many/Marxan_result_100.shp")

ph_25 <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/ph_step_batch/run__1_many/marxan_result_new_sporder.csv")
sp_25 <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/sp_step/run_many_spf_16/marxan_result.csv")
ph_10 <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/ph_10_batch/run__1_many/marxan_result.csv")
sp_10 <- read.csv("C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/sp_10/run_many_10pc/marxan_result.csv")

combo.shp@data[order(combo.shp@data$marxan_pu),13] <- ph_25$freq
combo.shp@data[order(combo.shp@data$marxan_pu),14] <- sp_25$freq
combo.shp@data[order(combo.shp@data$marxan_pu),17] <- ph_10$freq
combo.shp@data[order(combo.shp@data$marxan_pu),18] <- sp_10$freq

combo.shp@data$sum_phsp25 <- combo.shp@data$freq_ph25 + combo.shp@data$freq_sp25
combo.shp@data$diff_25 <- combo.shp@data$freq_ph25 - combo.shp@data$freq_sp25
combo.shp@data$sum_phsp10 <- combo.shp@data$freq_ph10 + combo.shp@data$freq_sp10
combo.shp@data$diff_10 <- combo.shp@data$freq_ph10 - combo.shp@data$freq_sp10

shapefile(combo.shp,"C:/Users/u3579238/Work/Software/Marxan/Marxan_runs/global_mammals/marxan_result_combined4.shp", overwrite=T)
