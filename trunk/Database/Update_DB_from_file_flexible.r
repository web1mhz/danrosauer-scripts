rm(list=ls())

library(RMySQL)

user   <- 'danr'
cat("\nPassword for ",user,": ",sep="")
passwd <- scan(what=character(), n=1, quiet = T)

con <- dbConnect(RMySQL::MySQL(), dbname="moritz_specimen", host = "moritzdb-moritzlab.azva.dotcloud.net", port=2600, user=user, password=passwd)
rm(passwd)

dbInfo   <- dbGetInfo(con)
dbTables <- dbListTables(con)

# load in a table and compare to the existing records
new.filename   <- "~/My Dropbox/ARC Laureate/sample info/Database/ForUploading/PMO_latlong_corrections.csv"
db.table       <- "specimen"
matching_field <- "specimen_local_id"

fields_to_update <- c("lineage_from_mtDNA")

new.data       <- read.csv(new.filename)
new.fields     <- names(new.data)

tbl.fields     <- dbListFields(con,db.table)

new.ids <- new.data[,matching_field]

for (id in new.ids) {
  
  # sql to select a matching record from the database
  select.SQL <- paste("SELECT * from ", db.table, " WHERE ", matching_field, "=", id, ";", sep="")
  result    <- dbSendQuery(con, select.SQL)
  db_row.df <- dbFetch(result)
  dbClearResult(result)
  
  rows <- which(new.data[,matching_field]==id)
  new_row.df <- new.data[rows,]
  
  cat("\n", matching_field, ":\t", id)
  cat("\nDatabase:", db_row.df$latitude, db_row.df$longitude, sep="\t")
  cat("\nNew data:", new_row.df$latitude, new_row.df$longitude, sep="\t")

  update.SQL <- paste("UPDATE ", db.table, " SET latitude = ", new_row.df$latitude, ", longitude = ", new_row.df$longitude,
                      " WHERE ",  matching_field, " = ", id, ";", sep="")
  result <- dbSendQuery(con, update.SQL)
  rows.affected <- dbGetRowsAffected(result)
  
  if (rows.affected == 1) {
    cat("\nUpdated\n\n")
  }
  
  dbClearResult(result)
}

dbDisconnect(con)


