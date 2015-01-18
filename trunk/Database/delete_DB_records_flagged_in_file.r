rm(list=ls())

library(RMySQL)
source("~/Work/Software/danrosauer-scripts/Database/DB_utilities.r")

user   <- 'danr'
cat("\nPassword for ",user,": ",sep="")
passwd <- scan(what=character(), n=1, quiet = T)
con <- moritz_db_login(user,passwd)

rm(passwd)

dbInfo   <- dbGetInfo(con)
dbTables <- dbListTables(con)

# load in a table and compare to the existing records
new.filename   <- "~/Dropbox/ARC Laureate/sample info/Database/ForUploading/Carlia_20141128_CN.csv"
db.table       <- "specimen"
matching_field <- "specimen_local_id"

delete_flag_field <- "Delete"

new.data       <- read.csv(new.filename, stringsAsFactors=F)
new.fields     <- names(new.data)

# confirm that the delete flag field is in new data
if (! delete_flag_field %in% new.fields) {
  err.text <- paste("\nDelete flag not in the new data:", new.fields)
  stop(err.text)
}

# reduce the new data to the matching field and the delete flag field
new.data <- subset(new.data, select = c(matching_field, delete_flag_field))
delete.ids <- new.data[which(new.data[,delete_flag_field]==1), matching_field]

# start a transaction
dbBegin(con)
records.deleted <- 0

cat("\nAbout to delete", length(delete.ids), "records")

# loop through records which have an id
for (id in delete.ids) {

  # create sql to delete flagged records

  delete.SQL <- paste("DELETE from ", db.table, " WHERE ", matching_field, " = ", id, ";", sep="")
  result    <- dbSendQuery(con, delete.SQL)
  dbClearResult(result)

  records.deleted <- records.deleted + 1
  cat("\n", id, "\t", records.deleted, "\n")

  dbClearResult(result)
}

dbCommit(con)
cat("\nDelete committed for", length(delete.ids), "records")

dbDisconnect(con)
