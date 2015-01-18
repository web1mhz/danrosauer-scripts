rm(list=ls())

library(RMySQL)
library(gdata)  # for reading excel files
source("~/Work/Software/danrosauer-scripts/Database/DB_utilities.r")

user   <- 'danr'
cat("\nPassword for ",user,": ",sep="")
passwd <- scan(what=character(), n=1, quiet = T)
con <- moritz_db_login(user,passwd)

rm(passwd)

dbInfo   <- dbGetInfo(con)
dbTables <- dbListTables(con)

# load in a table and compare to the existing records
#new.filename   <- "~/Dropbox/ARC Laureate/sample info/Database/ForUploading/Catalano_data_upload_Eremiascincus_14i15.csv"
new.filename   <- "~/Dropbox/ARC Laureate/sample info/Database/ForUploading/Carlia_20141128_CN.csv"
db.table       <- "specimen"
matching_field <- "specimen_local_id"

delete_flag_field <- "delete"

#new.data       <- read.csv(new.filename, stringsAsFactors=F)
new.data       <- read.xls(new.filename, sheet = "Carlia_20141128_CN", header=T, stringsAsFactors=F)
new.fields     <- names(new.data)
new.ids        <- new.data[,matching_field]

# confirm that the dlete flag field is in new data
if (! delete_flag_field %in% new.fields) {
  err.text <- paste("\nDelete flag not in the new data:", missing_fields_new)
  stop(err.text)
}

# reduce the new data to the matching field and the delete flag field
new.data <- subset(new.data, select = c(matching_field, delete_flag_field))
delete.ids <- new.data[which(new.data[,delete_flag_field]==1), matching_field]

# start a transaction
#dbBegin(con)
records.updated <- 0

# loop through records which have an id
for (id in extant.new.ids) {

  # create sql to select a matching record from the database - with fields_to_update
  field_list_txt <- com_sep(fields_to_update)

  select.SQL <- paste("SELECT ",  field_list_txt, " from ", db.table, " WHERE ", matching_field, " = ", id, ";", sep="")
  result    <- dbSendQuery(con, select.SQL)
  db_row.df <- dbFetch(result)
  column_info <- dbColumnInfo(result)
  dbClearResult(result)

  rows <- which(new.data[,matching_field]==id)
  new_row.df <- new.data[rows,]

  # create sql for update statement
  field_value_txt <- value_list_sql(new_row.df, column_info, matching_field)

  update.SQL <- paste("UPDATE ", db.table, " SET ", field_value_txt,
                      " WHERE ",  matching_field, " = ", id, ";", sep="")
  result <- dbSendQuery(con, update.SQL)

  records.updated <- records.updated + 1
  cat("\n", id, "\t", records.updated, "\n\n")

  dbClearResult(result)
}

dbDisconnect(con)


