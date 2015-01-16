rm(list=ls())

library(RMySQL)
source("~/Work/Software/danrosauer-scripts/Database/DB_utilities.r")

user   <- 'danr'
cat("\nPassword for ",user,": ",sep="")
passwd <- "Portl89d" #scan(what=character(), n=1, quiet = T)
con <- moritz_db_login(user,passwd)

rm(passwd)

dbInfo   <- dbGetInfo(con)
dbTables <- dbListTables(con)

# load in a table and compare to the existing records
new.filename   <- "~/Dropbox/ARC Laureate/sample info/Database/ForUploading/Catalano_data_upload_Eremiascincus_14i15.csv"
db.table       <- "specimen"
matching_field <- "specimen_local_id"

fields_to_update <- c("catalog_number", "ABTC_number", "lineage_from_mtDNA", "ND2_done", "notes")

tbl.fields     <- dbListFields(con,db.table)

new.data       <- read.csv(new.filename, stringsAsFactors=F)
new.fields     <- names(new.data)
new.ids        <- new.data[,matching_field]

# confirm that fields_to_update are in new data
missing_fields_new <- setdiff(fields_to_update, new.fields)
if (length(missing_fields_new) > 0) {
  err.text <- paste("\nField to update not in the new data:", missing_fields_new)
  stop(err.text)
}

# confirm that fields_to_update are in the database table
missing_fields_db <- setdiff(fields_to_update, tbl.fields)
if (length(missing_fields_db) > 0) {
  err.text <- paste("\nField to update not in database table", db.table, ":", missing_fields_db)
  stop(err.text)
}

rm(missing_fields_new, missing_fields_db)

# reduce the new data to the matching field and fields to update
new.data <- subset(new.data,select = c(matching_field, fields_to_update))
extant.new.ids <- new.ids[which(!is.na(new.ids))]

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


