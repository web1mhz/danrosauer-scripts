# This script inserts new records into the database.  Note - it assumes that records are new, and does not check
# Dan Rosauer - Feb 2015

rm(list=ls())

library(RMySQL)
library(stringr)
source("~/Work/Software/danrosauer-scripts/Database/DB_utilities.r")

#### PARAMETERS ####
user   <- 'danr'

new.filename   <- "~/Dropbox/ARC Laureate/sample info/Database/ForUploading/crenas_and_Pseudos.csv"
db.table       <- "specimen"

#fields_to_insert <- c("catalog_number", "ABTC_number", "lineage_from_mtDNA", "ND2_done", "notes")
fields_to_NOT_insert <- c("Delete", "specimen_local_id", "family", "last_change_date", "last_change_user") # update all fields, except as specified

# specify lookup fields
lookup <- data.frame(field="institution_full_name", lookup_table="institution", matching_field="institution_full_name", lookup_ID="institution_id", stringsAsFactors = F)
lookup[2,] <- c("genus", "genus", "genus", "genus_id")
lookup[3,] <- c("state_short", "state", "state_short", "state_ID")

# START FROM (to start not from the beginning)
#matching_field_start = 30855

####################

cat("\nPassword for ",user,": ",sep="")
passwd <- scan(what=character(), n=1, quiet = T)
con <- moritz_db_login(user,passwd)
rm(passwd)

dbInfo   <- dbGetInfo(con)
dbTables <- dbListTables(con)

tbl.fields     <- dbListFields(con,db.table)

new.data       <- read.csv(new.filename, stringsAsFactors=F)
new.fields     <- names(new.data)

# replace lookup fields with the corresponding ID
for (i in 1:nrow(lookup)) {
  new_column_pos <- which(tolower(new.fields) == lookup$field[i])
  #new.data[new_column_pos] <- str_trim(new.data[,new_column_pos])  # remove leading or trailing spaces
  id_column <- replace_lookup(input_vector=new.data[new_column_pos], connection=con, lookup_table=lookup$lookup_table[i], lookup_field=lookup$matching_field[i], id_field=lookup$lookup_ID[i])
  new.data[new_column_pos] <- id_column[,1]
  names(new.data)[new_column_pos] <- names(id_column[1])
  new.fields[new_column_pos] <- names(id_column[1])
  rm(id_column, new_column_pos)
}

if (! exists("fields_to_insert")) {fields_to_insert <- new.fields}

if (exists("fields_to_NOT_insert")) {
  fields_to_NOT_insert <- base::intersect(fields_to_insert, fields_to_NOT_insert)
  if(length(fields_to_NOT_insert) > 0) {
    fields_to_insert <- fields_to_insert[- which(fields_to_insert %in% fields_to_NOT_insert)]
  }
}

# confirm that fields_to_insert are in new data
missing_fields_new <- setdiff(fields_to_insert, new.fields)
if (length(missing_fields_new) > 0) {
  err.text <- paste("\nField to update not in the new data:", missing_fields_new)
  stop(err.text)
}

# confirm that fields_to_insert are in the database table
missing_fields_db <- setdiff(fields_to_insert, tbl.fields)
if (length(missing_fields_db) > 0) {
  err.text <- paste("\nField to update not in database table", db.table, ":", missing_fields_db)
  stop(err.text)
}

rm(missing_fields_new, missing_fields_db)

field_list_txt <- com_sep(fields_to_insert)

# reduce the new data to the matching fields
new.data <- subset(new.data,select = fields_to_insert)

# check for a start beyond the first record
# if (exists("matching_field_start")) {
#   starting_pos <- which(extant.new.ids==matching_field_start)
#   if (starting_pos > 0) {
#     extant.new.ids <- extant.new.ids[- (1:(starting_pos-1))]
#   }
# }

# create sql to select a  record from the database to get the correct column info
select.SQL <- paste("SELECT ",  field_list_txt, " from ", db.table, " LIMIT 1;", sep="")
result    <- dbSendQuery(con, select.SQL)
column_info <- dbColumnInfo(result)
dbClearResult(result)

# start a transaction
#dbBegin(con)
records.inserted <- 0

# loop through records which have an id
for (i in 1:nrow(new.data)) {

  # create sql for insert statement
  field_value_txt <- value_list_sql(new.data[i,], column_info, "", purpose="INSERT")

  insert.SQL <- paste("INSERT INTO ", db.table, "(", field_list_txt, ") VALUES (", field_value_txt, ");", sep="")
  result <- dbSendQuery(con, insert.SQL)

  records.inserted <- records.inserted + 1
  cat("\n", records.inserted, "\n\n")

  dbClearResult(result)

}

dbDisconnect(con)

