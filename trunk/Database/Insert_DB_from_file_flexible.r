# This script inserts new records into the database.  Note - it assumes that records are new, and does not check
# Dan Rosauer - Feb 2015

rm(list=ls())

library(RMySQL)
library(stringr)
source("~/Work/Software/danrosauer-scripts/Database/DB_utilities.r")

#### PARAMETERS ####
user   <- 'danr'

new.filename   <- "~/Dropbox/ARC Laureate/sample info/Database/ForUploading/Oliver_Nov14_Upload.xls"
db.table       <- "specimen"
matching_field <- "specimen_local_id"  # for inserts this is only used to avoid insert records which are already in the db

#fields_to_insert <- c("catalog_number", "ABTC_number", "lineage_from_mtDNA", "ND2_done", "notes")
fields_to_NOT_insert <- c("Delete", "delete", "cat_num", "specimen_local_id", "family", "last_change_date", "last_change_user") # update all fields, except as specified

#fields_to_check      <- c("field_number", "catalog_number", "ABTC_number", "ALA_record_ID")
#fields_to_check      <- c("field_number", "ABTC_number", "ALA_record_ID")
fields_to_check      <- c("field_number", "ABTC_number", "catalog_number")

# specify lookup fields
lookup <- data.frame(field="institution_full_name", lookup_table="institution", matching_field="institution_full_name", lookup_ID="institution_id", stringsAsFactors = F)
lookup[2,] <- c("genus", "genus", "genus", "genus_id")
lookup[3,] <- c("state_short", "state", "state_short", "state_ID")

# START FROM (to start not from the beginning)
#matching_field_start = 30855

####################

cat("\nPassword for ",user,": ",sep="")
passwd <- scan(what=character(), n=1, quiet = T)
#passwd <- 'Gehyra'
con <- moritz_db_login(user,passwd)
rm(passwd)

dbInfo   <- dbGetInfo(con)
dbTables <- dbListTables(con)
tbl.fields     <- dbListFields(con,db.table)

# check the new data type
suffix <- last(unlist(strsplit(new.filename,".", fixed=T)))

if (suffix=="xls") {
  library(gdata)
  new.data       <- read.xls(new.filename, sheet = 1, stringsAsFactors=F)
} else {
  new.data       <- read.csv(new.filename, stringsAsFactors=F)
}

new.fields     <- names(new.data)

# trim whitespace
for (i in 1:ncol(new.data)) {
  colname <- names(new.data)[i]
  nrows    <- nrow(new.data)
  if (class(new.data[1:nrows,colname]) == "character") {
    new.data[,i] <- str_trim(new.data[,i])
  }
}
rm(nrows,colname,i)

# replace lookup fields with the corresponding ID
for (i in 1:nrow(lookup)) {
  new_column_pos <- which(tolower(new.fields) == lookup$field[i])
  #new.data[new_column_pos] <- str_trim(new.data[,new_column_pos])  # remove leading or trailing spaces
  id_column <- replace_lookup(input_vector=new.data[new_column_pos], connection=con, lookup_table=lookup$lookup_table[i], lookup_field=lookup$matching_field[i], id_field=lookup$lookup_ID[i], verbose=T)
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
  err.text <- paste("\nField to insert not in the new data:", missing_fields_new)
  stop(err.text)
}

# confirm that fields_to_insert are in the database table
missing_fields_db <- setdiff(fields_to_insert, tbl.fields)
if (length(missing_fields_db) > 0) {
  err.text <- paste("\nField to insert not in database table", db.table, ":", missing_fields_db)
  stop(err.text)
}

rm(missing_fields_new, missing_fields_db)

field_list_txt <- paste(matching_field, ",", com_sep(fields_to_insert))

# reduce the data to records which don't already have a unique ID in the database
matching_records <- which(new.data[,matching_field] > 0)
if (length(matching_records) > 0) {
  new.data <- new.data[-matching_records,]
}

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

records_not_inserted <- new.data[0,] # a blank data frame which matches the format of new.data

# loop through records which have an id
for (i in 1:nrow(new.data)) {

  # check for a matching record already in the database
  # NOTE: only check for matches on non-blank fields

  # create sql to select a matching record from the database - with fields_to_update
  blank_fields <- which(new.data[i,fields_to_check] == "" | is.na(new.data[i,fields_to_check]))
  fields_to_check_this_time <- fields_to_check[-blank_fields]
  check_data <- as.data.frame(new.data[i,fields_to_check_this_time], stringsAsFactors = F)
  names(check_data) <- fields_to_check_this_time
  where_clause_text <- value_list_sql(check_data, column_info, "", purpose="WHERE_OR")

  select.SQL <- paste("SELECT ", field_list_txt, " from ", db.table, " WHERE ", where_clause_text, ";", sep="")
  result    <- dbSendQuery(con, select.SQL)
  db_row.df <- dbFetch(result)
  #column_info <- dbColumnInfo(result)
  dbClearResult(result)

  if (nrow(db_row.df) > 0) { # if there are matching records
    records_not_inserted <- rbind(records_not_inserted, new.data[i,])  #rbind will be slow if there are lots of records not inserted
    cat("\nRecord", i, "skipped due to a match to existing records\n\n")

  } else {

    # create sql for insert statement
    blank_fields <- which(new.data[i,] == "" | is.na(new.data[i,]))
    fields_to_insert_this_time <- fields_to_insert[-blank_fields]
    field_value_txt <- value_list_sql(new.data[i,fields_to_insert_this_time], column_info, "", purpose="INSERT")
    field_list_this_time_txt <- com_sep(fields_to_insert_this_time)
    insert.SQL <- paste("INSERT INTO ", db.table, "(", field_list_this_time_txt, ") VALUES (", field_value_txt, ");", sep="")
    result <- dbSendQuery(con, insert.SQL)

    records.inserted <- records.inserted + 1
    cat("\nRecord", i, records.inserted, "\n\n")

    dbClearResult(result)
  }
}

dbDisconnect(con)

if (nrow(records_not_inserted) > 0) {
  cat("\n", nrow(records_not_inserted), "records were not inserted due to a match with an existing record.  Please check the 'records_not_inserted' data frame which lists them.\n")
}

