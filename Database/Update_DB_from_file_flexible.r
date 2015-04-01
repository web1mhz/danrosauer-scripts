rm(list=ls())

library(RMySQL)
library(stringr)
source("~/Work/Software/danrosauer-scripts/Database/DB_utilities.r")

#### PARAMETERS ####
user   <- 'danr'

new.filename   <- "~/Work/AMT/Data/Cryptoblepharus/dan_endemism_modeling/CryptoDB_20141104_final duplicates_MB.csv"
db.table       <- "specimen"
matching_field <- "specimen_local_id"

#fields_to_update <- c("catalog_number", "ABTC_number", "lineage_from_mtDNA", "ND2_done", "notes")
fields_to_NOT_update <- c("delete", "cat_num", "family", "last_change_date", "last_change_user") # update all fields, except as specified

# specify lookup fields
lookup <- data.frame(field="institution_full_name", lookup_table="institution", matching_field="institution_full_name", lookup_ID="institution_id", stringsAsFactors = F)
lookup[2,] <- c("genus", "genus", "genus", "genus_id")
lookup[3,] <- c("state_short", "state", "state_short", "state_ID")

# START FROM (to start not from the beginning)
matching_field_start = 51169


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
  id_column <- replace_lookup(input_vector=new.data[new_column_pos], connection=con, lookup_table=lookup$lookup_table[i], lookup_field=lookup$matching_field[i], id_field=lookup$lookup_ID[i], verbose = T)
  new.data[new_column_pos] <- id_column[,1]
  names(new.data)[new_column_pos] <- names(id_column[1])
  new.fields[new_column_pos] <- names(id_column[1])
  rm(id_column, new_column_pos)
}

if (! exists("fields_to_update")) {
  fields_to_update <- new.fields
  fields_to_update <- fields_to_update[- which(fields_to_update==matching_field)]
}

if (exists("fields_to_NOT_update")) {
  fields_to_update <- fields_to_update[- which(fields_to_update %in% fields_to_NOT_update)]
}

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

field_list_txt <- com_sep(fields_to_update)

# remove rows where Delete=1
if ("delete" %in% new.fields) {
  new.data <- new.data[- which(new.data$delete==1),]
}

# reduce the new data to the matching field and fields to update
new.data <- subset(new.data,select = c(matching_field, fields_to_update))
extant.new.ids <- new.data[which(!is.na(new.data[,matching_field])), matching_field]

# check for a start beyond the first record
if (exists("matching_field_start")) {
  starting_pos <- which(extant.new.ids==matching_field_start)
  if (starting_pos > 0) {
    extant.new.ids <- extant.new.ids[- (1:(starting_pos-1))]
  }
}

# start a transaction
#dbBegin(con)
records.updated <- 0

# loop through records which have an id
for (id in extant.new.ids) {

  # create sql to select a matching record from the database - with fields_to_update
  select.SQL <- paste("SELECT ",  field_list_txt, " from ", db.table, " WHERE ", matching_field, " = ", id, ";", sep="")
  result    <- dbSendQuery(con, select.SQL)
  db_row.df <- dbFetch(result)
  column_info <- dbColumnInfo(result)
  dbClearResult(result)

  rows <- which(new.data[,matching_field]==id)
  new_row.df <- new.data[rows,]

  changed <- are_there_changes(new_row.df, db_row.df)

  ############################
  # add code here (or as a separate function) to identify records (and fields)
  # which have changed, and only update them
  ############################

  if (changed[[1]]) {
    # create sql for update statement
    field_value_txt <- value_list_sql(new_row.df, column_info, matching_field)

    update.SQL <- paste("UPDATE ", db.table, " SET ", field_value_txt,
                        " WHERE ",  matching_field, " = ", id, ";", sep="")
    result <- dbSendQuery(con, update.SQL)

    records.updated <- records.updated + 1
    cat("\n", id, "\t", records.updated, "\n\n")

    dbClearResult(result)
  } else {
      cat("\n", id, "\tNot changed\n")
  }
}

dbDisconnect(con)

