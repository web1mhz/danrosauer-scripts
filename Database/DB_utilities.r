
library(RMySQL)
library(stringr)
library(dplyr)

source("~/Work/Software/danrosauer-scripts/Database/DB_details.r")

moritz_db_login <- function(user, passwd) {

  # check for an existing connection with the same name
  if (exists("con")) {
    if (class(con) == "MySQLConnection") {
      dbDisconnect(con)
    }
    rm(con)
  }

  db_details <- moritz_db_details()
  con <- dbConnect(RMySQL::MySQL(), dbname = db_details$dbname, host = db_details$host, port = db_details$port, user=user, password=passwd)

  # check that the connection was successful
  success=F
  if (exists("con")) {
    if (class(con) == "MySQLConnection") { success=T }
  }

  if (success) {
    cat("\nConnected to moritz_specimen")
    return(con)
  } else {
    cat("\nConnection to moritz_specimen failed")
    return(0)
  }
}

com_sep <- function(vec, separator = ", ") {
  # this function takes each item in a vector and returns a single charater string in
  # which each item is separated by the separator, eg a comma
  result <- ""
  n <- length(vec)
  for (i in 1:n) {
    item <- vec[i]
    if (i < n) {
      result <- paste(result, item, separator, sep="")
    } else {
      result <- paste(result, item, sep="")
    }
  }

  result <- str_trim(result)
  return(result)
}

value_list_sql <- function(data, column_info, matching_field = NULL, purpose = "UPDATE") {
  # takes a data frame with a single row, for the data to be updated or queried in the database
  # the column names must match the corresponding database columns
  # column_info is generated by dbColumnInfo, and is used to identify the data types, to
  # get quoting right

  # returns a text string in the format required for
  # an UPDATE statement, eg:
    # name = "gouldii", svl = 7, is_observation = 0
  # or an INSERT statement, eg:
    # "gouldii", 7, 0
  # or a WHERE clause, eg:
    # name = "gouldii OR svl = 7 or is_observation = 0

  col_count <- ncol(data)

  result <- ""

  if (purpose == "UPDATE") {
    for (i in 1:col_count) {
      col_name <- names(data)[i]

      # ignore the matching field if specified
      if (col_name == matching_field) {next}

      col_class <- as.character(column_info[which(column_info$name==col_name), "Sclass"])

      if (as.character(data[i]) == "NA") {
        result <- paste(result, names(data)[i], "= NULL", sep="")
      } else if (col_class=="character") {
        result <- paste(result, names(data)[i], ' = "', data[i], '"', sep="")
      } else { # assumes other column classes are unquoted
        result <- paste(result, names(data)[i], " = ", data[i], sep="")
      }

      if (i < col_count) {result <- paste(result, ", ", sep="") }
    }

  } else if (purpose == "INSERT") {
    for (i in 1:col_count) {
      col_name <- names(data)[i]
      col_class <- as.character(column_info[which(column_info$name==col_name), "Sclass"])

      if (is.na(data[1,i])) {
        result <- paste(result, "NULL", sep="")
      } else if (col_class=="character") {
        result <- paste(result, '"', data[1,i], '"', sep="")
      } else { # assumes other column classes are unquoted
        result <- paste(result, data[1,i], sep="")
      }

      if (i < col_count) {result <- paste(result, ", ", sep="") }
    }

  } else if (purpose == "WHERE_OR") {
    for (i in 1:col_count) {
      col_name <- names(data)[i]
      col_class <- as.character(column_info[which(column_info$name==col_name), "Sclass"])

      if (as.character(data[1,i]) == "NA") {
        result <- paste(result, names(data)[i], "= NULL", sep="")
      } else if (col_class=="character") {
        result <- paste(result, names(data)[i], ' = "', data[1,i], '"', sep="")
      } else { # assumes other column classes are unquoted
        result <- paste(result, names(data)[i], " = ", data[1,i], sep="")
      }

      if (i < col_count) {result <- paste(result, " OR ", sep="") }
    }

  }

  return(result)
}

replace_lookup <- function(input_vector, connection, lookup_table, lookup_field, id_field, verbose=F) {
  # for a vector of input values, returns the value, and its unique identifier
  # this is used to match a real value to the ID stored in a relational structure

  input <- data.frame(input_vector, stringsAsFactors=F)  # make the input vector a single column data frame
  names(input) <- lookup_field

  select.SQL <- paste("SELECT ", id_field, ", ", lookup_field, " FROM ", lookup_table, sep="")
  lookup    <- dbGetQuery(con, select.SQL)

  output <- left_join(x=input, y=lookup)

  if (verbose) {
    # print unmatched values to screen
    unmatched <- unique(output[is.na(output[,2]),1])

    if (length(unmatched) > 0) {
      cat("\nThe following values did not have a match in the lookup table", lookup_table, "\n")
      print(unmatched)
      cat("\n")
    } else {
      cat("\nAll values matched to the lookup table", lookup_table, "\n")
    }
  }

  return(output[,2:1])
}

are_there_changes <- function(new_row, db_row) {
  # returns true if the fields found in both the new_now and the db_row (matched by name) are different
  # if differences are found, returns FALSE, and the name of the first non-matching field
  # new_row and db_row should be a single row data frame


  n <- ncol(new_row)
  new_names <- names(new_row)
  db_names  <- names(db_row)

  blanks <- which(new_row[1,]=="")
  new_row[1,blanks] <- NA
  blanks <- which(db_row[1,]=="")
  db_row[1,blanks] <- NA

  same <- TRUE
  i    <- 0

  while (i <= n & same==TRUE) {
    i <- i+1
    db_colnum <- which(db_names==new_names[i])

    if (!length(db_colnum)>0) {next}  # is there a matching database column (in the selection)
    if (is.na(new_row[1,i]) & is.na(db_row[1,db_colnum])) {next} # values are equal because both are NA

    if (is.na(new_row[1,i]) | is.na(db_row[1,db_colnum])) { # values are equal because just one is NA
      # this preceding line shouldn't be reached if both are NA
      same <- FALSE
      next
    }

    if (new_row[1,i] == db_row[1,db_colnum]) {next} # values are different (but not NA)
    same <- FALSE
  }
return(list(!same, new_names[i]))

}
