Add2DB <- function(myData,
	dbFile,
	tblName="Seqs",
	clause="",
	verbose=TRUE) {
	
	# error checking
	if (!is.data.frame(myData))
		stop("myData must be a data frame object.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (startsWith(tblName, "_"))
		stop("Invalid tblName.")
	if (grepl("temp", tblName, ignore.case=TRUE))
		stop(tblName, " is a reserved tblName.")
	if (!is.character(clause))
		stop("clause must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (ncol(myData) == 0)
		stop("myData contains no columns.")
	if (any(grepl(".", names(myData), fixed=TRUE)))
		stop("Column names cannot contain periods.")
	
	# initialize database
	if (is.character(dbFile)) {
		if (!requireNamespace("RSQLite", quietly=TRUE))
			stop("Package 'RSQLite' must be installed.")
		dbConn <- dbConnect(dbDriver("SQLite"), dbFile)
		on.exit(dbDisconnect(dbConn))
	} else {
		dbConn <- dbFile
		if (!dbIsValid(dbConn))
			stop("The connection has expired.")
	}
	
	if (verbose)
		time.1 <- Sys.time()
	
	result  <- dbListTables(dbConn)
	w <- which(result == tblName)
	if (length(w) == 0) { # need to make the table
		dbWriteTable(dbConn,
			tblName,
			myData,
			row.names=TRUE,
			overwrite=TRUE,
			append=FALSE)
		if (verbose)
			cat("Created table '", tblName, "'.\n", sep="")
	} else {
		result <- dbListFields(dbConn, tblName)
		colIDs <- names(myData)
		
		if (is.na(match("row_names", colIDs)))
			myData$row_names <- as.double(row.names(myData))
		
		x <- dbGetQuery(dbConn,
			paste("select ",
				dbQuoteIdentifier(dbConn, tblName),
				".",
				dbQuoteIdentifier(dbConn, "row_names"),
				" from ",
				dbQuoteIdentifier(dbConn, tblName),
				sep=""))
		m <- match(myData$row_names, x$row_names)
		if (any(is.na(m)))
			stop("row.names of myData are missing from '", tblName, "'.")
		
		dbWriteTable(dbConn, "temp", myData, overwrite=TRUE)
		for (i in seq_along(colIDs)) {
			if (is.na(match(colIDs[i], result))) { # add missing column
				expression <- paste("alter table",
					dbQuoteIdentifier(dbConn, tblName),
					"add column",
					dbQuoteIdentifier(dbConn, colIDs[i]),
					dbDataType(dbConn, myData[1, i]))
				if (verbose)
					cat("Expression:\n",
						paste(strwrap(expression,
								width=getOption("width") - 1L),
							collapse="\n"),
						"\n\n",
						sep="")
				rs <- dbSendStatement(dbConn, expression)
				dbClearResult(rs)
			}
			
			expression <- paste("update ",
				dbQuoteIdentifier(dbConn, tblName),
				" set ",
				dbQuoteIdentifier(dbConn, colIDs[i]),
				" = (select ",
				dbQuoteIdentifier(dbConn, "temp"),
				".",
				dbQuoteIdentifier(dbConn, colIDs[i]),
				" from ",
				dbQuoteIdentifier(dbConn, "temp"),
				" where ",
				dbQuoteIdentifier(dbConn, "temp"),
				".",
				dbQuoteIdentifier(dbConn, "row_names"),
				" = ",
				dbQuoteIdentifier(dbConn, tblName),
				".",
				dbQuoteIdentifier(dbConn, "row_names"),
				") where exists (select ",
				dbQuoteIdentifier(dbConn, "temp"),
				".",
				dbQuoteIdentifier(dbConn, colIDs[i]),
				" from ",
				dbQuoteIdentifier(dbConn, "temp"),
				" where ",
				dbQuoteIdentifier(dbConn, "temp"),
				".",
				dbQuoteIdentifier(dbConn, "row_names"),
				" = ",
				dbQuoteIdentifier(dbConn, tblName),
				".",
				dbQuoteIdentifier(dbConn, "row_names"),
				")",
				sep="")
			if (verbose)
				cat("Expression:\n",
					paste(strwrap(expression,
							width=getOption("width") - 1L),
						collapse="\n"),
					"\n\n",
					sep="")
			rs <- dbSendStatement(dbConn, expression)
			dbClearResult(rs)
		}
		expression <- "drop table temp"
		rs <- dbSendStatement(dbConn, expression)
		dbClearResult(rs)
		
		if (verbose) {
			cat("Added to table ",
				tblName,
				':  "',
				paste(names(myData)[-match("row_names",
					names(myData))],
				collapse='", "'),
				'".\n\n',
				sep="")
		}
	}
	
	if (verbose) { # print the elapsed time to update table
		time.2 <- Sys.time()
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	invisible(TRUE)
}
