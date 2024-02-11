IdentifyByRank <- function(dbFile,
	tblName="Seqs",
	level=0,
	add2tbl=FALSE,
	verbose=TRUE) {
	
	# error checking:
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.numeric(level) || floor(level) != level)
		stop("level must be an integer.")
	
	if (verbose)
		time.1 <- Sys.time()
	
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
	
	if (is.na(match("organism",
		dbListFields(dbConn,
			tblName))))
		stop("No 'organism' column in ", tblName, ".")
	
	searchExpression <- paste("select distinct",
		dbQuoteIdentifier(dbConn, "organism"),
		"from",
		dbQuoteIdentifier(dbConn, tblName))
	rs <- dbSendQuery(dbConn, searchExpression)
	x <- dbFetch(rs, n=-1, row.names=FALSE)
	dbClearResult(rs)
	
	.change <- function(id) {
		id <- .Call("replaceChar", id, '"', "", PACKAGE="DECIPHER")
		id <- .Call("replaceChar", id, "'", "", PACKAGE="DECIPHER")
		id <- gsub("^\\s+|\\s+$", "", id) # trim flanking white space
		id <- gsub("\\.+$", "", id)
		return(id)
	}
	
	z <- x
	if (level == 0) {
		x <- strsplit(x$organism, "\n", fixed=TRUE)
		z$origin <- unlist(lapply(x,
			function (x) {
				x <- paste(x[-1], collapse=" ")
			}))
		z$identifier <- .change(unlist(lapply(x, `[`, 1L)))
	} else {
		x$organism <- unlist(lapply(strsplit(x$organism,
				"\n",
				fixed=TRUE),
			function (x) {
				x <- paste(x[-1], collapse=" ")
			}))
		for (j in seq_along(x$organism)) {
			a <- strsplit(x$organism[j], ";")[[1]]
			l <- length(a)
			if (level < 0) {
				temp_level <- l + level + 1L
			} else {
				temp_level <- level
			}
			
			if (temp_level > l) {
				id <- as.character(a[l])
			} else if (temp_level < 1) {
				id <- as.character(a[1])
			} else {
				id <- as.character(a[temp_level])
			}
			z$origin[j] <- unlist(strsplit(as.character(x$organism[j]),
				id,
				fixed=TRUE))[1]
			
			z$identifier[j] <- .change(id)
		}
	}
	
	if (is.character(add2tbl) || add2tbl) {
		dbWriteTable(dbConn, "temp", z, overwrite=TRUE)
		
		searchExpression <- paste("update ",
			ifelse(is.character(add2tbl),
				dbQuoteIdentifier(dbConn, add2tbl),
				dbQuoteIdentifier(dbConn, tblName)),
			" set ",
			dbQuoteIdentifier(dbConn, "identifier"),
			" = (select ",
			dbQuoteIdentifier(dbConn, "identifier"),
			" from ",
			dbQuoteIdentifier(dbConn, "temp"),
			" where ",
			ifelse(is.character(add2tbl),
				dbQuoteIdentifier(dbConn, add2tbl),
				dbQuoteIdentifier(dbConn, tblName)),
			".",
			dbQuoteIdentifier(dbConn, "organism"),
			" = ",
			dbQuoteIdentifier(dbConn, "temp"),
			".",
			dbQuoteIdentifier(dbConn, "organism"),
			")",
			sep="")
		rs <- dbSendStatement(dbConn, searchExpression)
		dbClearResult(rs)
		
		searchExpression <- paste("drop table",
			dbQuoteIdentifier(dbConn, "temp"))
		rs <- dbSendStatement(dbConn, searchExpression)
		dbClearResult(rs)
	}
	
	if (verbose) {
		cat("\nFormed",
			length(unique(z$identifier)),
			"distinct groups.")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to table ",
				ifelse(is.character(add2tbl), add2tbl, tblName),
				": \"identifier\".",
				sep="")
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	names(z)[1] <- "organism"
	return(z)
}
