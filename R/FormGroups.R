FormGroups <- function(dbFile,
	tblName="Seqs",
	goalSize=50,
	minGroupSize=25,
	maxGroupSize=5000,
	includeNames=FALSE,
	add2tbl=FALSE,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.numeric(goalSize))
		stop("goalSize must be a numeric.")
	if (goalSize <= 0)
		stop("goalSize must be greater than zero.")
	if (!is.numeric(minGroupSize) || minGroupSize < 0)
		stop("minGroupSize must be a numeric.")
	if (minGroupSize < 0)
		stop("minGroupSize must be greater than or equal to zero.")
	if (!is.numeric(maxGroupSize))
		stop("maxGroupSize must be a numeric.")
	if (minGroupSize > goalSize)
		stop("goalSize must be at least minGroupSize.")
	if (maxGroupSize < goalSize)
		stop("maxGroupSize must be at least goalSize.")
	if (!is.logical(includeNames))
		stop("includeNames must be a logical.")
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
	
	searchExpression <- paste("select",
		dbQuoteIdentifier(dbConn, "organism"),
		"from",
		dbQuoteIdentifier(dbConn, tblName))
	rs <- dbSendQuery(dbConn, searchExpression)
	allranks <- dbFetch(rs, n=-1, row.names=FALSE)$organism
	dbClearResult(rs)
	
	if (includeNames) {
		taxonomy <- unlist(lapply(strsplit(allranks,
				"\n",
				fixed=TRUE),
			function (x) {
				paste(paste(x[-1],
						collapse=" "),
					x[1],
					sep=ifelse(grepl("; *$", x[length(x)]),
						"",
						";"))
			}))
	} else {
		taxonomy <- unlist(lapply(strsplit(allranks,
				"\n",
				fixed=TRUE),
			function (x) {
				paste(x[-1], collapse=" ")
			}))
	}
	taxonomy <- gsub(" *; *",
		";",
		taxonomy)
	organism <- sort(table(taxonomy))
	
	searchResult <- data.frame(organism=names(organism),
		count=as.integer(organism),
		counts=as.integer(-organism),
		origin="",
		identifier="",
		stringsAsFactors=FALSE)
	
	organism <- names(organism)
	lineages <- strsplit(as.character(organism),
		";",
		fixed=TRUE)
	
	if (verbose)
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	
	.change <- function(id) {
		id <- .Call("replaceChar", id, '"', "_", PACKAGE="DECIPHER")
		id <- .Call("replaceChar", id, "'", "_", PACKAGE="DECIPHER")
		id <- .Call("replaceChar", id, ".", "_", PACKAGE="DECIPHER")
		id <- gsub("^\\s+|\\s+$", "", id) # trim flanking white space
		return(id)
	}
	
	o <- order(searchResult$count,
		lengths(lineages),
		decreasing=TRUE)
	for (i in seq_along(o)) {
		if (searchResult$identifier[o[i]] == "") {
			lineage <- lineages[[o[i]]]
			for (j in rev(seq_along(lineage))) {
				w <- startsWith(organism,
					paste(lineage[seq_len(j)],
						collapse=";"))
				w <- which(w)
				w <- w[searchResult$counts[w] < 0]
				counts <- sum(abs(searchResult$counts[w]))
				
				if (counts >= goalSize) {
					if (counts > maxGroupSize &&
						j < length(lineage)) {
						j <- j + 1 # go down one organism
						w <- startsWith(organism,
							paste(lineage[seq_len(j)],
								collapse=";"))
						w <- which(w)
						w <- w[searchResult$counts[w] < 0]
						counts <- sum(abs(searchResult$counts[w]))
					}
					
					if (j > 1) {
						origin <- paste(lineage[seq_len(j - 1)], collapse=";")
					} else {
						origin <- ""
					}
					
					if (counts < minGroupSize) { # mark for later inclusion
						searchResult$counts[w] <- -counts
					} else {
						searchResult$counts[w] <- counts
					}
					
					searchResult$origin[w] <- origin
					id <- .change(lineage[j])
					if (id %in% searchResult$identifier[-w])
						id <- paste(id, o[i], sep="_")
					searchResult$identifier[w] <- id
					break
				} else if (j == 1) { # create singleton group
					searchResult$origin[o[i]] <- ""
					id <- .change(lineage[j])
					if (id %in% searchResult$identifier[-o[i]])
						id <- paste(id, o[i], sep="_")
					searchResult$identifier[o[i]] <- id
				}
			}
		}
		if (verbose)
			setTxtProgressBar(pBar, i/length(o))
	}
	
	w <- which(!duplicated(allranks))
	allranks <- allranks[w]
	taxonomy <- taxonomy[w]
	m <- match(taxonomy, organism)
	searchResult <- searchResult[m,]
	searchResult$organism <- allranks
	searchResult$counts <- NULL
	
	if (is.character(add2tbl) || add2tbl) {
		dbWriteTable(dbConn, "temp", searchResult, overwrite=TRUE)
		
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
		
		if (is.na(match("origin",
			dbListFields(dbConn,
				ifelse(is.character(add2tbl), add2tbl, tblName))))) {
			searchExpression <- paste("alter table",
			ifelse(is.character(add2tbl),
				dbQuoteIdentifier(dbConn, add2tbl),
				dbQuoteIdentifier(dbConn, tblName)),
				"add column",
				dbQuoteIdentifier(dbConn, "origin"),
				dbDataType(dbConn, ""))
			rs <- dbSendStatement(dbConn, searchExpression)
			dbClearResult(rs)
		}
		
		searchExpression <- paste("update ",
			ifelse(is.character(add2tbl),
				dbQuoteIdentifier(dbConn, add2tbl),
				dbQuoteIdentifier(dbConn, tblName)),
			" set ",
			dbQuoteIdentifier(dbConn, "origin"),
			" = (select ",
			dbQuoteIdentifier(dbConn, "origin"),
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
		cat("\n\nFormed",
			length(unique(searchResult$identifier)),
			"distinct groups.")
		if (is.character(add2tbl) || add2tbl)
			cat("\nAdded to table ",
				ifelse(is.character(add2tbl), add2tbl, tblName),
				": \"identifier\", \"origin\".\n",
				sep="")
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(searchResult)
}
