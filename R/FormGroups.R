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
	
	# initialize databases
	driver = dbDriver("SQLite")
	if (is.character(dbFile)) {
		dbConn = dbConnect(driver, dbFile)
		on.exit(dbDisconnect(dbConn))
	} else {
		dbConn = dbFile
		if (!inherits(dbConn,"SQLiteConnection")) 
			stop("'dbFile' must be a character string or SQLiteConnection.")
		if (!dbIsValid(dbConn))
			stop("The connection has expired.")
	}
	
	if (is.na(match("rank",
		dbListFields(dbConn,
			tblName))))
		stop("No rank column in table.")
	
	searchExpression <- paste("select rank from",
		tblName)
	rs <- dbSendQuery(dbConn, searchExpression)
	allranks <- dbFetch(rs, n=-1, row.names=FALSE)$rank
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
	rank <- sort(table(taxonomy))
	
	searchResult <- data.frame(rank=names(rank),
		count=as.integer(rank),
		counts=as.integer(-rank),
		origin="",
		identifier="",
		stringsAsFactors=FALSE)
	
	rank <- names(rank)
	lineages <- strsplit(as.character(rank),
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
	for (i in seq_along(o)) {		if (searchResult$identifier[o[i]] == "") {			lineage <- lineages[[o[i]]]
			for (j in rev(seq_along(lineage))) {
				w <- startsWith(rank,
					paste(lineage[seq_len(j)],
						collapse=";"))
				w <- which(w)
				w <- w[searchResult$counts[w] < 0]
				counts <- sum(abs(searchResult$counts[w]))				
				if (counts >= goalSize) {
					if (counts > maxGroupSize &&
						j < length(lineage)) {						j <- j + 1 # go down one rank						w <- startsWith(rank,
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
										if (counts < minGroupSize) { # mark for later inclusion						searchResult$counts[w] <- -counts					} else {
						searchResult$counts[w] <- counts
					}
										searchResult$origin[w] <- origin
					id <- .change(lineage[j])
					if (id %in% searchResult$identifier[-w])
						id <- paste(id, o[i], sep="_")
					searchResult$identifier[w] <- id
					break				} else if (j == 1) { # create singleton group
					searchResult$origin[o[i]] <- ""
					id <- .change(lineage[j])
					if (id %in% searchResult$identifier[-o[i]])
						id <- paste(id, o[i], sep="_")
					searchResult$identifier[o[i]] <- id
				}
			}		}
		if (verbose)
			setTxtProgressBar(pBar, i/length(o))	}
	
	w <- which(!duplicated(allranks))
	allranks <- allranks[w]
	taxonomy <- taxonomy[w]
	m <- match(taxonomy, rank)
	searchResult <- searchResult[m,]
	searchResult$rank <- allranks
	searchResult$counts <- NULL
	
	if (is.character(add2tbl) || add2tbl) {		dbWriteTable(dbConn, "taxa", searchResult)		
		if (verbose)
			cat("\nUpdating column: \"identifier\"...")		searchExpression <- paste("update ",
			tblName,
			" set identifier = (select identifier from taxa where ",
			ifelse(is.character(add2tbl), add2tbl, tblName),
			".rank = taxa.rank)",
			sep="")		rs <- dbSendStatement(dbConn, searchExpression)		dbClearResult(rs)
		
		if (is.na(match("origin",
			dbListFields(dbConn,
				ifelse(is.character(add2tbl), add2tbl, tblName))))) {			searchExpression <- paste("alter table ",
				ifelse(is.character(add2tbl), add2tbl, tblName),
				" add column origin",
				sep="")			rs <- dbSendStatement(dbConn, searchExpression)
			dbClearResult(rs)		}
		
		if (verbose)
			cat("\nUpdating column: \"origin\"...")		searchExpression <- paste("update ",
			ifelse(is.character(add2tbl), add2tbl, tblName),
			" set origin = (select origin from taxa where ",
			ifelse(is.character(add2tbl), add2tbl, tblName),
			".rank = taxa.rank)",
			sep="")		rs <- dbSendStatement(dbConn, searchExpression)
		dbClearResult(rs)				searchExpression <- "drop table taxa"		rs <- dbSendStatement(dbConn, searchExpression)
		dbClearResult(rs)	}
	
	if (verbose) {
		cat("\nFormed",
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
