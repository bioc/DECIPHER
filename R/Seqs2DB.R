Seqs2DB <- function(seqs,
	type,
	dbFile,
	identifier,
	tblName="Seqs",
	chunkSize=1e7,
	replaceTbl=FALSE,
	fields=c(accession="ACCESSION", organism="ORGANISM"),
	processors=1,
	verbose=TRUE,
	...) {
	
	# initialize variables
	time.1 <- Sys.time()
	
	# error checking
	if (length(seqs)==0)
		stop("seqs is zero length.")
	SEQTYPES <- c("FASTA", "FASTQ", "GenBank", "XStringSet", "DNAStringSet", "RNAStringSet", "AAStringSet", "BStringSet", "QualityScaledXStringSet", "QualityScaledDNAStringSet", "QualityScaledRNAStringSet", "QualityScaledAAStringSet", "QualityScaledBStringSet")
	type <- pmatch(type, SEQTYPES)
	if (is.na(type))
		stop("Invalid seqs type.")
	if (type==-1)
		stop("Ambiguous seqs type.")
	if (length(type) > 1)
		stop("type must be length 1.")
	if (!is.character(identifier))
		stop("Identifier must be a character.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (startsWith(tblName, "_"))
		stop("Invalid tblName.")
	if (grepl("temp", tblName, ignore.case=TRUE))
		stop(tblName, " is a reserved tblName.")
	if (!is.numeric(chunkSize))
		stop("chunkSize must be a numeric.")
	if (floor(chunkSize) != chunkSize)
		stop("chunkSize must be a whole number.")
	if (chunkSize <= 0)
		stop("chunkSize must be greater than zero.")
	if (!is.logical(replaceTbl))
		stop("replaceTbl must be a logical.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors) != processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- .detectCores()
	} else {
		processors <- as.integer(processors)
	}
	if (!is(seqs, "XStringSet") &&
		!is(seqs, "QualityScaledXStringSet") &&
		type > 3)
		stop("seqs must be an XStringSet or QualityScaledXStringSet.")
	if (type==3 && length(fields) > 0) {
		if (!is.character(fields))
			stop("fields must be a character vector.")
		if (any(is.na(fields)))
			stop("fields cannot contain NA.")
		if (is.null(names(fields)) || any(is.na(names(fields))))
			stop("fields must be a named character vector.")
		if (any(names(fields) == ""))
			stop("The names of fields cannot be empty.")
		if (any(fields == ""))
			stop("fields cannot be empty.")
		if (any(duplicated(names(fields))))
			stop("The names of fields must be unique.")
		if (any(duplicated(fields)))
			stop("fields must be unique.")
		if ("DEFINITION" %in% fields)
			stop("fields cannot contain 'DEFINITION'.")
		if (any(names(fields) %in% c("description", "identifier", "row_names")))
			stop("The names of fields contains a reserved name.")
		if (any(nchar(fields) > 12))
			stop("fields may be at most 12 characters wide.")
		fields <- c(description="DEFINITION", fields)
	} else {
		fields <- c(description="DEFINITION")
	}
	
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
	result <- dbListTables(dbConn)
	w <- which(result == tblName)
	if (length(w) == 1L && length(which(result==paste("_", tblName, sep=""))) != 1)
		stop("Table is corrupted")
	f <- character()
	fts1 <- c("BIGINT PRIMARY KEY", rep(dbDataType(dbConn, list(raw())), 2))
	names(fts1) <- c("row_names", "sequence", "quality")
	if (type != 3) {
		fts2 <- c("BIGINT PRIMARY KEY", rep(dbDataType(dbConn, ""), 2))
		names(fts2) <- c("row_names", "identifier", "description")
	} else {
		fts2 <- c("BIGINT PRIMARY KEY",
			rep(dbDataType(dbConn, ""), length(fields) + 1L))
		names(fts2) <- c("row_names", "identifier", names(fields))
	}
	if (length(w) == 1L && !replaceTbl) {
		searchExpression <- paste("select max(",
			dbQuoteIdentifier(dbConn, "row_names"),
			") from ",
			dbQuoteIdentifier(dbConn, tblName),
			sep="")
		numSeq <- as.double(dbGetQuery(dbConn, searchExpression)[[1]])
		
		searchExpression <- paste("select max(",
			dbQuoteIdentifier(dbConn, "row_names"),
			") from ",
			dbQuoteIdentifier(dbConn, paste("_", tblName, sep="")),
			sep="")
		if (as.double(dbGetQuery(dbConn, searchExpression)[[1]]) != numSeq)
			stop("Table is corrupted.")
		
		# make sure all the necessary fields exist
		f <- dbListFields(dbConn, paste("_", tblName, sep=""))
		for (i in seq_along(fts1)) {
			if (is.na(match(names(fts1)[i], f))) {
				# first add the column if it does not already exist
				expression1 <- paste("alter table",
					dbQuoteIdentifier(dbConn, paste("_", tblName, sep="")),
					"add column",
					dbQuoteIdentifier(dbConn, names(fts1)[i]),
					fts1[i])
				if (verbose)
					cat("Expression:  ", expression1, "\n", sep="")
				rs <- dbSendStatement(dbConn, expression1)
				dbClearResult(rs)
			}
		}
		
		# make sure all the necessary fields exist
		f <- dbListFields(dbConn, tblName)
		for (i in seq_along(fts2)) {
			if (is.na(match(names(fts2)[i], f))) {
				# first add the column if it does not already exist
				expression1 <- paste("alter table",
					dbQuoteIdentifier(dbConn, tblName),
					"add column",
					dbQuoteIdentifier(dbConn, names(fts2)[i]),
					fts2[i])
				if (verbose)
					cat("Expression:  ", expression1, "\n", sep="")
				rs <- dbSendStatement(dbConn, expression1)
				dbClearResult(rs)
			}
		}
	} else {
		numSeq <- 0
		replaceTbl <- TRUE # necessary for field types
	}
	
	if (verbose) {
		it <- 0L
		priorWidth <- 0L
		.cat <- function(text, ...) {
			width <- nchar(text)
			if (width < priorWidth) {
				text <- paste(text,
					paste(rep(" ", priorWidth - width),
						collapse=""),
					sep="")
			}
			priorWidth <<- width
			
			cat(text, sep="")
			flush.console()
		}
	}
	
	if (type==1) { # FASTA
		if (is.character(seqs)) {
			if (substr(seqs, 1, 7)=="http://" ||
				substr(seqs, 1, 8)=="https://" ||
				substr(seqs, 1, 7)=="ftps://" ||
				substr(seqs, 1, 6)=="ftp://") {
				con <- gzcon(url(seqs, "rb"))
			} else {
				con <- gzfile(seqs, "rb")
			}
			on.exit(close(con), add=TRUE)
		} else if (inherits(seqs, "connection")) {
			con <- seqs
			if (!isOpen(con)) {
				open(con, type="rb")
				on.exit(close(con), add=TRUE)
			}
		} else {
			stop("seqs must be a file path or connection.")
		}
		
		newSeqs <- 0
		buffer <- ""
		enter <- TRUE
		newline <- FALSE
		while (enter) {
			if (verbose) {
				it <- it + 1L
				.cat(paste("\rReading FASTA file chunk ",
					it,
					sep=""))
			}
			
			# read in chunkSize characters
			r <- readChar(con=con,
				nchars=chunkSize,
				useBytes=TRUE)
			if (length(r)==0L) {
				# reached the end of the connection
				if (buffer=="") {
					break
				} else {
					enter <- FALSE
				}
			} else if (nchar(r) < chunkSize) {
				# need to ensure that the end of the connection was reached
				# because reading can terminate unexpectedly before chunkSize
				r2 <- readChar(con=con,
					nchars=chunkSize - nchar(r),
					useBytes=TRUE)
				while (length(r2) > 0 && nchar(r) < chunkSize) {
					r <- paste(r, r2, sep="")
					r2 <- readChar(con=con,
						nchars=chunkSize - nchar(r),
						useBytes=TRUE)
				}
				
				if (nchar(r) < chunkSize)
					enter <- FALSE
			}
			
			if (buffer != "")
				r <- paste(buffer, r, sep=ifelse(newline, "\n", ""))
			newline <- substr(r, nchar(r), nchar(r))=="\n"
			r <- gsub("\r\n", "\n", r, fixed=TRUE)
			r <- gsub("\r", "\n", r, fixed=TRUE)
			r <- strsplit(r, "\n", fixed=TRUE)[[1]]
			
			# descriptions contains the line index of each sequence
			descriptions <- which(substr(r, 1L, 1L)==">")
			if (length(descriptions)==0)
				stop("No FASTA records found.")
			
			# create new vectors by adding or subtracting 1
			# from each line index in description
			dp <- descriptions + 1L
			dm <- descriptions - 1L
			end <- c(dm[-1], length(r))
			
			# remove the last incomplete sequence unless it is the end of file
			if (enter) {
				buffer <- paste(r[descriptions[length(descriptions)]:length(r)],
					collapse="\n")
				length(descriptions) <- length(descriptions) - 1L
			}
			
			# numF contains the number of sequences for this iteration
			numF <- length(descriptions)
			if (numF==0)
				next
			
			# build the data frame sequence by sequence
			myData <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
				identifier=identifier)
			
			myData$description <- substr(r[descriptions],
				2L,
				nchar(r[descriptions]))
			
			sequence <- .Call("collapse",
				r,
				dp[seq_len(numF)],
				end[seq_len(numF)],
				PACKAGE="DECIPHER")
			sequence <- gsub(" ", "", sequence, fixed=TRUE)
			myData_ <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
				sequence=I(Codec(sequence,
					processors=processors,
					...)),
				quality=I(rep(list(raw()), length(sequence))))
			
			# add database columns to the data frame
			if (length(f) > 0) {
				for (i in 1:length(f)) {
					if (is.na(match(f[i], names(myData)))) {
						d <- data.frame(rep(NA, numF))
						names(d) <- f[i]
						myData <- data.frame(myData, d)
					}
				}
			}
			
			# numSeq contains the total number of sequences so far
			numSeq <- numSeq + length(descriptions)
			newSeqs <- newSeqs + length(descriptions)
			
			if (replaceTbl) {
				ft <- fts2
				ft_ <- fts1
			} else {
				ft <- ft_ <- NULL
			}
			
			dbWriteTable(dbConn,
				tblName,
				myData,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft)
			dbWriteTable(dbConn,
				paste("_", tblName, sep=""),
				myData_,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft_)
			
			replaceTbl <- FALSE
		}
	} else if (type==2) { # FASTQ
		# scan in the FASTQ file
		if (is.character(seqs)) {
			if (substr(seqs, 1, 7)=="http://" ||
				substr(seqs, 1, 8)=="https://" ||
				substr(seqs, 1, 7)=="ftps://" ||
				substr(seqs, 1, 6)=="ftp://") {
				con <- gzcon(url(seqs, "rb"))
			} else {
				con <- gzfile(seqs, "rb")
			}
			on.exit(close(con), add=TRUE)
		} else if (inherits(seqs, "connection")) {
			con <- seqs
			if (!isOpen(con)) {
				open(con, type="rb")
				on.exit(close(con), add=TRUE)
			}
		} else {
			stop("seqs must be a file path or connection.")
		}
		
		newSeqs <- 0
		buffer <- ""
		enter <- TRUE
		newline <- FALSE
		while (enter) {
			# scan piece of file into memory
			if (verbose) {
				it <- it + 1L
				.cat(paste("\rReading FASTQ file chunk ",
					it,
					sep=""))
			}
			
			# read in chunkSize characters
			r <- readChar(con=con,
				nchars=chunkSize,
				useBytes=TRUE)
			if (length(r)==0L) {
				# reached the end of the connection
				if (buffer=="") {
					break
				} else {
					enter <- FALSE
				}
			} else if (nchar(r) < chunkSize) {
				# need to ensure that the end of the connection was reached
				# because reading can terminate unexpectedly before chunkSize
				r2 <- readChar(con=con,
					nchars=chunkSize - nchar(r),
					useBytes=TRUE)
				while (length(r2) > 0 && nchar(r) < chunkSize) {
					r <- paste(r, r2, sep="")
					r2 <- readChar(con=con,
						nchars=chunkSize - nchar(r),
						useBytes=TRUE)
				}
				
				if (nchar(r) < chunkSize)
					enter <- FALSE
			}
			
			if (buffer != "")
				r <- paste(buffer, r, sep=ifelse(newline, "\n", ""))
			newline <- substr(r, nchar(r), nchar(r))=="\n"
			r <- gsub("\r\n", "\n", r, fixed=TRUE)
			r <- gsub("\r", "\n", r, fixed=TRUE)
			r <- strsplit(r, "\n", fixed=TRUE)[[1]]
			
			# descriptions contains the line index of each sequence
			descriptions <- which(substr(r, 1L, 1L)=="@")
			if (length(descriptions)==0)
				stop("No FASTQ records found.")
			descriptions <- descriptions[which((descriptions - descriptions[1]) %% 4==0)]
			
			# remove the last incomplete sequence unless it is the end of file
			if (enter) {
				buffer <- paste(r[descriptions[length(descriptions)]:length(r)],
					collapse="\n")
				length(descriptions) <- length(descriptions) - 1L
			}
			
			# numF contains the number of sequences for this iteration
			numF <- length(descriptions)
			if (numF==0)
				next
			
			# build the data frame sequence by sequence
			myData <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
				identifier=identifier)
			
			myData$description <- substr(r[descriptions],
				2L,
				nchar(r[descriptions]))
			
			myData_ <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
				sequence=I(Codec(r[descriptions + 1L],
					processors=processors,
					...)),
				quality=I(Codec(r[descriptions + 3L],
					compression=c("qbit", "gzip"),
					processors=processors)))
			
			# add database columns to the data frame
			if (length(f) > 0) {
				for (i in 1:length(f)) {
					if (is.na(match(f[i], names(myData)))) {
						d <- data.frame(rep(NA, numF))
						names(d) <- f[i]
						myData <- data.frame(myData, d)
					}
				}
			}
			
			# numSeq contains the total number of sequences so far
			numSeq <- numSeq + length(descriptions)
			newSeqs <- newSeqs + length(descriptions)
			
			if (replaceTbl) {
				ft <- fts2
				ft_ <- fts1
			} else {
				ft <- ft_ <- NULL
			}
			
			dbWriteTable(dbConn,
				tblName,
				myData,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft)
			dbWriteTable(dbConn,
				paste("_", tblName, sep=""),
				myData_,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft_)
			
			replaceTbl <- FALSE
		}
	} else if (type==3) { # GenBank
		# scan in the GenBank file
		if (is.character(seqs)) {
			if (substr(seqs, 1, 7)=="http://" ||
				substr(seqs, 1, 8)=="https://" ||
				substr(seqs, 1, 7)=="ftps://" ||
				substr(seqs, 1, 6)=="ftp://") {
				con <- gzcon(url(seqs, "rb"))
			} else {
				con <- gzfile(seqs, "rb")
			}
			on.exit(close(con), add=TRUE)
		} else if (inherits(seqs, "connection")) {
			con <- seqs
			if (!isOpen(con)) {
				open(con, type="rb")
				on.exit(close(con), add=TRUE)
			}
		} else {
			stop("seqs must be a file path or connection.")
		}
		
		newSeqs <- 0
		buffer <- ""
		enter <- TRUE
		newline <- FALSE
		while (enter) {
			# scan piece of file into memory
			if (verbose) {
				it <- it + 1L
				.cat(paste("\rReading GenBank file chunk ",
					it,
					sep=""))
			}
			
			# read in chunkSize characters
			r <- readChar(con=con,
				nchars=chunkSize,
				useBytes=TRUE)
			if (length(r)==0L) {
				# reached the end of the connection
				if (buffer=="") {
					break
				} else {
					enter <- FALSE
				}
			} else if (nchar(r) < chunkSize) {
				# need to ensure that the end of the connection was reached
				# because reading can terminate unexpectedly before chunkSize
				r2 <- readChar(con=con,
					nchars=chunkSize - nchar(r),
					useBytes=TRUE)
				while (length(r2) > 0 && nchar(r) < chunkSize) {
					r <- paste(r, r2, sep="")
					r2 <- readChar(con=con,
						nchars=chunkSize - nchar(r),
						useBytes=TRUE)
				}
				
				if (nchar(r) < chunkSize)
					enter <- FALSE
			}
			
			if (buffer != "")
				r <- paste(buffer, r, sep=ifelse(newline, "\n", ""))
			newline <- substr(r, nchar(r), nchar(r))=="\n"
			r <- gsub("\r\n", "\n", r, fixed=TRUE)
			r <- gsub("\r", "\n", r, fixed=TRUE)
			r <- strsplit(r, "\n", fixed=TRUE)[[1]]
			
			# descriptions contains the line index of each sequence
			descriptions <- which(substr(r, 1L, 5L)=="LOCUS")
			if (length(descriptions)==0)
				stop("No GENBANK records found.")
			seq_start <- which(substr(r, 1L, 6L)=="ORIGIN")
			seq_end <- which(substr(r, 1L, 2L)=="//") - 1L
			temp <- integer(length(seq_end))
			for (i in seq_along(seq_end)) {
				w <- which(seq_start > descriptions[i] & seq_start <= seq_end[i])
				if (length(w)==0) {
					temp[i] <- seq_end[i] + 1L # no sequence
				} else if (length(w)==1) {
					temp[i] <- seq_start[w] + 1L
				} else {
					stop("More than one ORIGIN found in GenBank record.")
				}
			}
			seq_start <- temp
			
			# remove the last incomplete sequence unless it is the end of file
			if (enter) {
				buffer <- paste(r[descriptions[length(descriptions)]:length(r)],
					collapse="\n")
				length(descriptions) <- length(descriptions) - 1L
				length(seq_start) <- length(seq_end) <- length(descriptions)
			}
			
			# numF contains the number of sequences for this iteration
			numF <- length(descriptions)
			if (numF==0)
				next
			
			# parse fields
			x <- .Call("extractFields",
				r,
				fields,
				descriptions,
				seq_start,
				PACKAGE="DECIPHER")
			names(x) <- names(fields)
			
			# build the data frame sequence by sequence
			myData <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
				identifier=identifier,
				x)
			
			sequence <- .Call("collapse",
				substr(r, 11, nchar(r)),
				seq_start[seq_len(numF)],
				seq_end[seq_len(numF)],
				PACKAGE="DECIPHER")
			sequence <- gsub(" ", "", sequence, fixed=TRUE)
			
			myData_ <- data.frame(row_names=seq(from=(numSeq + 1),
					to=(numSeq + length(descriptions))),
				sequence=I(Codec(sequence,
					processors=processors,
					...)),
				quality=I(rep(list(raw()), length(sequence))))
			
			# add database columns to the data frame
			if (length(f) > 0) {
				for (i in 1:length(f)) {
					if (is.na(match(f[i], names(myData)))) {
						d <- data.frame(rep(NA, numF))
						names(d) <- f[i]
						myData <- data.frame(myData, d)
					}
				}
			}
			
			# numSeq contains the total number of sequences so far
			numSeq <- numSeq + length(descriptions)
			newSeqs <- newSeqs + length(descriptions)
			
			if (replaceTbl) {
				ft <- fts2
				ft_ <- fts1
			} else {
				ft <- ft_ <- NULL
			}
			
			dbWriteTable(dbConn,
				tblName,
				myData,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft)
			dbWriteTable(dbConn,
				paste("_", tblName, sep=""),
				myData_,
				row.names=FALSE,
				overwrite=replaceTbl,
				append=!replaceTbl,
				field.types=ft_)
			
			replaceTbl <- FALSE
		}
	} else { # XStringSet or QualityScaledXStringSet
		# add the sequences to the database
		newSeqs <- length(seqs)
		if (length(newSeqs)==0)
			stop("No sequences in seqs.")
		
		if (verbose) {
			cat("Adding", newSeqs, "sequences to the database.")
			flush.console()
		}
		
		# give the seqs names if they do not have any
		if (is.null(names(seqs)))
			names(seqs) <- (numSeq + 1):(numSeq + newSeqs)
		
		# build the data frame sequence by sequence
		myData <- data.frame(row_names=seq(from=(numSeq + 1),
			to=(numSeq + newSeqs)),
			identifier=identifier,
			description=names(seqs))
		
		if (type > 8) {
			quality <- Codec(as.character(quality(seqs)),
				compression=c("qbit", "gzip"),
				processors=processors)
		} else {
			quality <- rep(list(raw()), length(seqs))
		}
		
		myData_ <- data.frame(row_names=seq(from=(numSeq + 1),
			to=(numSeq + newSeqs)),
			sequence=I(Codec(as.character(seqs),
				processors=processors,
				...)),
			quality=I(quality))
		
		# add database columns to the data frame
		if (length(f) > 0) {
			for (i in 1:length(f)) {
				if (is.na(match(f[i], names(myData)))) {
					d <- data.frame(rep(NA, newSeqs))
					names(d) <- f[i]
					myData <- data.frame(myData, d)
				}
			}
		}
		
		if (replaceTbl) {
			ft <- fts2
			ft_ <- fts1
		} else {
			ft <- ft_ <- NULL
		}
		
		dbWriteTable(dbConn,
			tblName,
			myData,
			row.names=FALSE,
			overwrite=replaceTbl,
			append=!replaceTbl,
			field.types=ft)
		dbWriteTable(dbConn,
			paste("_", tblName, sep=""),
			myData_,
			row.names=FALSE,
			overwrite=replaceTbl,
			append=!replaceTbl,
			field.types=ft_)
	}
	
	searchExpression <- paste("select count(*) from ",
		dbQuoteIdentifier(dbConn, tblName),
		sep="")
	numSeq <- as.double(dbGetQuery(dbConn, searchExpression)[[1]])
	
	if (verbose) { # print the elapsed time to import
		time.2 <- Sys.time()
		cat("\n")
		if (newSeqs != numSeq)
			cat("\nAdded ",
				newSeqs,
				" new sequence",
				ifelse(newSeqs != 1, "s", ""),
				" to table ",
				tblName,
				".",
				sep="")
		cat("\n",
			numSeq,
			" total sequence",
			ifelse(numSeq != 1, "s", ""),
			" in table ",
			tblName,
			".",
			"\n",
			sep="")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
		flush.console()
	}
	invisible(numSeq)
}
