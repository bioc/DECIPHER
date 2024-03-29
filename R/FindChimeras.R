FindChimeras <- function(dbFile,
	tblName="Seqs",
	identifier="",
	dbFileReference,
	tblNameReference="Seqs",
	batchSize=100,
	minNumFragments=20000,
	tb.width=5,
	multiplier=20,
	minLength=30,
	minCoverage=0.6,
	overlap=100,
	minSuspectFragments=4,
	showPercentCoverage=FALSE,
	add2tbl=FALSE,
	maxGroupSize=-1,
	minGroupSize=25,
	excludeIDs=NULL,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (!is.character(tblNameReference))
		stop("tblNameReference must be a character string.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.null(excludeIDs) && !is.character(excludeIDs))
		stop("excludeIDs must be a character string.")
	if (!is.numeric(batchSize))
		stop("batchSize must be a numeric.")
	if (!is.numeric(minNumFragments))
		stop("minNumFragments must be a numeric.")
	if (!is.numeric(tb.width))
		stop("tb.width must be a numeric.")
	if (tb.width >= 15)
		stop("tb.width must be less than 15.")
	if (tb.width <= 0)
		stop("tb.width must be greater than zero.")
	if (!is.numeric(multiplier))
		stop("multiplier must be a numeric.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.logical(add2tbl) && !is.character(add2tbl))
		stop("add2tbl must be a logical or table name.")
	if (!is.numeric(minLength))
		stop("minLength must be a numeric.")
	if (minLength < 0)
		stop("minLength must be positive or zero.")
	if (!is.numeric(overlap))
		stop("overlap must be a numeric.")
	if (overlap < 0)
		stop("overlap must be positive or zero.")
	if (!is.numeric(minCoverage))
		stop("minCoverage must be a numeric.")
	if (minCoverage < 0 || minCoverage > 1)
		stop("minCoverage must be between 0 and 1 inclusive.")
	if (!is.numeric(minSuspectFragments))
		stop("minSuspectFragments must be a numeric.")
	if (minSuspectFragments < 0)
		stop("minSuspectFragments must be positive or zero.")
	if (!is.numeric(maxGroupSize))
		stop("maxGroupSize must be a numeric.")
	if (maxGroupSize == 0)
		stop("maxGroupSize must be non-zero.")
	if (maxGroupSize != floor(maxGroupSize))
		stop("maxGroupSize must be an integer.")
	if (!is.numeric(minGroupSize))
		stop("minGroupSize must be a numeric.")
	if (minGroupSize < 0)
		stop("minGroupSize must be at least zero.")
	if (minGroupSize != floor(minGroupSize))
		stop("minGroupSize must be an integer.")
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
	
	if (verbose)
		time.1 <- Sys.time()
	
	# initialize databases
	if (is.character(dbFile)) {
		if (!requireNamespace("RSQLite", quietly=TRUE))
			stop("Package 'RSQLite' must be installed.")
		dbConn1 <- dbConnect(dbDriver("SQLite"), dbFile)
		on.exit(dbDisconnect(dbConn1))
	} else {
		dbConn1 <- dbFile
		if (!dbIsValid(dbConn1))
			stop("The connection has expired.")
	}
	if (is.character(dbFileReference)) {
		if (!requireNamespace("RSQLite", quietly=TRUE))
			stop("Package 'RSQLite' must be installed.")
		dbConn2 <- dbConnect(dbDriver("SQLite"), dbFileReference)
		on.exit(dbDisconnect(dbConn2), add=TRUE)
	} else {
		dbConn2 <- dbFileReference
		if (!dbIsValid(dbConn2))
			stop("The connection has expired.")
	}
	
	searchExpression <- paste("select distinct ",
		dbQuoteIdentifier(dbConn2, "identifier"),
		", ",
		dbQuoteIdentifier(dbConn2, "origin"),
		" from ",
		dbQuoteIdentifier(dbConn2, tblNameReference),
		sep="")
	rs <- dbSendQuery(dbConn2, searchExpression)
	searchResult <- dbFetch(rs, n=-1, row.names=FALSE)
	groups <- searchResult$identifier
	origins <- searchResult$origin
	dbClearResult(rs)
	
	searchExpression <- paste("select distinct",
		dbQuoteIdentifier(dbConn1, "identifier"),
		"from",
		dbQuoteIdentifier(dbConn1, tblName))
	rs <- dbSendQuery(dbConn1, searchExpression)
	searchResult <- dbFetch(rs, n=-1, row.names=FALSE)
	myGroups <- searchResult$identifier
	dbClearResult(rs)
	
	if (identifier[1] != "") {
		w <- which(!(identifier %in% groups))
		if (length(w) > 0) {
			warning("Identifier(s) not found in dbFile: ", paste(unique(identifier[w]), sep=", "))
			identifier <- identifier[-w]
		}
		myGroups <- unique(identifier)
	}
	w <- which(!(myGroups %in% groups))
	if (length(w) > 0) {
		warning("Identifier(s) not found in dbFileReference: ", paste(unique(myGroups[w]), sep=", "))
		myGroups <- myGroups[-w]
	}
	w <- which(myGroups %in% excludeIDs)
	if (length(w) > 0)
		myGroups <- myGroups[-w]
	
	# determine if the reference db contains a chimera column
	f <- dbListFields(dbConn2, tblNameReference)
	firstRound2 <- is.na(match("chimera", f))
	
	# determine if the dbFile contains a chimera column
	f <- dbListFields(dbConn1, tblName)
	firstRound1 <- is.na(match("chimera", f))
	
	# remove groups with too few sequences
	remove <- character()
	for (group in myGroups) {
		if (firstRound2) {
			numG <- SearchDB(dbConn2,
				tblName=tblNameReference,
				identifier=group,
				processors=processors,
				verbose=FALSE,
				clause=paste(dbQuoteIdentifier(dbConn2, "nonstandard"), "< 20"),
				countOnly=TRUE)
		} else {
			numG <- SearchDB(dbConn2,
				identifier=group,
				tblName=tblNameReference,
				processors=processors,
				verbose=FALSE,
				clause=paste(dbQuoteIdentifier(dbConn2, "chimera"),
					"is NULL and",
					dbQuoteIdentifier(dbConn2, "nonstandard"),
					"< 20"),
				countOnly=TRUE)
		}
		if (numG < minGroupSize) # too small
			remove <- c(remove, group)
	}
	w <- which(myGroups %in% remove)
	if (length(w) > 0)
		myGroups <- myGroups[-w]
	
	# initialize fragment variables
	offset <- 5
	min_hits_in <- 5
	fragments <- DNAStringSet()
	all_hits_in <- integer()
	all_results <- NULL
	percentCoverage <- ""
	m <- integer()
	rns <- integer()
	
	if (verbose)
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	
	for (group in myGroups) {
		if (firstRound2) {
			group_dna <- SearchDB(dbConn2,
				tblName=tblNameReference,
				identifier=group,
				processors=processors,
				verbose=FALSE,
				type="DNAStringSet",
				limit=maxGroupSize,
				removeGaps="all",
				clause=paste(dbQuoteIdentifier(dbConn2, "nonstandard"), "< 20"))
		} else {
			group_dna <- SearchDB(dbConn2,
				tblName=tblNameReference,
				identifier=group,
				processors=processors,
				verbose=FALSE,
				type="DNAStringSet",
				limit=maxGroupSize,
				removeGaps="all",
				clause=paste(dbQuoteIdentifier(dbConn2, "chimera"),
					"is NULL and",
					dbQuoteIdentifier(dbConn2, "nonstandard"),
					"< 20"))
		}
		numG <- length(group_dna)
		
		if (firstRound1) {
			dna <- SearchDB(dbConn1,
				tblName=tblName,
				identifier=group,
				type="DNAStringSet",
				replaceChar="",
				removeGaps="all",
				processors=processors,
				verbose=FALSE)
		} else {
			dna <- SearchDB(dbConn1,
				tblName=tblName,
				identifier=group,
				type="DNAStringSet",
				replaceChar="",
				removeGaps="all",
				processors=processors,
				verbose=FALSE,
				clause=paste(dbQuoteIdentifier(dbConn1, "chimera"), "is NULL"))
		}
		
		# reduce to the set of unique dna sequences
		u_dna <- unique(dna)
		names(u_dna) <- names(dna)[match(u_dna, dna)]
		m <- c(m, as.integer(names(u_dna)[match(dna, u_dna)]))
		rns <- c(rns, as.integer(names(dna)))
		dna <- u_dna
		numF <- length(dna)
		
		batches <- seq(1, numF, batchSize)
		
		# for all sequences
		for (j in batches) {
			# subset of the whole group's sequence set
			seqs <- j:(j + batchSize - 1)
			# remove sequences beyond length of sequence set
			# (only applies to the final iteration)
			w <- which(seqs > numF)
			if (length(w) > 0)
				seqs <- seqs[-w]
			
			# split each sequence into fragments of length 30
			temp_fragments <- DNAStringSet()
			count <- 0
			widths <- width(dna[seqs])
			numFragments <- floor((max(widths) - 30)/offset) + 1
			for (i in seq_len(numFragments)) {
				start <- ((i - 1)*offset + 1)
				end <- ((i - 1)*offset + 30)
				w <- which(widths >= end)
				l <- length(w)
				temp_fragments <- c(temp_fragments,
					subseq(dna[seqs[w]], start, end))
				names(temp_fragments)[seq((count + 1), (count + l))] <-
					paste(group, ".", names(dna[seqs[w]]),
						".", rep(i, l), ".", widths[w], sep="")
					count <- count + l
			}
			
			# remove fragments with a trusted-band containing ambiguities
			# first 5 nucleotide trusted-band
			w <- which(rowSums(alphabetFrequency(subseq(temp_fragments,
				1,
				2*tb.width))[,5:15]) == 0)
			if (length(w) > 0)
				temp_fragments <- temp_fragments[w]
			
			if (length(temp_fragments) > 0) {
				# remove fragments with a middle containing Ns
				w <- which(alphabetFrequency(subseq(temp_fragments,
					2*tb.width + 1,
					30))[,15] == 0)
				if (length(w) > 0)
					temp_fragments <- temp_fragments[w]
			}
			
			if (length(temp_fragments) > 0) {
				
				u_temp_fragments <- unique(temp_fragments)
				m_fragments <- match(temp_fragments, u_temp_fragments)
				pdict <- PDict(u_temp_fragments,
					algorithm="ACtree2",
					tb.end=tb.width)
				
				if (numG <= 3000) { # search whole group
					# count number of hits in group
					counts <- vwhichPDict(pdict,
						group_dna,
						max.mismatch=1,
						fixed=FALSE)
					counts <- unlist(counts)
					
					# tabulate counts
					if (length(counts) > 0) { # hits
						t <- tabulate(counts)
						if (length(t) < length(u_temp_fragments))
							t[(length(t) + 1):length(u_temp_fragments)] <- 0
						hits_in <- t[m_fragments]
						
						# remove fragments with too many hits in group
						w <- which(hits_in < min_hits_in)
					} else { # no hits
						hits_in <- integer(length(temp_fragments))
						w <- seq_along(temp_fragments)
					}
				} else { # search part of group first for speed
					# count number of hits in part of group
					counts <- vwhichPDict(pdict,
						group_dna[1:2999],
						max.mismatch=1,
						fixed=FALSE)
					counts <- unlist(counts)
					
					# tabulate counts
					if (length(counts) > 0) { # hits
						t <- tabulate(counts)
						if (length(t) < length(u_temp_fragments))
							t[(length(t) + 1):length(u_temp_fragments)] <- 0
						hits_in <- t[m_fragments]
						
						# remove fragments with too many hits in group
						w <- which(hits_in < min_hits_in)
					} else { # no hits
						hits_in <- integer(length(temp_fragments))
						w <- seq_along(temp_fragments)
					}
					
					hits_in <- hits_in[w]
					temp_fragments <- temp_fragments[w]
					
					if (length(temp_fragments) > 0) {
						u_temp_fragments <- unique(temp_fragments)
						m_fragments <- match(temp_fragments, u_temp_fragments)
						
						# count number of hits in group
						pdict <- PDict(u_temp_fragments,
							algorithm="ACtree2",
							tb.end=tb.width)
						counts <- vwhichPDict(pdict,
							group_dna[3000:numG],
							max.mismatch=1,
							fixed=FALSE)
						counts <- unlist(counts)
						
						# tabulate counts
						if (length(counts) > 0) { # hits
							t <- tabulate(counts)
							if (length(t) < length(u_temp_fragments))
								t[(length(t) + 1):length(u_temp_fragments)] <- 0
							hits_in2 <- t[m_fragments]
							hits_in <- hits_in + hits_in2
						
							# remove fragments with too many hits in group
							w <- which(hits_in < min_hits_in)
						} else { # no hits
							w <- seq_along(temp_fragments)
						}
					}
				}
			} else { # no temp_fragments
				w <- integer()
			}
			
			if (length(w) > 0) {
				temp_fragments <- temp_fragments[w]
				u_temp_fragments <- unique(temp_fragments)
				m_fragments <- match(temp_fragments, u_temp_fragments)
				hits_in <- hits_in[w]
				
				# count number of lenient hits in group
				pdict <- PDict(u_temp_fragments,
					algorithm="ACtree2",
					tb.start=tb.width + 1,
					tb.end=2*tb.width)
				counts <- vwhichPDict(pdict,
					group_dna,
					max.mismatch=2,
					fixed=FALSE)
				counts <- unlist(counts)
				
				# tabulate counts
				if (length(counts) == 0) { # no hits
					w <- integer() # should never occurr
				} else { # hits
					t <- tabulate(counts)
					if (length(t) < length(u_temp_fragments))
						t[(length(t) + 1):length(u_temp_fragments)] <- 0
					lenient_hits_in <- t[m_fragments]
					
					# remove fragments with too many lenient hits in group
					w <- which(lenient_hits_in < 10)
				}
			}
			
			if (length(w) > 0) {
				# remove fragments that would not meet requirements for a chimera
				name.position.width <- unlist(strsplit(names(temp_fragments[w]),
					".",
					fixed=TRUE))
				l <- length(name.position.width)
				
				# get the sequences of origin
				name <- name.position.width[seq(2, l, 4)]
				unique_name <- unique(name)
				
				# get the fragment positions
				position <- as.numeric(name.position.width[seq(3, l, 4)])
				
				# get the original sequence widths
				widths <- as.numeric(name.position.width[seq(4, l, 4)])
				
				remove <- character()
				for (i in unique_name) {
					w1 <- which(name == i)
					
					if (length(w1) < minSuspectFragments) {
						remove <- c(remove, i)
						next # chimeras need at least minSuspectFragments fragments
					}
					
					start <- ((position[w1[1]] - 1)*offset + 1)
					end <- ((position[w1[length(w1)]] - 1)*offset + 30)
					
					# chimeric section needs to be at least minLength long
					if ((end - start + 1) < minLength) {
						remove <- c(remove, i)
						next
					}
					
					# chimeric section must overlap ends
					if (!(end > (widths[w1[1]] - overlap)) &&
						!(start < overlap))
						remove <- c(remove, i)
				}
				
				if (length(remove) > 0)
					w <- w[-which(name %in% remove)]
				
				if (length(w) > 0) {
					fragments <- c(fragments, temp_fragments[w])
					all_hits_in <- c(all_hits_in, hits_in[w])
				}
			}
			
			# keep adding fragments until above minimum
			if (length(unique(fragments)) < minNumFragments &&
				(!(j == batches[length(batches)] && group == myGroups[length(myGroups)])
				|| length(fragments) == 0)) # and not the last batch
				next
			
			# search through other groups for fragments
			ns <- unlist(strsplit(names(fragments), ".", fixed=TRUE))
			groupsOfOrigin <- ns[seq(1, length(ns), 4)]
			n <- ns[seq(2, length(ns), 4)]
			un <- unique(n)
			result <- data.frame(chimera=I(character(length(un))),
				row.names=as.character(un))
			
			# for each group where the fragment did not originate
			if (verbose) {
				beg <- j - 1
				fin <- j + batchSize
				if (fin > numF)
					fin <- numF
				count <- 0
			}
			for (other in groups) {
				if (verbose)
					count <- count + 1
				if (other %in% excludeIDs)
					next # don't allow chimeras
				
				# remove groups that are from the other group and
				# remove groups that share unclassified's line of descent
				uGroupsOfOrigin <- unique(groupsOfOrigin)
				# find origin_other (line of descent of other)
				origin_other <- origins[match(other, groups)]
				g2 <- grepl("unclassified", other, fixed=TRUE)
				for (g in uGroupsOfOrigin) {
					if (other == g) {
						uGroupsOfOrigin <- uGroupsOfOrigin[-match(g, uGroupsOfOrigin)]
					} else if (g2 ||
						grepl("unclassified", g, fixed=TRUE)) {
						# find origin (line of descent of group)
						origin <- origins[match(g, groups)]
						if (grepl(origin, origin_other, fixed=TRUE) ||
							grepl(origin_other, origin, fixed=TRUE))
							uGroupsOfOrigin <- uGroupsOfOrigin[-match(g, uGroupsOfOrigin)]
					}
				}
				
				# do not search the group of origin
				h <- which(groupsOfOrigin %in% uGroupsOfOrigin)
				if (length(h) == 0)
					next
				
				u_fragments <- unique(fragments[h])
				m_fragments <- match(fragments[h], u_fragments)
				pdict <- PDict(u_fragments,
					algorithm="ACtree2",
					tb.end=tb.width)
				
				# search for the fragments
				if (firstRound2) {
					other_dna <- SearchDB(dbConn2,
						tblName=tblNameReference,
						identifier=other,
						processors=processors,
						verbose=FALSE,
						type="DNAStringSet",
						limit=maxGroupSize,
						removeGaps="all",
						clause=paste(dbQuoteIdentifier(dbConn2, "nonstandard"), "< 20"))
				} else {
					other_dna <- SearchDB(dbConn2,
						tblName=tblNameReference,
						identifier=other,
						processors=processors,
						verbose=FALSE,
						type="DNAStringSet",
						limit=maxGroupSize,
						removeGaps="all",
						clause=paste(dbQuoteIdentifier(dbConn2, "chimera"),
							"is NULL and",
							dbQuoteIdentifier(dbConn2, "nonstandard"),
							"< 20"))
				}
				
				if (length(other_dna) <= multiplier)
					next # not enough sequences
				
				counts <- vwhichPDict(pdict,
					other_dna,
					max.mismatch=1,
					fixed=FALSE)
				counts <- unlist(counts)
				
				# tabulate counts
				if (length(counts) == 0) # no hits
					next
				t <- tabulate(counts)
				if (length(t) < length(u_fragments))
					t[(length(t) + 1):length(u_fragments)] <- 0
				hits_out <- t[m_fragments]
				
				# mark fragments with hits out of group greater than
				# multiplier times the number of hits in group
				# as possible chimeras from the other group
				x <- all_hits_in[h]*multiplier < hits_out
				y <- hits_out > multiplier
				w <- which(x & y)
				if (length(w) > 0) {
					# find all fragments corresponding to the unique fragment
					name.position.width <- unlist(strsplit(names(fragments[h[w]]),
						".",
						fixed=TRUE))
					l <- length(name.position.width)
					
					# get the sequences of origin
					name <- name.position.width[seq(2, l, 4)]
					unique_name <- unique(name)
					
					# get the fragment positions
					position <- as.numeric(name.position.width[seq(3, l, 4)])
					
					# get the original sequence widths
					widths <- as.numeric(name.position.width[seq(4, l, 4)])
					
					for (i in unique_name) {
						w1 <- which(name == i)
						
						if (length(w1) < minSuspectFragments)
							next # chimeras need at least minSuspectFragments fragments
						
						# needs to have a high density of fragments
						start <- 0
						end <- 0
						coverage <- 0
						for (p in position[w1]) {
							newCoverage <- (p - 1)*offset + 30 - end
							if (newCoverage > 30)
								coverage <- coverage + 30
							else
								coverage <- coverage + newCoverage
							end <- ((p - 1)*offset + 30)
						}
						
						start <- ((position[w1[1]] - 1)*offset + 1)
						end <- ((position[w1[length(w1)]] - 1)*offset + 30)
						
						# chimeric section needs > % coverage of fragments
						if (coverage <= minCoverage*(end - start + 1))
							next
						
						# chimeric section needs to be at least minLength long
						if ((end - start + 1) < minLength)
							next
						
						# chimeric section must overlap ends
						if (!(end > (widths[w1[1]] - overlap)) &&
							!(start < overlap))
							next
						
						index <- match(i,row.names(result))
						
						# append this group to the result
						w2 <- which(result$chimera[index] == "")
						if (length(w2) > 0) {
							result$chimera[index[w2]] <- other
							
							if (length(index[-w2]) > 0)
								result$chimera[index[-w2]] <- paste(result$chimera[index[-w2]],
									other,
									sep=", ")
						} else {
							result$chimera[index] <- paste(result$chimera[index],
								other,
								sep=", ")
						}
						
						if (showPercentCoverage)
							percentCoverage <- paste(":",
								round(100*coverage/(end - start + 1),
									digits=2),
								"%",
								sep="")
						
						# append chimeric positions to the result
						result$chimera[index] <- paste(result$chimera[index],
							" (",
							paste(as.character(start),
								"-",
								as.character(end),
								percentCoverage,
								")",
								sep=""),
							sep="")
					}
				}
				if (verbose)
					setTxtProgressBar(pBar,
						(beg + (fin - beg)*count/length(groups))/numF)
			}
			
			# initialize fragment variables
			fragments <- DNAStringSet()
			all_hits_in <- integer()
			
			# update the result in the database
			result <- subset(result, result$chimera != "")
			if (dim(result)[1] > 0) {
				if (is.null(all_results))
					all_results <- result
				else
					all_results <- rbind(all_results, result)
			}
		}
	}
	
	# expand results to include duplicated sequences
	rows <- match(m, as.integer(rownames(all_results)))
	w <- which(!is.na(rows))
	if (length(w) > 0) {
		rows <- rows[w]
		if (length(rows) > dim(all_results)[1])
			all_results <- data.frame(chimera=I(all_results$chimera[rows]),
				row.names=as.character(rns[w]))
	}
	
	if ((is.character(add2tbl) || add2tbl) && !is.null(all_results))
		Add2DB(all_results,
			dbConn1,
			tblName=ifelse(is.character(add2tbl), add2tbl, tblName),
			verbose=FALSE)
	
	if (verbose) {
		setTxtProgressBar(pBar, 1)
		close(pBar)
		if (is.null(all_results)) {
			cat("\nNo chimeras found.")
		} else {
			cat("\nFound",
				dim(all_results)[1],
				ifelse(dim(all_results)[1] > 1,
					"possible chimeras.",
					"possible chimera."))
			if ((is.character(add2tbl) || add2tbl) && !is.null(all_results))
				cat("\nAdded to table ",
					ifelse(is.character(add2tbl), add2tbl, tblName),
					":  \"chimera\".",
					sep="")
		}
		
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(all_results)
}
