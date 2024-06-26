BrowseDB <- function(dbFile,
	htmlFile=tempfile(fileext=".html"),
	openURL=interactive(),
	tblName="Seqs",
	identifier="",
	limit=-1,
	orderBy="row_names",
	maxChars=50,
	title="",
	clause="") {
	
	# error checking
	if (!is.character(htmlFile))
		if (!inherits(htmlFile, "connection"))
			stop("htmlFile must be a character string or connection.")
	if (!is.logical(openURL) || is.na(openURL))
		stop("openURL must be TRUE or FALSE.")
	if (!is.character(tblName))
		stop("tblName must be a character string.")
	if (startsWith(tblName, "_"))
		stop("Invalid tblName.")
	if (!is.character(identifier))
		stop("identifier must be a character string.")
	if (!is.character(orderBy))
		stop("orderBy must be a character string.")
	if (is.numeric(limit)) {
		if (floor(limit) != limit)
				stop("limit must be a whole number or two comma-separated whole numbers specifying offset,limit.")
	} else {
		if (!grepl("[0-9],[0-9]", limit, perl=TRUE)) {
			limit <- as.numeric(limit)
			if (floor(limit) != limit)
				stop("limit must be a whole number or two comma-separated whole numbers specifying offset,limit.")
		}
	}
	if (!is.numeric(maxChars))
		stop("maxChars must be a numeric.")
	if (floor(maxChars) != maxChars)
		stop("maxChars must be a whole number.")
	if (maxChars <= 0)
		stop("maxChars must be greater than zero.")
	if (!is.character(clause))
		stop("clause must be a character string.")
	if (length(title) != 1L)
		stop("title must be length one.")
	
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
	
	# find the number of rows
	searchExpression <- paste('select count(*) from',
		dbQuoteIdentifier(dbConn, tblName))
	if (identifier != "")
		searchExpression <- paste(searchExpression,
			' where ',
			dbQuoteIdentifier(dbConn, "identifier"),
			' = "',
			identifier,
			'"',
			sep="")
	if (clause != "")
		searchExpression <- paste(searchExpression,
			ifelse(identifier == "", "where", "and"),
			clause)
	if (orderBy != "row_names") # default ordering is row_names
		searchExpression <- paste(searchExpression,
			' order by ',
			dbQuoteIdentifier(dbConn, tblName),
			'.',
			dbQuoteIdentifier(dbConn, orderBy),
			sep='')
	
	rs <- dbSendQuery(dbConn, searchExpression)
	count <- as.double(dbFetch(rs, n=-1, row.names=FALSE)[[1]])
	dbClearResult(rs)
	
	# count is the numer of rows in the table body
	if (is.character(limit)) {
		temp <- as.numeric(strsplit(limit, ",", fixed=TRUE)[[1]])
		count <- min(count - temp[1], temp[2])
	} else if (count > limit && # use limit
		limit > 0) { # limit is positive
		count <- limit
	}
	if (is.na(count) || count <= 0)
		stop("No results matched the specified constraints.")
	html <- character(count)
	
	# gives the table rows alternating colors
	for (j in 1:count) {
		if ((j %% 2) == 1) {
			html[j] <- "<tr class=\"row1\">"
		} else {
			html[j] <- "<tr class=\"row2\">"
		}
	}
	
	# gets all the fields available
	field <- dbListFields(dbConn, tblName)
	header <- ""
	tds <- ""
	tableWidth <- 0L
	for (i in seq_along(field)) {
		w <- width(field[i])
		
		# build the search expression
		searchExpression <- paste('select',
			dbQuoteIdentifier(dbConn, field[i]),
			'from',
			dbQuoteIdentifier(dbConn, tblName))
		if (identifier != "")
			searchExpression <- paste(searchExpression,
				' where ',
				dbQuoteIdentifier(dbConn, "identifier"),
				' = "',
				identifier,
				'"',
				sep="")
		if (clause != "")
			searchExpression <- paste(searchExpression,
				ifelse(identifier == "", " where ", " and "),
				clause,
				sep="")
		if (orderBy != "row_names") # default ordering is row_names
			searchExpression <- paste(searchExpression,
				' order by ',
				dbQuoteIdentifier(dbConn, tblName),
				'.',
				dbQuoteIdentifier(dbConn, orderBy),
				sep='')
		if (limit > 0)
			searchExpression <- paste(searchExpression,
				'limit',
				limit)
		
		rs <- dbSendQuery(dbConn, searchExpression)
		searchResult <- dbFetch(rs, n=-1, row.names=FALSE)
		dbClearResult(rs)
		
		# NAs in R are NULLs in SQL: replace missing values
		searchResult[is.na(searchResult)] <- "NULL"
		
		# replace newlines with html breaks
		searchResult[, 1] <- gsub("\n",
			"<br>",
			substring(searchResult[, 1], 1, maxChars),
			fixed=TRUE)
		
		if (max(width(searchResult[,1])) > w)
			w <- max(width(searchResult[,1]))
		
		html <- paste(html,
			"<td class=\"td",
			i,
			"\">",
			searchResult[,1],
			"</td>",
			sep="")
		
		# use the multiplier to convert chars to pixels
		# there are two multipliers for better looks
		if (w > 10) {
			multiplier <- 9
		} else {
			multiplier <- 10
		}
		
		# build table header row
		header <- paste(header,
			"<td class=\"td",
			i,
			"\" style=\"font-weight:bold;\">",
			field[i],
			"</td>",
			sep="")
		
		# build the table column heights
		tds <- paste(tds,
			"\n.td",
			i,
			" {width:",
			w*multiplier,
			"px;}",
			sep="")
		tableWidth <- tableWidth + w*multiplier
	}
	html <- paste(html,
		"</tr>",
		sep="")
	
	# the html is arranged into two tables
	# the first table is for the header
	# the second is for the data and scrolls
	# they have the same column widths for alignment
	
	styles <- paste("\n<style type=text/css>",
		# creates the alternate row background coloring
		"\ntr.row1 td {background-color: #edf3fe; color: #000;}",
		"\ntr.row2 td {background-color: #FFFFFF; color: #000;}",
		# aligns the text in the center of the cell
		# and sets the cell height
		"\ntr td {height: 20px; text-align:center;}",
		# causes the scroll bar to always show for alignment reasons
		# sets the table at 90% height because that prevents overflow
		"\n.tbody {overflow-y: scroll; height: 90%; width: ",
		# sets the table width to the sum of column widths
		tableWidth,
		# gives the headers an invisible padding
		"px; border-style: solid; border-color: #EEEEEE}",
		# removes the cell borders and sets the font
		"\ntable {border-collapse: collapse; font-size: 14px; font-family: monospace;}",
		"\n.thead{text-align: center; width:",
		tableWidth - 16,
		"px; border-style: solid; border-color: #FFFFFF}",
		# makes the cell not have any padding that throws off alignment
		# and positions the header so that it will be fixed
		"\n.thead div{padding: 0px; float:left;}",
		# positions the table so that it can scroll independently
		"\n.table_div {float: left; height: 100%}",
		"\ntd {padding: 0px;}",# border: 1px solid red;
		tds,
		"\n</style>")
	
	html <- c("<html>",
		styles,
		"<head><title>",
		title,
		"</title></head>",
		title,
		# creates a table just for the header row so that it aligns
		"\n<div class=\"table_div\"><div class=\"thead\"><table>",
		header,
		# clears styles so that the table will be on a new line
		# without this the table and header are on one line
		"</table></div><div style=\"clear: both;\"></div>",
		# makes a new div for the main table
		"\n<div class=\"tbody\"><table>",
		html,
		"\n</table></div></div></html>")
	writeLines(html, htmlFile)
	
	if (openURL)
		browseURL(path.expand(htmlFile))
	
	invisible(htmlFile)
}
