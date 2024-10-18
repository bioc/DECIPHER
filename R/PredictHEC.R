PredictHEC <- function(myAAStringSet,
	type="states",
	windowSize=8,
	background=c(H=-0.2, E=-0.6, C=0),
	HEC_MI1=NULL,
	HEC_MI2=NULL) {
	
	# error checking
	if (!is(myAAStringSet, "AAStringSet"))
		stop("myAAStringSet must be an AAStringSet.")
	TYPES <- c("states", "scores", "probabilities")
	if (length(type) == 0)
		stop("No type specified.")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type == -1)
		stop("Ambiguous type.")
	if (!is.numeric(windowSize))
		stop("windowSize must be a numeric.")
	if (any(floor(windowSize) != windowSize))
		stop("windowSize must be an integer number.")
	if (any(windowSize < 1))
		stop("windowSize must be at least 1.")
	if (length(windowSize) == 1L) {
		windowSize <- rep(windowSize, 2L)
	} else if (length(windowSize) != 2L) {
		stop("windowSize must be one or two numbers.")
	}
	windowSize <- as.integer(windowSize)
	if (!is.double(background))
		stop("background must be a numeric.")
	if (is(background, "matrix")) {
		N <- nrow(background)
		nms <- rownames(background)
		if (ncol(background) != length(myAAStringSet))
			stop("background must have as many columns as sequences in myAAStringSet.")
	} else {
		N <- length(background)
		nms <- names(background)
	}
	if (N < 2L)
		stop("The number of states in background must be at least 2.")
	if (length(nms) != N)
		stop("background must have state names.")
	if (sum(nchar(names(background)) != 1L) > 0L)
		stop("background must have single letter state names.")
	if (is.null(HEC_MI1)) {
		data("HEC_MI1", envir=environment(), package="DECIPHER")
	} else {
		if (!is.double(HEC_MI1))
			stop("HEC_MI1 must be an array of numerics.")
		if (length(dim(HEC_MI1)) != 3L)
			stop("HEC_MI1 must be a three dimensional array.")
		if (dim(HEC_MI1)[1L] != 20L ||
			windowSize[1L] > ((dim(HEC_MI1)[2L] - 1)/2) ||
			(dim(HEC_MI1)[2L] %% 2L) != 1L ||
			dim(HEC_MI1)[3L] != N)
			stop("HEC_MI1 must have dimensions 20 x (2*windowSize + 1) x ", N, ".")
		if (sum(dimnames(HEC_MI1)[[3L]] != nms) > 0L)
			stop("HEC_MI1 states must be in the same order as background.")
	}
	if (is.null(HEC_MI2)) {
		data("HEC_MI2", envir=environment(), package="DECIPHER")
	} else {
		if (!is.double(HEC_MI2))
			stop("HEC_MI2 must be an array of numerics.")
		if (length(dim(HEC_MI2)) != 5L)
			stop("HEC_MI2 must be a three dimensional array.")
		if (dim(HEC_MI2)[1L] != 20L ||
			dim(HEC_MI2)[2L] != 20L ||
			windowSize[2L] > ((dim(HEC_MI2)[3L] - 1)/2) ||
			dim(HEC_MI2)[3L] != dim(HEC_MI2)[4L] ||
			(dim(HEC_MI2)[3L] %% 2L) != 1L ||
			dim(HEC_MI2)[5L] != N)
			stop("HEC_MI2 must have dimensions 20 x 20 x (2*windowSize + 1) x (2*windowSize + 1) x ", N, ".")
		if (sum(dimnames(HEC_MI2)[[5L]] != nms) > 0L)
			stop("HEC_MI2 states must be in the same order as background.")
	}
	
	states <- .Call("predictHEC",
		myAAStringSet,
		windowSize,
		background,
		HEC_MI1,
		HEC_MI2,
		type,
		paste(nms, collapse=""),
		PACKAGE="DECIPHER")
	
	if (type > 1) {
		states <- lapply(states, function(x) {
			rownames(x) <- nms
			return(x)
		})
	}
	
	return(states)
}
