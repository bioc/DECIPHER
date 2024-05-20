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
	if (length(background) < 3L || length(background) %% 3L != 0L)
		stop("The length of background must be evenly divisible by 3.")
	if (is.null(HEC_MI1)) {
		data("HEC_MI1", envir=environment(), package="DECIPHER")
	} else {
		if (!is.double(HEC_MI1))
			stop("HEC_MI1 must be an array of numerics.")
		if (length(dim(HEC_MI1)) != 3)
			stop("HEC_MI1 must be a three dimensional array.")
		if (dim(HEC_MI1)[1] != 20 ||
			windowSize[1L] > ((dim(HEC_MI1)[2] - 1)/2) ||
			(dim(HEC_MI1)[2] %% 2) != 1 ||
			dim(HEC_MI1)[3] != 3)
			stop("HEC_MI1 must have dimensions 20 x (2*windowSize + 1) x 3.")
	}
	if (is.null(HEC_MI2)) {
		data("HEC_MI2", envir=environment(), package="DECIPHER")
	} else {
		if (!is.double(HEC_MI2))
			stop("HEC_MI2 must be an array of numerics.")
		if (length(dim(HEC_MI2)) != 5)
			stop("HEC_MI2 must be a three dimensional array.")
		if (dim(HEC_MI2)[1] != 20 ||
			dim(HEC_MI2)[2] != 20 ||
			windowSize[2L] > ((dim(HEC_MI2)[3] - 1)/2) ||
			dim(HEC_MI2)[3] != dim(HEC_MI2)[4] ||
			(dim(HEC_MI2)[3] %% 2) != 1 ||
			dim(HEC_MI2)[5] != 3)
			stop("HEC_MI2 must have dimensions 20 x 20 x (2*windowSize + 1) x (2*windowSize + 1) x 3.")
	}
	
	states <- .Call("predictHEC",
		myAAStringSet,
		windowSize,
		background,
		HEC_MI1,
		HEC_MI2,
		type)
	
	if (type > 1) {
		states <- lapply(states, function(x) {
			rownames(x) <- c("H", "E", "C")
			return(x)
		})
	}
	
	return(states)
}
