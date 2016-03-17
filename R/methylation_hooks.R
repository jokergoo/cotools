

# == title
# Hook functions to extract methylation
#
# == param
# -... Arguments for the parameters, see "details" section
# -RESET reset to default values
# -READ.ONLY whether only return read-only options
# -LOCAL switch local mode
#
# == detail
# there are following hooks:
#
# -set set a chromosome. The function accepts a single chromosome name and 
#      returns an object which is used as the first argument in other functions
# -meth how to extract methylation value. The function should have three arguments,
#       the object returned from ``set()``, index of rows and index of columns
# -raw how to extract row methylation value
# -site how to extract site
# -coverage how to extract coverage
# -GRanges howt to extract sites as a GRanges object
#
methylation_hooks = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}
methylation_hooks = setGlobalOptions(
	set = list(.value = function(chr) stop("you need to define `set`"),
	                             .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 2
	                             	}),
	meth  = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `meth`"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	raw  = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `raw`"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	site         = list(.value = function(obj)  stop("you need to define `site`"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 3
	                             	}),
	coverage     = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `coverage`"),
		                       .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	GRanges      = list(.value =function(obj)  stop("you need to define `GRanges`"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 2
	                             	}),
	obj = NULL
)

.obj_is_set = function() {
	!is.null(methylation_hooks$obj)
}
