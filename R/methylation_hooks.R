

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
# Methylation from whole genome bisulfite seuqencing is always huge and it does not
# make sense to read them all into the memory. This hook sets how to read the methylation
# data and how to return methylation value (e.g. CpG coverage, methylation rate...)
# 
# There are following hooks:
#
# -set set a chromosome. The function accepts a single chromosome name and 
#      returns an object which is used as the first argument in other functions
# -meth how to extract methylation value. The function should have three arguments:
#       the object returned from ``set()``, index of rows and index of columns. Normally,
#       the first argument (``obj``) can be ignored when calling this hook. Note the methylation
#       matrix should be column names (used as sample id in other functions)
# -raw how to extract raw methylation value, same setting as ``meth``
# -site the function should return a vector of CpG sites
# -coverage how to extract CpG coverage, same setting as ``meth``.
# -GRanges howt to extract CpG sites as a `GenomicRanges::GRanges` object.
#
# Note: positions of CpG sites in a chromosome should be sorted.
#
methylation_hooks = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE) {}
methylation_hooks = setGlobalOptions(
	set = list(.value = function(chr) stop("you need to define `set`"),
	                             .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 2
	                             	}),
	meth  = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `meth` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	raw  = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `raw` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	site         = list(.value = function(obj)  stop("you need to define `site` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 3
	                             	}),
	coverage     = list(.value = function(obj, row_index = NULL, col_index = NULL)  stop("you need to define `coverage` hook"),
		                       .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 4
	                             	}),
	GRanges      = list(.value =function(obj)  stop("you need to define `GRanges` hook"),
		                        .class = "function",
	                             .validate = function(f) {
	                             	length(as.list(args(f))) == 2
	                             	}),
	obj = NULL
)

# .obj_is_set = function() {
# 	!is.null(methylation_hooks$obj)
# }
