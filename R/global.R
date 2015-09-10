
co_opt = setGlobalOptions(
	tmpdir = "/icgc/dkfzlsdf/analysis/B080/guz/temp",
	wd = getwd(),
	email = "z.gu@dkfz.de",
	"__flag__" = NULL,
	qsub_enforce = FALSE

)

dir.create(co_opt$tmpdir, showWarnings = FALSE, recursive = TRUE, mode = "0755")
#dir.create(co_opt$wd, showWarnings = FALSE, recursive = TRUE, mode = "0755")

initialize_workspace = function() {
	dir.create(paste0(co_opt$wd, "/.flag"), showWarnings = FALSE)
	co_opt("__flag__" = paste0(co_opt$wd, "/.flag"))

}