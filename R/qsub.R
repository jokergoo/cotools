


QSUBENV = new.env()

qsub_reset = function() {
	if(!is.null(QSUBENV$image_file)) {
		file.remove(QSUBENV$image_file, )
	}
	QSUBENV$image_file = NULL
}

qsub = function(name, options, code, share = FALSE, private = NULL, dependency = NULL, 
	tmpdir = co_opt("tmpdir")) {

	if(missing(name)) {
		name = qq("R_chunk_@{paste(sample(c(0:9, letters, LETTERS), 6, replace = TRUE), collapse = '')}")
	}
	if(grepl("\\b-N\\b", options)) {
		stop("job name should be specified by `name` argument.")
	}
	if(!identical(parent.frame(), .GlobalEnv)) {
		stop("qsub can only be called in '.GlobalEnv'")
	}

	flag_file = paste0(co_opt("__flag__"), "/", name, ".success.flag")
	if(!co_opt$qsub_enforce && file.exists(flag_file)) {
		qqcat("Job @{name} is already finished, skip.\n")
		return(NULL)
	}

	suppressWarnings(file.remove(flag_file))

	image_file = tempfile(tmpdir = tmpdir, fileext = ".RData")
	wd = getwd()

	if(share) {
		if(is.null(QSUBENV$image_file)) {
			save.image(file = image_file)
			QSUBENV$image_file = image_file
		}
		image_file = QSUBENV$image_file
	} else {
		save.image(file = image_file)
		QSUBENV$image_file = NULL
	}

	if(!is.null(private)) {
		private_file = tempfile(tmpdir = tmpdir, fileext = ".RData")
		save(private, file = private_file)
	}

code = deparse(substitute(code))
code = paste(code, collapse = "\n")

	# generate a temp R script
	if(is.null(private)) {
	Rscript = qq("
setwd('@{wd}')
load('@{image_file}')

@{code}
")	
	} else {
Rscript = qq("
setwd('@{wd}')
load('@{image_file}')
load('@{private_file}')

@{code}

file.remove('@{private_file}')
")
	}

	if(!share) {
		Rscript = qq("@{Rscript}\n\nfile.remove('@{image_file}')")
	}

	r_file = tempfile(tmpdir = tmpdir, fileext = ".R")
	writeLines(Rscript, con = r_file)

	options = qq("@{options} -N @{name}")

	if(!is.null(dependency)) {
		options = qq("@{options} @{paste(c('-W depend=afterok', dependency), sep = ':', collapse = ':')}")
	}

	sh_file = tempfile(tmpdir = tmpdir, fileext = ".sh")
	writeLines(qq(
"
#!/bin/sh

#PBS -j oe
#PBS -o @{tmpdir}
#PBS -M @{co_opt('email')}

Rscript-3.1.2 @{r_file}

rm @{r_file}

@{if(share) paste('rm', image_file)}
touch @{flag_file}

"), con = sh_file)

	cat("submit job:", name, "\n")

	
	con = pipe(qq("qsub @{options} @{sh_file}"))
	x = readLines(con)
	close(con)
	strsplit(x, "\\t")[[1]][1]
}

