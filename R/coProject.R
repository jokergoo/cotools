
cp = setRefClass("coProject",
	fields = list(
		sample_id = "character",
		type = "character",
		col = "character",
		chromosome = "character",
		species = "character",
		txdb = "ANY",
		methylation_hooks = "ANY",
		WGBS = "ANY",
		expr = "ANY",
		RNASEQ = "ANY",
		genomic_features = "ANY",
		CHIPSEQ = "ANY",
		get_hm = "ANY",
		output = "list"
	)
)

initialize_project = function(sample_id, chromosome, species,
	type, col, output,
	chipseq, 
	txdb, 
	methylation_hooks, expr,
	get_hm,
	genomic_features, dir = ".", gencode_gtf_file = NULL) {

	cp = new("coProject")

	if(missing(sample_id)) {
		stop("'sample_id' is mandatory.")
	}
	cp$sample_id = sample_id

	if(!missing(type)) {
		if(length(type) != length(sample_id)) {
			stop("'type' should have same length as 'sample_id'.")
		}
		cp$type = type
	}

	if(!missing(col)) {
		if(is.numeric(col)) col = rgb(t(col2rgb(col)), maxColorValue = 255)
		if(!missing(type)) {
			if(!is.null(names(col)) && length(col) == length(unique(type))) {
				cp$col = col
			} else if(length(col) == length(sample_id)) {
				cp$col = structure(unique(col), names = unique(type))
			} else {
				stop("'col' should have same length as 'sample_id', or a named vector corresponding to 'type'.")
			}
		} else {
			stop("You should define 'type' if you want to define 'col'.")
		}
		
	}

	if(missing(chromosome)) {
		stop("'chromosome' must be defined.")
	}
	cp$chromosome = as.character(chromosome)

	if(missing(species)) {
		stop("'species' must be defined.")
	}
	cp$species = species

	if(!missing(txdb)) {
		cp$txdb = txdb
	}

	if(!missing(methylation_hooks)) {
		if(!inherits(methylation_hooks, "GlobalOptionsFun")) {
			stop("You should use 'GlobalOptions::setGlobalOptions' to generate 'methylation_hooks'.")
		}
		cat("validate methylation (random pick one chromosome)...\n")
		methylation_hooks$set(sample(chromosome, 1))
		cn = intersect(colnames(methylation_hooks$meth(row_index = 1:2)), sample_id)
		if(length(cn) == 0) {
			stop("samples names for methylation should correspond to 'sample_id'.")
		}
		methylation_hooks$raw(row_index = 1, col_index = 1)
		methylation_hooks$coverage(row_index = 1, col_index = 1)
		methylation_hooks$site(index = 1)
		methylation_hooks$GRanges()

		cp$methylation_hooks = methylation_hooks
		cp$WGBS = sample_id %in% cn
	}

	if(!missing(expr)) {
		if(is.data.frame(expr)) {
			stop("I don't like expression stored in a data frame, change to matrix.")
		}
		cat("validate expression...\n")
		cn = intersect(sample_id, colnames(expr))
		if(length(cn) == 0) {
			stop("column names of 'expr' should correspond to 'sample_id'.")
		}
		expr = expr[, cn, drop = FALSE]
		
		# test row names in 'expr' and names in 'txdb'
		if(!missing(txdb)) {
			genes = genes(txdb)
			gn = intersect(names(genes), rownames(expr))
			if(length(gn) == 0) {
				stop("row names of 'expr' should correspond to 'txdb'.")
			}
		}

		cp$expr = expr
		cp$RNASEQ = sample_id %in% cn
	}

	if(!missing(genomic_features)) {
		if(!all(sapply(genomic_features, inherits, "GRanges"))) {
			stop("'genomic_features' should be a list of 'GRanges' objects.")
		}
		if(is.null(names(genomic_features))) {
			stop("'genomic_features' should be a named list.")
		}
		cp$genomic_features = genomic_features
	}

	if(!missing(chipseq)) {
		cat("validate chipseq (random pick one sample)...\n")
		histome_marks = names(chipseq)
		if(is.null(histome_marks)) {
			stop("'chipseq' should be a list with names of histome marks.")
		}
		lt = lapply(seq_along(chipseq), function(i) {
			cn = intersect(chipseq[[i]], sample_id)
			if(length(cn) == 0) {
				stop(qq("sample ids for '@{histome_marks[i]}' should correspond to 'sample_id'."))
			}
			sample_id %in% cn
		})
		names(lt) = histome_marks
		cp$CHIPSEQ = lt

		if(missing(get_hm)) {
			stop("since you specified 'chipseq', 'get_hm' should be defined as well.")
		}
		# test get_hm
		lapply(seq_along(lt), function(i) {
			gr = get_hm(histome_marks[i], sample_id[lt[[i]]][1])
			if(!inherits(gr, "GRanges")) {
				stop("'get_hm' should return a 'GRanges' object.")
			}
			if(is.null(gr$density)) {
				stop("the 'Granges' object should contain a column 'density'.")
			}
		})
		cp$get_hm = get_hm
	}

	## initialize output folder structure
	dir.create(dir, showWarnings = FALSE, recursive = TRUE)
	dir.create(qq("@{dir}/rds/"), showWarnings = FALSE, recursive = TRUE)
	dir.create(qq("@{dir}/image/"), showWarnings = FALSE, recursive = TRUE)
	dir.create(qq("@{dir}/temp/"), showWarnings = FALSE, recursive = TRUE)
	dir.create(qq("@{dir}/table/"), showWarnings = FALSE, recursive = TRUE)

	if(!missing(output)) {
		cp$output = list(base = dir,
			             rds = qq("@{dir}/rds/"),
			             image = qq("@{dir}/image/"),
			             temp = qq("@{dir}/temp/"),
			             table = qq("@{dir}/table/"))
	}

	return(cp)
}

cp$methods(


show = function() {
		qqcat(
"Summary of the coProject:
| samples: @{length(.self$sample_id)}
| chromosomes: @{length(.self$chromosome)}
| species: @{.self$species}
")
		if(!is.null(.self$type)) {
			qqcat("| types: @{length(unique(.self$type))}\n")
		}
		if(!is.null(.self$txdb)) {
			qqcat("| txdb: @{.self$txdb@.xData$conn@dbname}\n")
		}
		if(!is.null(.self$WGBS)) {
			qqcat("| WGBS: @{sum(.self$WGBS)}\n")
		}
		if(!is.null(.self$RNASEQ)) {
			qqcat("| RNASEQ: @{sum(.self$RNASEQ)}\n")
		}
		if(!is.null(.self$CHIPSEQ)) {
			histome_marks = names(.self$CHIPSEQ)
			qqcat("| CHIPSEQ:\n")
			for(hm in histome_marks) {
				qqcat("|   @{hm}: @{sum(.self$CHIPSEQ[[hm]])}\n")
			}
		}
		if(!is.null(.self$genomic_features)) {
			qqcat("| genomic features:\n")
			gf = names(.self$genomic_features)
			for(i in seq_along(.self$genomic_features)) {
				qqcat("|   @{gf[i]}: @{length(.self$genomic_features[[i]])}\n")
			}
		}
	if(length(.self$output) > 0) {
		qqcat("| output: @{.self$output$base}\n")
	}

	qqcat("\nYou can save this object to get rid of spending a lot of time to re-create it.\n\n")
},

correlated_regions = function(sample_id, chr, ...) {

	l = rep(TRUE, length(.self$sample_id))
	if(!missing(sample_id)) {
		l = .self$sample_id %in% sample_id
	}
	l = l & .self$WGBS & .self$RNASEQ

	sample_id = .self$sample_id[l]
	cat("samples that both have WGBS and RNASEQ:", sum(l), "\n")

	if(is.null(.self$type)) {
		factor = NULL
	} else {
		factor = .self$type[l]
	}

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$correlated_regions(sample_id, 
		expr = .self$expr[, sample_id, drop = FALSE], 
		txdb = .self$txdb, 
		chr = chr,
		factor = factor,
		...)
	
},

filter_correlated_regions = function(template, chromosome = .self$chromosome, ...) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$filter_correlated_regions(chromosome = chromosome, template = template, ...)
},

reduce_cr = function(cr, ...) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$reduce_cr(cr = cr, expr = .self$expr, txdb = .self$txdb, ...)
},

cr_qc = function(template, chromosome = .self$chromosome) {
	
	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$cr_qc(chromosome = chromosome, template = template)	
},

cr_correlated_to_genomic_features = function(cr) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$cr_correlated_to_genomic_features(cr, gf_list = .self$genomic_features, species = .self$species)
},

cr_hilbert = function(cr, template, chromosome = .self$chromosome) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$cr_hilbert(cr, template, txdb = .self$txdb, chromosome = chromosome)
},

compare_meth = function(chr, start, end, x = NULL, x2 = NULL) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$compare_meth(chr, start, end, x, x2)
},

cr_gviz = function(cr, gi, gf_list = NULL, ...) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$cr_gviz(cr, gi, expr = .self$expr, txdb = .self$txdb, 
		species = .self$species, gf_list = .self$genomic_features[gf_list])
},

cr_scatterplot_me = function(cr, gi = NULL, ...) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$cr_scatterplot_me(cr, expr = .self$expr, 
		gi = gi, annotation_color = .self$col)
},

cr_enriched_at_tss = function(cr) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$cr_enriched_at_tss(cr, .self$txdb)
},

enriched_heatmap_list_on_gene = function(cr, cgi, hm_name = NULL, on = "tss", by = "gene", ...) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	if(is.null(.self$CHIPSEQ) || is.null(hm_name)) {
		parent.env(parent.env(e))$enriched_heatmap_list_on_gene(cr, cgi, txdb = .self$txdb, expr = .self$expr, on = no, by = by, ...)
	} else {
		sid = .self$sample_id[.self$CHIPSEQ[[hm_name]]]
		hm_list = lapply(sid, function(x) .self$get_hm(hm_name, x))
		names(hm_list) = sid
		parent.env(parent.env(e))$enriched_heatmap_list_on_gene(cr, cgi, txdb = .self$txdb, expr = .self$expr, hm_list = hm_list, hm_name = hm_name, on = no, by = by, ...)
	}
},

enriched_heatmap_list_on_tss_cgi = function(cr, cgi, hm_name = NULL, by = "gene", ...) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	if(is.null(.self$CHIPSEQ) || is.null(hm_name)) {
		parent.env(parent.env(e))$enriched_heatmap_list_on_tss_cgi(cr, cgi, txdb = .self$txdb, expr = .self$expr, by = by, ...)
	} else {
		sid = .self$sample_id[.self$CHIPSEQ[[hm_name]]]
		hm_list = lapply(sid, function(x) .self$get_hm(hm_name, x))
		names(hm_list) = sid
		parent.env(parent.env(e))$enriched_heatmap_list_on_tss_cgi(cr, cgi, txdb = .self$txdb, expr = .self$expr, hm_list = hm_list, hm_name = hm_name, by = by, ...)
	}
},

wgbs_qcplot = function(sample_id, chromosome = .self$chromosome) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$wgbs_qcplot(sample_id, chromosome)
},

plot_coverage_and_methylation_on_genome = function(sid, chromosome = .self$chromosome,	nw = 10000, ...) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$plot_coverage_and_methylation_on_genome(sid, chromosome, .self$species, nw, ...)

},

plot_multiple_samples_methylation_on_genome = function(sample_id, chromosome = .self$chromosome, nw = 1000, ...) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$plot_multiple_samples_methylation_on_genome(sample_id, annotation = .self$type, 
		annotation_color = chromosome, species = .self$species, nw = nw, ...)

},

global_methylation_distribution = function(sample_id, chromosome = .self$chromosome, by_chr = FALSE) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$global_methylation_distribution(sample_id, annotation = .self$type, 
		annotation_color = .self$col, chromosome = chromosome, by_chr = by_chr)

},

get_mean_methylation_in_genomic_features = function(sample_id, chromosome = .self$chromosome, gf_list = NULL) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$get_mean_methylation_in_genomic_features(sample_id, gf_list = .self$genomic_features[gf_list], 
		chromosome = chromosome)

},

heatmap_diff_methylation_in_genomic_features = function(gr, ...) {

	methylation_hooks = .self$methylation_hooks
	e = environment()

	parent.env(parent.env(e))$heatmap_diff_methylation_in_genomic_features(gr, annotation = .self$type, 
		annotation_color = .self$col, txdb = .self$txdb, gf_list = .self$genomic_features, ...)

},

pipeline = function() {

}

)

