######################################################################################
# this file contains functions that wraps functions with same name in other pacakges
######################################################################################

getChromInfoFromUCSC = function(species) {
	dir = tempdir()

	op = qq.options(READ.ONLY = FALSE)
	qq.options(code.pattern = "@\\{CODE\\}")

	filename = qq("@{species}_getChromInfoFromUCSC")
	if(file.exists(qq("@{dir}/@{filename}"))) {
		df = read.table(qq("@{dir}/@{filename}"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	} else {
		suppressMessages(df <- GenomicFeatures::getChromInfoFromUCSC(species))
		write.table(df, file = qq("@{dir}/@{filename}"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	}
	
	qq.options(op)

	return(df)
}

return2 = function(expr, invisible = FALSE) {
	env = parent.frame()
	
	if(identical(env, .GlobalEnv)) {
		base::return(NULL)
	}
	
	# check on.exit
	if(exists(".on.exit.expression", envir = env)) {
		.on.exit.expression = get(".on.exit.expression", envir = env)
		for(i in seq_along(.on.exit.expression)) {
			eval(.on.exit.expression[[i]], envir = env)
		}
	}

	value = eval(substitute(expr), envir = env)
	
	obj = ls(envir = env)
	
	all_obj = ls(envir = env, all.names = TRUE)
	rm(list = all_obj, envir = env)
	gc(verbose = FALSE)
	
	if(invisible) {
		base::return(invisible(value))
	} else {
		base::return(value)
	}
}

set_counter = function(n) {

	n = as.integer(n)
	i = 1

	f = function() {

		i = as.integer(i)
		pct = sprintf("%.1f", i/n*100)
		cat(paste(rep("\b", 100), collapse=""))
		cat(i, "/", n, " (", pct, "%)", sep = "")

		if(i == n) cat("\n")

		i <<- i + 1
	}
}

diameter = function(x) {
	max(x) - min(x)
}

is.file = function(path) {
	if(length(path) == 1) {
		is.atomic(path) && is.character(path) && file.exists(path)
	} else {
		return(FALSE)
	}
}

sleep = function(time) {
	pb = txtProgressBar(style = 3)
	for(i in seq_len(time)/time) {Sys.sleep(1); setTxtProgressBar(pb, i)}
	close(pb)
}

check_system_command = function(cmd) {
	if(Sys.which(cmd) == "") {
		warning(paste0("Cannot find system command ", cmd, ". Please install it or add the path to PATH.\n"))
	}
}


sort_chr = function(x) {
	y = gsub("^chr(\\d)$", "chr0\\1", x)
	y = gsub("^chr(\\d)_", "chr0\\1_", y)
	x[order(y)]
}

order_chr = function(x) {
	y = gsub("^chr(\\d)$", "chr0\\1", x)
	y = gsub("^chr(\\d)_", "chr0\\1_", y)
	order(y)
}


subset_txdb = function(txdb, chromosome = "chr1") {

	txdump = as.list(txdb)
	txdump$transcripts = txdump$transcripts[txdump$transcripts$tx_chrom %in% chromosome, , drop = FALSE]
	txdump$splicings = txdump$splicings[txdump$splicings$tx_id %in% txdump$transcripts$tx_id, , drop = FALSE]
	txdump$genes = txdump$genes[txdump$genes$tx_id %in% txdump$transcripts$tx_id, , drop = FALSE]
	txdump$chrominfo = txdump$chrominfo[txdump$chrominfo$chrom %in% chromosome, , drop = FALSE]
	txdb2 = do.call(makeTranscriptDb, txdump)

	return(txdb2)
}



findNeighbours = function(gr1, gr2, upstream = 1000, downstream = 1000) {
	gr_extended = gr1
	strd = strand(gr1)
	start(gr_extended) = ifelse(strd == "-", start(gr1) - downstream, start(gr1) - upstream)
	end(gr_extended) = ifelse(strd == "-", end(gr1) + upstream, end(gr1) + downstream)

	mtch = findOverlaps(gr_extended, gr2)

	mtch = as.matrix(mtch)
	neighbours = gr2[mtch[,2]]
	neighbours$distance = distance(gr2[mtch[,2]], gr1[mtch[,1]])

	neighbours$host_id = names(gr1[mtch[,1]])
	if(!is.null(names(gr2))) {
		neighbours = neighbours[names(neighbours) != neighbours$host_id]
	}
	neighbours
}



top_largest_objects = function() {
	envir = parent.frame()
	s = sapply(ls(envir = envir),function(x) eval(parse(text = qq("object.size(@{x})")), envir = envir))
	print(head(sort(s, decreasing = TRUE)))
}


set_proper_seqlengths = function(gr, species) {
	
	chr_len_df = getChromInfoFromUCSC(species)
	
	chr = as.character(chr_len_df[[1]])
	chr_len = chr_len_df[[2]]
	names(chr_len) = chr
	
	slev = seqlevels(gr)
	slev = chr[chr %in% slev]
	seqlevels(gr) = slev
	slen = chr_len[slev]
	seqlengths(gr) = slen
	return(gr)
}


cor_cols = function (x, abs_cutoff = 0.5, size = 1000, mc = 1, ...) {
    
    split_by_block = function(n, size) {
    	size = min(c(n, size))
    	REST <- n%%size
	    LARGE <- n - REST
	    NBLOCKS <- n%/%size
	    GROUP <- rep(1:NBLOCKS, each = size)
	    if (REST > 0)
	        GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
	    split(1:n, GROUP)
    }	

    NCOL <- ncol(x)
    SPLIT = split_by_block(NCOL, size)
    COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
    COMBS <- t(apply(COMBS, 1, sort))
    COMBS <- unique(COMBS)

    nr = nrow(COMBS)

    count_list = mclapply(split_by_block(nr, floor(nr/mc)), function(ind) {

	    count = matrix(0, nrow = NCOL, ncol = length(abs_cutoff))
	    for (i in ind) {
	        COMB <- COMBS[i, ]

	        qqcat("block @{COMB[1]}(row @{(COMB[1]-1)*size+1}~@{COMB[1]*size})/@{max(COMBS[,1])} and @{COMB[2]}(row @{(COMB[2]-1)*size+1}~@{COMB[2]*size})/@{max(COMBS[,2])}\n")
	        G1 <- SPLIT[[COMB[1]]]
	        G2 <- SPLIT[[COMB[2]]]
	        
	        RES <- cor(x[, G1], x[, G2], ...)
	        for(k in seq_along(abs_cutoff)) {
	        	tmp_mat = RES
		        tmp_mat[abs(tmp_mat) > abs_cutoff[k]] = 1
		        tmp_mat[abs(tmp_mat) < abs_cutoff[k]] = 0

		        count[G1, k] = count[G1, k] + rowSums(tmp_mat)
		        if(COMB[1] != COMB[2]) {
		        	count[G2, k] = count[G2, k] + colSums(tmp_mat)
		        }
		    }
	    }
	    return(count)
	}, mc.cores = mc)

    count = matrix(0, nrow = NCOL, ncol = length(abs_cutoff))
    for(i in seq_along(count_list)) {
    	count = count + count_list[[i]]
    }
	dim(count) = c(NCOL, length(abs_cutoff))
	rownames(count) = colnames(x)
	colnames(count) = abs_cutoff
	return(count)
}
