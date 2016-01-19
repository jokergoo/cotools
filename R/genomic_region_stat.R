########################################################
# this file should contain functions which visualize general 
# distribution of genomic regions 


# == title
# Basic statistics on genomic regions
#
# == param
# -gr_list a list of `GenomicRanges::GRanges`.
# -annotation a vector which contains levels of samples, better have names which correspond to the names of ``gr_list``
# -annotation_color colors corresponding to levels of annotations
# -main title of the plot
# -species species, necessary if ``type`` equals to ``proportion``.
# -type type of statistics
# -by_chr take all chromosomes as a whole or calculate statistics for every chromosome?
#
# == details
# For ``type`` settings:
#
# -proportion proportion of total length of regions compared to the whole genome
# -number number of regions
# -median_width median width of regions
#
# == value
# A data frame which contains statistics for each chromosome in each sample.
#
basic_genomic_regions_stat = function(gr_list, annotation = NULL, annotation_color = NULL, 
	main = NULL, species = "hg19", type = c("proportion", "number", "median_width"),
	by_chr = FALSE) {

	type = match.arg(type)[1]

	if(is.null(names(gr_list))) {
		names(gr_list) = seq_along(gr_list)
	}
	sid = names(gr_list)

	# just need the chromosome information
	chromInfo = read.chromInfo(species = species)
            chr_len = sort(chromInfo$chr.len, decreasing = TRUE)

            # sometimes there are small scaffold
            i = which(chr_len[seq_len(length(chr_len)-1)] / chr_len[seq_len(length(chr_len)-1)+1] > 5)[1]
            if(length(i)) {
                chromosome = chromInfo$chromosome[chromInfo$chromosome %in% names(chr_len[chr_len >= chr_len[i]])]
            } else {
                chromosome = chromInfo$chromosome
            }

            category = chromosome
        
    chromInfo = read.chromInfo(species = species, chromosome.index = category)
	
	if(!is.null(annotation)) {
		if(is.null(names(annotation))) names(annotation) = sid
		annotation = annotation[sid]
	}

	n_gr = length(gr_list)
	if(by_chr) {
		n_chr = length(chromInfo$chromosome)
		if(type %in% "proportion") {
			# proportion in genome
			lt = lapply(chromInfo$chromosome, function(chr) {
					sapply(gr_list, function(gr) {
						gr = gr[ seqnames(gr) == chr ]
						if(length(gr) == 0) {
							return(0)
						} else {
							sum(width(gr))/chromInfo$chr.len[chromInfo$chromosome == chr]
						}
					})
				})
			x = unlist(lt)
			ylab = "proportion in genome"

		} else if(type %in% "number") {
			# number of regions
			# proportion in genome
			lt = lapply(chromInfo$chromosome, function(chr) {
					sapply(gr_list, function(gr) {
						gr = gr[ seqnames(gr) == chr ]
						length(gr)
					})	
				})
			x = unlist(lt)
			ylab = "numbers of regions"

		} else if(type %in% "median_width") {
			# median length
			# proportion in genome
			lt = lapply(chromInfo$chromosome, function(chr) {
					sapply(gr_list, function(gr) {
						gr = gr[ seqnames(gr) == chr ]
						if(length(gr) == 0) {
							return(0)
						} else {
							median(width(gr))
						}
					})	
				})
			x = unlist(lt)
			ylab = "median width of regions"

		}	
	} else {
		if(type %in% "proportion") {
			# proportion in genome
			genome.len = sum(chromInfo$chr.len)
			x = sapply(gr_list, function(gr) sum(width(gr))/genome.len)
			ylab = "proportion in genome"

		} else if(type %in% "number") {
			# number of regions
			x = sapply(gr_list, length)
			ylab = "numbers of regions"

		} else if(type %in% "median_width") {
			# median length
			x = sapply(gr_list, function(gr) median(width(gr)))
			ylab = "median width of regions"
			
		}
	}

	if(by_chr) {   # use lattice plot
		df = data.frame(x = x, 
			sid = rep(sid, times = length(chromInfo$chromosome)), 
			chr = rep(chromInfo$chromosome, each = length(sid)),
			stringsAsFactors = FALSE)
		df$chr = factor(df$chr, levels = chromInfo$chromosome)
		df$sid = factor(df$sid, levels = sid)	
	} else {   # use barplot
		sid = names(x)
		df = data.frame(x = x, sid = factor(sid, levels = sid))
	}

	
	if(is.null(annotation)) {
		gp = ggplot(df, aes(y = x, x = sid)) + geom_bar(stat = "identity")
	} else {
		df = cbind(df, annotation = annotation[df$sid])
		gp = ggplot(df, aes(y = x, x = sid)) + geom_bar(stat = "identity", aes(fill = annotation))
	}

	if(by_chr) gp = gp + facet_wrap( ~ chr)

	if(!is.null(annotation_color)) {	
		gp = gp + scale_fill_manual(values = annotation_color)
	}

	gp = gp + xlab("") + ylab(ylab) + theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + ggtitle(main)
	print(gp)

	return2(df, invisible = TRUE)
}
