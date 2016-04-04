
# == title 
# Extract field from gencode GTF file
#
# == param
# -file the input gtf file
# -level level of the annotation (e.g. gene, transcript, exon, ...)
# -primary_key primary field
# -field field to be retrieved
#
# == details
# Although gtf file can be imported by `GenomicFeatures::makeTranscriptDbFromGFF`, some information
# in the original gtf file will not be imported. This function aims to extract additionally information
# from gtf file.
#
# The function calls external perl script, so you need to have perl installed.
#
# == value
# A vector in which 'primary_key' corresponds to the name and 'field' corresponds to the value
#
extract_field_from_gencode = function(file, level = "gene", primary_key = "gene_id", field = "gene_name") {
	df = read.table(pipe(qq("perl \"@{system.file(package = 'cotools')}/perl_scripts/extract_field_from_gencode.pl\" @{file} @{level} @{primary_key} @{field}")), 
		stringsAsFactors = FALSE)
	return(structure(df[[2]], names = df[[1]]))
}

# == title
# Returns all supported fields in gtf data
#
# == param
# -file the input gtf file
# -level level of the annotation (e.g. gene, transcript, exon, ...)
#
# == details
# These fields are stored in the 9th column in the gtf file.
#
# == value
# A vector of available fields
#
available_gencode_field = function(file, level = "gene") {
	if(grepl("\\.gz$", file)) {
		line = read.table(pipe(qq("zcat @{file} | awk '$3==\"@{level}\"' | head -n 1")), sep = "\t", stringsAsFactors = FALSE)
	} else {
		line = read.table(pipe(qq("awk '$3==\"@{level}\"' @{file} | head -n 1")), sep = "\t", stringsAsFactors = FALSE)
	}
	pair = strsplit(line[1, 9], ";")[[1]]
	gsub("^\\s*(\\w+)\\s.*$", "\\1", pair)
}
