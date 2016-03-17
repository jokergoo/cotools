
# == title 
# extract field from gencode GTF file
#
# == param
# -file path of the gencode GTF file
# -level level 
# -primary_key primary key
# -field field
#  
extract_field_from_gencode = function(file, level = "gene", primary_key = "gene_id", field = "gene_name") {
	df = read.table(pipe(qq("perl /home/guz/project/development/cotools/script/extract_field_from_gencode.pl @{file} @{level} @{primary_key} @{field}")), 
		stringsAsFactors = FALSE)
	return(structure(df[[2]], names = df[[1]]))
}
