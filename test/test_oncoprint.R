files1 = scan(pipe("ls /icgc/dkfzlsdf/analysis/hipo/hipo_047/results_per_pid_pcawg/*/mpileup/*_somatic_functional_snvs_conf_8_to_10.vcf"), what = "character")
samples = gsub("^.*results_per_pid_pcawg/(.*?)/mpileup.*$", "\\1", files1)

snv_oncoprint(files1, samples)

files2 = scan(pipe("ls /icgc/dkfzlsdf/analysis/hipo/hipo_047/results_per_pid_pcawg/*/platypus_indel/*_somatic_functional_indels_conf_8_to_10.vcf"), what = "character")
samples = gsub("^.*results_per_pid_pcawg/(.*?)/platypus_indel.*$", "\\1", files2)

indel_oncoprint(files2, samples)

snv_indel_oncoprint(files1, files2, samples)

