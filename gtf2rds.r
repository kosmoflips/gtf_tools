library(tidyverse)
library(optparse)


opt_parser = OptionParser(option_list=list(
	# make_option(
	# 	c("-t", "--type"),
	# 	type="character",
	# 	default=NULL,
	# 	help="output dir path",
	# 	metavar="ensembl or stringtie"
	# ),
	make_option(
		c("-i", "--input_gtf"),
		type="character",
		default=NULL,
		metavar="path/to/input.gtf"
	),
	make_option(
		c("-o", "--output"),
		type="character",
		default=NULL,
		help="output dir path",
		metavar="path/to/output.gtf.rds"
	)
))
opt = parse_args(opt_parser)


if (is.null(opt$input_gtf)) {
	stop ('need input gtf [-i]')
}
if (is.null(opt$output)) {
	stop ('need output file path [-o]')
} else {
	odir = dirname(opt$output)
	if (!dir.exists(odir)) {
		dir.create(odir, recursive = T)
	}
}
# if (is.null(opt$type)) {
# 	stop ('need input gtf type [-t] ensembl or stringtie')
# }


message('loading ', opt$input_gtf, ' . . .')
d0=read.delim(opt$input_gtf, header = F, sep = '\t', comment.char = '#')
colnames(d0) = c('chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')

d0 = d0 %>% dplyr::select(-source, -score, -frame)

# i=1

keepattr=c("gene_name","gene_id","gene_version","transcript_id","transcript_version","exon_number","gene_biotype","transcript_biotype")

message('  - process data . . .')

d1=list()
for (i in 1:nrow(d0)) {
# for (i in 1:100) {
	attr1=d0$attribute[i]
	ax=str_split(attr1, ';\\s*')[[1]]
	nlist=vector()
	vlist=vector()
	for (a1 in ax) {
		if (nchar(a1)==0) {
			next()
		}
		b1=str_split(a1, ' ')[[1]]
		if (b1[1] %in% keepattr) {
			nlist=append(nlist,b1[1])
			vlist=append(vlist,b1[2])
		}
	}
	names(vlist) = nlist
	v2=as.list(vlist)
	d1[[i]] = v2
}
d1b = bind_rows(d1)
d0b = d0 %>% select(-attribute)

d2 = bind_cols(d0b,d1b)

message('  - save to ', opt$output)
saveRDS(d2, opt$output)
