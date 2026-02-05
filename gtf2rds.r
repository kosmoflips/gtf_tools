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
		metavar="path/to/output_dir"
	),
	make_option(
		c("--testline"),
		type="integer",
		default=-1,
		help="only read first N lines",
		metavar="N"
	)
))
opt = parse_args(opt_parser)


if (is.null(opt$input_gtf)) {
	stop ('need input gtf [-i]')
}
if (is.null(opt$output)) {
	stop ('need output file path [-o]')
} else {
	odir = opt$output
	if (!dir.exists(odir)) {
		dir.create(odir, recursive = T)
	}
}
# if (is.null(opt$type)) {
# 	stop ('need input gtf type [-t] ensembl or stringtie')
# }



# for saving chr data
odir2 = file.path(odir, 'chr_rds')
if (!dir.exists(odir2)) {
	dir.create(odir2,recursive = T)
}



# ----- load gtf as df -----
chkfile1 = file.path(odir2, 'tmp.data_d2.rds')
if (file.exists(chkfile1)) {
	message('loading existing intermediate GTF data . . .')
	message('  ', chkfile1)
	d2 = readRDS(chkfile1)
} else {
message('loading ', opt$input_gtf, ' . . .')
d0=read.delim(opt$input_gtf, header = F, sep = '\t', comment.char = '#', nrows = opt$testline)
colnames(d0) = c('chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute')

d0 = d0 %>% dplyr::select(-source, -score, -frame)

# i=1

keepattr=c("gene_name","gene_id","gene_version","transcript_id","transcript_version","exon_number","gene_biotype","transcript_biotype")

message('processing data . . .')

d1=list()
# d1b = data.frame()
# for (i in 1:nrow(d0)) {
for (i in 1:nrow(d0)) {
	print (i)
	attr1=d0$attribute[i]
	ax=str_split(attr1, ';\\s*')[[1]]
	nlist=vector()
	vlist=vector()
	for (a1 in ax) {
		if (nchar(a1)==0) {
			next()
		}
		b1=str_split(a1, ' ')[[1]]
		
		if (b1[1] == 'tag') {
			# it's found that many ensembl GTF has info string like "tag CCDS, tag basic ...", so unless they are needed in the future, they're ignored for now.
			next()
		}
		
		# if (b1[1] %in% keepattr) {
			nlist=append(nlist,b1[1])
			vlist=append(vlist,b1[2])
		# }
	}
	
	# names(vlist) = nlist
	# v2=as.list(vlist)
	v2=as.data.frame(t(as.data.frame(vlist)))
	colnames(v2) = nlist
	

	# if ('transcript_id' %in% names(v2)) {
	# 	if (v2$transcript_id != 'ENST00000400921') {
	# 		next()
	# 	}
	# }
	
	d1[[i]] = v2
	# d1b = bind_rows(d1b, v2)
}
d1b = bind_rows(d1)
d0b = d0 %>% select(-attribute)

d2 = bind_cols(d0b,d1b)
gc()
saveRDS(d2, chkfile1)
} # close processing 'd2'

# stop()

# ----- calc rel pos -----
message('calculating relative positions . . .')

allcols1 = colnames(d2)
nestedcols = c('feature','start','end', 'exon_number')
groupedcols=allcols1[! allcols1 %in% nestedcols]

# process by chr to save time & add chkpoint
chrlist0 = unique(d2$chr)
chrlist = chrlist0 [ nchar(chrlist0) <= 6 ] # for chr1, 1, chrMT, MT ... nchar=6 for safe. should work at least for human/mouse

add_exinfo_for_one_chr <- function(initial_gtf, chr_name, saveto, reverse=F) {
	# reverse, get chr NOT in chr_name
	if (isFALSE(reverse)) {
		gtf2 = initial_gtf %>% filter(chr %in% chr_name)
	} else {
		gtf2 = initial_gtf %>% filter(! chr %in% chr_name)
	}
	
	rels0 = gtf2 %>% filter(feature %in% c('gene','transcript'))
	rels1 = gtf2 %>% filter(!feature %in% c('gene','transcript')) %>%
		group_by(chr, strand, gene_id, transcript_id, feature) %>%
		nest()
	# k=1
	for (k in 1:nrow(rels1)) {
		featblock=rels1$data[[k]] %>%
			arrange(start) %>%
			mutate(
				feature1 = NA,
				relative_start = -1,
				relative_end = -1
			)
		strand1=rels1$strand[k]
		if (strand1=='+') {
			posorder=1:nrow(featblock)
		} else {
			posorder=nrow(featblock):1
		}
		p=1
		q=0
		for (x in posorder) {
			len1=featblock$end[x] - featblock$start[x]
			q = p + len1
			featblock$relative_start[x] = p
			featblock$relative_end[x] = q
			p = q + 1
		}
		rels1$data[[k]] = featblock
	}
	
	# merge
	merge = rels1 %>% unnest( cols = c(data) ) %>%
		bind_rows(rels0) %>%
		mutate(feature=ifelse( is.na(feature1), feature, feature1)) %>%
		dplyr::select(-feature1) %>%
		arrange(start)
	
	# save
	wdir = dirname(saveto)
	if (!dir.exists(wdir)) {
		dir.create(wdir, recursive = T)
	}
	saveRDS(merge, saveto)
	
	return (TRUE)
}

# standard chr names
basefname = basename(opt$output)
chrfiles = list()

for (mychr in chrlist) {
	message('  > chr ', mychr, ' . . .')
	ofilechr=file.path(odir2, sprintf('%s.chr%s.rds', basefname,mychr))
	if (!file.exists(ofilechr)) {
		dorev=ifelse(mychr=='OTHER', T, F)
		add_exinfo_for_one_chr(initial_gtf = d2, chr_name = mychr, saveto = ofilechr, reverse = dorev)
		gc()
	} else {
		message('    file exists')
	}
	message('    ofilechr')
	chrfiles[[mychr]] = ofilechr
}


message('  - merging chr data . . .')
masterfinal = data.frame()
for (chrfilex in chrfiles) {
	message('  -- ', chrfilex)
	if (file.exists(chrfilex)) {
		datatmp=readRDS(chrfilex)
		masterfinal = bind_rows(masterfinal, datatmp)
		rm(datatmp)
	} else {
		message('    ! file not found!')
	}
}

for (c1 in c("chr", "feature", "strand", "gene_id", "gene_version", "transcript_id", "transcript_version", "gene_name", "gene_source", "gene_biotype", "transcript_name", "transcript_source", "transcript_biotype", "transcript_support_level", "protein_id", "ccds_id", "exon_number", "exon_version", "protein_version")) {
	masterfinal[[c1]] = factor(masterfinal[[c1]])
}

finalofile=file.path(opt$output,'gtf_all_chr.rds')
message('save to ', finalofile)
saveRDS(masterfinal, finalofile)
