use strict;
use warnings;

use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
use File::Spec;
# use File::Path;
# use File::Copy;
# use File::Temp;
use Getopt::Long;

use lib ".";
use GTFsupport;

# convert gtf to bed, can extract selected features

# has confirmed using flag [-cdna] following "bedtools getfasta" resulted correct sequence for both strands.


# note bed format is 0-based. [0-100) means 0-based [0-99], 1-based position [1-100] inclusive

# gtf format (1-based): https://uswest.ensembl.org/info/website/upload/gff.html
# bed format (0-based): https://genome.ucsc.edu/FAQ/FAQformat.html#format1

my @infiles;
my $get_all_feat; # if want to get only selected features (e.g. only get exon and transcript), set this to 0, meanwhile the $get_feat must be configured.
my @feats;
my $cdna;
my $chrname;
my $splitfile;

my $help;
GetOptions(
	"inputfile|gtf=s{1,}"=>\@infiles,
	"allfeat"=>\$get_all_feat,
	"feat|feat=s{1,}"=>\@feats,
	"chr"=>\$chrname,
	"cdna"=>\$cdna,
	"x|split"=>\$splitfile,
	"help"=>\$help
);

if ($help or !@infiles) {
die <<HELP;
--------------------------------------------------
** this script should be run from its storing directory **

extract info from GTF and covert to BED format

----- required -----
[-g INPUT_GTF] # full path preferred

----- optional -----
[-chr] # to write chromosome names as "chr1" instead of "1", default=OFF

--- feature filter ---
* use only one of [-cdna], [-a] or [-f], while [-cdna] > [-f] > [-a] if multiple are given (untested)
[-cdna] # get cdna sequence
  - extract rows with the 3rd column being "exon" and combine by transcript_id.
  - this option is mostly designed for stringtie assembled GTF (other GTF sources untested)
[-a] # get all features as in GTF's 3rd column
[-f feature1 2 ...] # specify one or more features as in GTF's 3rd column
  - e.g. exon, transcript, five_prime_utr, etc.
  - make sure given strings are matched in GTF, case sensitive
[-x] # save bed file for each feature part in a different BED file
  - this flag only works with [-a] or [-f]
  - make sure there are no existing files having the same output file name (e.g. from previous runs on the same source GTF), or the output will be screwed up

--------------------------------------------------

HELP
}

my $get_feat={};
if ($cdna) {
	$get_feat->{exon}=1; # don't consider UTR here, input gtf is assumed to be from assembly program that dont deal with UTR
} elsif (@feats) {
	$get_all_feat=0;
	foreach my $x (@feats) {
		$get_feat->{$x} =1;
	}
} else {
	$get_all_feat=1;
}

foreach my $infile (@infiles) {
	printf "> %s . . .\n", $infile;
	if (!-e $infile or -z $infile) {
		print "  file doesn't exist or is empty, skip!\n";
		next;
	}

	open (my $fh, $infile);
	my $fh2;

	my $ofile;
	if ($cdna) { # [-cdna] is used
		$ofile=$infile.'.bed12';
		open ($fh2, ">", $ofile);
		# print command info
		printf $fh2 "# cdna extracted from file %s\n", $infile;
	}
	elsif (!$splitfile) { # # [-f] or [-a] is used AND save all matching features into one file
	# code to open individual file for each feature is in later file-reading loop
		if (!$get_all_feat) {
			$ofile=$infile.'.sub_feats.bed';
			open ($fh2, ">", $ofile);
			printf $fh2 "# selected features `%s` extracted from file %s\n", (join ",", keys %$get_feat), $infile;
		} else { # all features
			$ofile=$infile.'.bed';
			open ($fh2, ">", $ofile);
			printf $fh2 "# bed converted from file %s\n", $infile;
		}
	}

	my $currtid= {
		trxid=>'',
		bedname=>'',
	}; # use for cdna mode only

	while (<$fh>) {
		next if /^#/;
		chomp;
		my @c=split /\t/;
		if (!$get_all_feat and !$get_feat->{$c[2]}) { # ignore this current line when not getting all features AND current feature isn't required
			next;
		}

		# now the current line feature is required, here it will process individual file opening when to save each feature data to individual bed file.
		if ($splitfile) {
			$ofile=sprintf "%s.%s.bed", $infile,$c[2];
			# make new file or append to existing one
			if (!-e $ofile) {
				open ($fh2, ">", $ofile);
				printf $fh2 "# feature `%s` extracted from file %s\n", $c[2], $infile;
			} else {
				open ($fh2, ">>", $ofile);
			}
		}

		my $attr=parse_gtf_attr($c[-1]);
		# process bedline name
		my $bname=mk_bed_name($attr, $cdna, $get_feat->{exon});
		# print $c[8],"\n",$bname;die;
		if ($cdna) {
			# my ($thisid)=$c[8]=~/transcript_id "(.+?)"/;
			my $thisid=$attr->{transcript_id};
			# die $thisid;
			if ($thisid ne $currtid->{trxid}) { # a new trxid
				if ($currtid->{exon}) { # have currently saved exons
					# print currently saved exon as one cdna
					merge_exons($currtid, $fh2);
				}
				# reset
				$currtid=undef;
				$currtid={
					trxid => $thisid,
					bedname => $bname,
					exon=>[]
				 };
				 # die Dumper $currtid;
			}
			# push this exon
			push @{$currtid->{exon}}, [$c[0], ($c[3]-1), $c[4], $c[6]]; # chrom, start, end, strand
			# print Dumper $currtid; # regardless strand, it's always from small to large coordinates
		}
		else {
			my $chrstring='';
			if ($chrname and $c[0]=~/^(\d+|X|Y|M|MT)$/i) { # add 'chr' before digits and X/Y/M
				$chrstring='chr'.$c[0];
			} else {
				$chrstring=$c[0];
			}
			# die $chrstring;
			printf $fh2 "%s\n", join ("\t",
				$chrstring, # 1. chrom - The name of the chromosome (e.g. chr3, chrY
				$c[3]-1, # 2. chromStart - The starting position. 0 based.
				$c[4], # 3. chromEnd - The ending position. The chromEnd base is NOT included
				$bname, # 4. name - Defines the name of the BED line
				$c[5], # 5. score - A score between 0 and 1000
				$c[6] # 6. strand - Defines the strand. Either "." (=no strand) or "+" or "-".
			)
		}
	}

	# process last exon-block if in 'cdna' mode
	if ($cdna) {
		merge_exons($currtid, $fh2);
	}

	printf "  - write to %s\n", $ofile;
}

# -----------------------------------------------------

sub mk_bed_name {
	# cdna, exon: if extracting features for cdna OR exon
	# cdna mode, will ignore exon number info
	my ($attr, $cdna, $exon)=@_; # attr, a H ref, parsed by &parse_gtf_attr()
# die Dumper $attr;
	my $line2=sprintf 'gene:%s/v%d', $attr->{gene_id}, $attr->{gene_version}||0; # every line should have it
	if ($exon and !$cdna) { # only print exon info when not getting bed12-formatted cdna info
		$line2.=' ';
		if ($attr->{exon_id}) {
			$line2.=sprintf 'exon:%s/v%d', $attr->{exon_id}||'NA', $attr->{exon_version}||0;
		}
		$line2.=sprintf ' exon_number:%d', $attr->{exon_number}||0;
	}
	# stringtie assembly GTF specific names
	if ($attr->{reference_id}) {
		# die $ref_gid;
		$line2.=sprintf ' transcript_reference:%s', $attr->{reference_id};
	}
	if ($attr->{ref_gene_id}) {
		$line2.=sprintf ' gene_reference:%s', $attr->{ref_gene_id};
	}
	if (!$cdna) { # non cDNA mode, not getting specific transcript, add trx name
		$line2.=sprintf ' transcript:%s/v%d', $attr->{transcript_id}||'NA', $attr->{transcript_version}||0;
	}
	else {
		$line2=sprintf '%s/v%d cdna %s ', $attr->{transcript_id}, $attr->{transcript_version}||0, $line2; # add one space at end of string, to separate from bedtools' additional info <<< I don't remember the exact string format design here when relook, will figure out when i encounter formatting bugs in the future probably :3
	}
	return $line2;
}

sub merge_exons {
	my ($data, $fh)=@_; # already adjusted to bed's 0-based position
	# die Dumper $data;
	$fh = *STDOUT if !$fh;
	my $enum=scalar(@{$data->{exon}})-1;
	my $min=$data->{exon}[0][1]; # since it's always from left to right on the positive strand, use 1st exon's start
	my $max=$data->{exon}[$enum][2]; # use last exon's end
	my $line_len='';
	my $line_p='';
	foreach my $e (@{$data->{exon}}) {
		my $len=$e->[2]-$e->[1];
		$line_len .= $len.',';
		$line_p .= ($e->[1] - $min).','; # exon start relative to chr-start, which is $min here
	}
	printf $fh "%s\n", join ("\t",
		$data->{exon}[0][0], # 1. chrom - The name of the chromosome (e.g. chr3, chrY
		$min, # 2. chromStart - The starting position. 0 based.
		$max, # 3. chromEnd - The ending position. The chromEnd base is NOT included
		# $data->{trxid}.' '.$data->{bedname}, # 4. name - Defines the name of the BED line
		$data->{bedname}, # 4. name - Defines the name of the BED line
		0, # 5. score - A score between 0 and 1000
		$data->{exon}[0][3], # 6. strand - Defines the strand. Either "." (=no strand) or "+" or "-".
		0, 0, # 7, 8 thick line. doesn't matter here
		0, # 9. rgb item
		($enum+1), # 10. number of exons
		$line_len, # 11. block length
		$line_p); # 12. block start, relative to chromStart
}
