#!/usr/bin/perl
use strict;
use warnings;

use Data::Dumper;
use File::Spec;
use File::Path;
use File::Copy;
use File::Temp;
use Getopt::Long;
use Storable qw/:DEFAULT nstore dclone/;

use lib '.';
use GTFsupport;


my $help;
my ($infile,$ofile);
my ($zero,$gid,$tid);
my ($gtf_ensembl_hash, $gtf_custom_hash, $gtf_type);
GetOptions(
	"input=s{1}"=>\$infile,
	"output=s{1}"=>\$ofile,
	"zero"=>\$zero,
	"gid"=>\$gid,
	"tid"=>\$tid,
	# "gpos"=>\$gpos,
	# "tpos"=>\$tpos,
	"ensembl=s{1}"=>\$gtf_ensembl_hash,
	"gtf=s{1}"=>\$gtf_custom_hash,
	"type=s{1}"=>\$gtf_type,
	"help"=>\$help,
);




if (!$infile or !-e $infile) {
	print "input file table undefined!\n";
	$help=1;
}
if (!$ofile) {
	print "output file save path undefined\n";
	$help=1;
}

if (!$gtf_ensembl_hash or $gtf_ensembl_hash!~/\.hash$/i or !-e $gtf_ensembl_hash) {
	print "Ensembl GTF undefined or isn't a perl hash! Run `index_gtf.pl` first to save GTF to hash!\n";
	$help=1;
}
if ($gtf_custom_hash) { # chk type, for now only supports "stringtie" so no check (yet!)
	if ($gtf_custom_hash!~/\.hash$/i or !-e $gtf_custom_hash) {
		print "Custom GTF undefined or isn't a perl hash! Run `index_gtf.pl` first to save GTF to hash!\n";
		$help=1;
	}
	$gtf_type='stringtie';
	# if (!$gtf_type) {
		# print "custom GTF type undefined!\n";
	# } else {
		# if ($gtf_type!~/stringtie/i) {
			# print "custom GTF type must be one of: stringtie";
		# }
	# }
}
if (($gid and $tid) or (!$gid and !$tid)) {
	print "define either of [-gid] or [-tid]\n";
	$help=1;
}
# if (($gpos and $tpos) or (!$gpos and !$tpos)) {
	# print "define either of [-gpos] or [-tpos]\n";
	# $help=1;
# }

if ($help) {die <<USAGE;

-----------------------------------------
# map given single NT positions to genomic coordinates, and what elements they belong to

-- input --
[-i INPUT.tsv] # Required. One tsv file with the first row being the header, and has two columns: (1) ID, (2) position
[-ensembl ENSEMBL_GTF.hash] # Required. Ensembl GTF path, or a custom GTF that's equivalent to Ensembl. This file is used as the mapping reference.
[-gtf CUSTOM_GTF.hash] # if IDs in ENSEMBL_GTF does NOT correspond to input files, provide a custom GTF that links both input files and ENSEMBL_GTF. One example is StringTie de novo assembly GTF that has its own IDs while linking to Ensembl coordinates/annotations.
[-type] # if CUSTOM_GTF is provided, specify its type. Currently only supports "stringtie".

-- configs --
[-gid OR -tid] # input file ID and position info is linked to gene or transcript level
[-zero] # input file positions are zero based.

-- output --
[-o OUTPUT.txt] # where to save output file


** this script MAY be modified in the future to deal with block of sequences

-----------------------------------------

USAGE
}

open (my $fh2, ">", $ofile) or die "can't open output file, maybe the containing folder doesn't exist\n";

my $gtf_ensembl=retrieve($gtf_ensembl_hash);
my $gtf_custom=retrieve($gtf_custom_hash) if $gtf_custom_hash;

open (my $fh, $infile);
my $curr_id='';
my $ginfo;
my $tinfo;
while (<$fh>) {
	chomp;
	if ($.==1) {
		printf $fh2 "%s\t%s\n", $_, join "\t", qw/gene_id trx_id gene_pos trx_pos chr strand ensembl_gene_id ensembl_trx_id ensembl_gene_name ensembl_feat ensembl_feat_pos ensembl_feat_len ensembl_exon_pos ensembl_exon_len/;
		next;
	}
	my ($id,$pos,@ex)=split /\t/;
	
	# remove ID if needed
	if (!$gtf_custom or $gtf_type!~/stringtie/i) { # input is of custom gtf. stringtie ID should NOT remove "version" like string
		# otherwise remove version
		$id=~s/\.\d+$//;
	}

	# --- find gtf info ---
	# if custom gtf is given, input must correspond to it. otherwise , input must correspond to ensembl GTF
	if ($curr_id ne $id) { # new id, find ref
		$curr_id = $id;
		($ginfo, $tinfo) = find_ref($gtf_custom||$gtf_ensembl, $gid?$id:undef, $tid?$id:undef);
	}

	# --- convert to gpos if tpos ---
	my $pos_gen;
	my $pos_trx;
	if ($tid) { # input is tid, so using trx-related pos
		if ($tinfo->{exon}) {
			$pos_gen=GTFsupport::convert_trx2gen($tinfo->{exon}, $pos, $ginfo->{strand}||'');
		}
		$pos_trx=$pos;
	} else {
		# can't convert to trx related pos unless exact Trx ID is known
		$pos_gen=$pos;
	}

	# --- find what elem it belongs to if annotation available ---
	# if not ensembl, get ensembl ref
	my $anno;
	my $ensembl_gid;
	my $ensembl_tid;
	if ($gtf_custom) {
		$ensembl_gid=$tinfo->{info}{ref_gene_id}||'';
		$ensembl_tid=$tinfo->{info}{ref_transcript_id}||'';
		$anno=$gtf_ensembl->{$ensembl_gid}{$ensembl_tid}||undef;
	} else {
		$anno=$tinfo;
	}

	my $feat = GTFsupport::map_pos_in_feat($anno, $pos_gen, $ginfo->{strand});
	my $exon_pos = GTFsupport::map_pos_in_feat_one_elem($anno, 'exon',  $pos_gen, $ginfo->{strand});

	# --- print output ---
	# die Dumper $ginfo, $tinfo;
	my $printset=[
		$ginfo->{gene_id}, $tinfo->{info}{transcript_id},
		$pos_gen, $pos_trx,
		$ginfo->{chr}, $ginfo->{strand},
		$gtf_custom?$ensembl_gid:'', $gtf_custom?$ensembl_tid:'', $gtf_custom?$tinfo->{info}{ref_gene_name}:$ginfo->{gene_name},
		$feat->[0],$feat->[2],$feat->[1],
		$exon_pos->[2],$exon_pos->[1]
	];
	# print current line content
	print $fh2 $_;
	foreach my $x (@$printset) {
		printf $fh2 "\t%s", $x||'';
	}
	print $fh2 "\n";
}


sub find_ref { # find trx level ref in hash
	my ($gtf, $gid, $tid) =@_;
	my ($info, $tref);
	if (!$tid) { # have gene ID only, only need info
		$info=$gtf->{$gid}{info};
	} else { # have trx ID, need both info and pos ref
		if (!$gid) { # find gid first
			foreach my $g1 (keys %$gtf) {
				foreach my $t1 (keys %{$gtf->{$g1}}) {
					if ($tid eq $t1) {
						$info=$gtf->{$g1}{info};
						$tref=$gtf->{$g1}{$t1};
						last;
					}
				}
			}
		} else {
			$info=$gtf->{$gid}{info};
			$tref=$gtf->{$gid}{$tid};
		}
	}
	return ($info, $tref);
}
