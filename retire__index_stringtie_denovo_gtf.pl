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


# my $x=retrieve('/home/elly/workcode/shirahagi/tmp/ensembl.hash');
# print Dumper $x->{ENSG00000160072}{info};
# print Dumper $x->{ENSG00000160072}{ENST00000308647};
# die;
# index denovo assembly
# input : denovo assembly generated GTF by stringtie
# output : a perl hash file, listing all start/end coordinates for each row

my @gtfs;
my $help;
GetOptions(
	"inputfile|gtf=s{1,}"=>\@gtfs,
	"help"=>\$help
);

if ($help or !@gtfs) {
die <<HELP;
--------------------------------------------------
** this script should be run from its storing directory **

index a Stringtie-generated GTF file for each exon.

----- required -----
[-g INPUT_GTF] # full path preferred
 - extension must be *.gtf

output file is a perl hash and can be read by Perl's Storable::retrieve,
access the hash by key `_sample_` to view the data structure.
--------------------------------------------------

HELP
}

foreach my $file (@gtfs) {
	printf ">>%s . . .\n", $file;
	if (!-e $file) {
		print "  file doesn't exist, skip!\n";
	}
	if ($file!~/\.gtf$/i) {
		print "  file must have a *.gtf extension, skip!\n";
	}
	my $idx=stringtie_gtf_indexer($file);
	my $ofile=$file.'.index.hash';
	printf "  writing indexed coordinate info to %s . . .\n", $ofile;
	# nstore($idx, $ofile);
}

print "\n\nall done.";



sub stringtie_gtf_indexer {
# IMPORTANT: stringtie GTF doesn't treat transcripts on reverse strand in biological order. so need to change them here.
# step [exon order correction] confirmed with ensembl.org
	my ($gtffile)=@_; # must be a stringtie-generated GTF
	open (my $fh, $gtffile);
	# with sample data structure
	my $idx={};
	while (<$fh>) {
		next if /^#/;
		next if !/\S/;
		chomp;
		my @c=split /\t/;
		# die Dumper \@c;
		my $attr=parse_gtf_attr($c[-1]);
		if ($c[2] eq 'transcript') {
			$idx->{$attr->{gene_id}}{info}={
				gene_id=>$attr->{gene_id},
				strand => $c[6] eq '+'? 1:2,
				chr => $c[0],
				start=>$c[3],
				end=>$c[4]
			};
			next;
		}
		if ($attr->{exon_number}==1) { # this row is a new transcript, add gene/trx info
			$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{info}={
				gene_id=>$attr->{gene_id},
				transcript_id => $attr->{transcript_id},
				ref_gene_id => $attr->{ref_gene_id}||undef, # ref gene id
				ref_gene_name => $attr->{ref_gene_name}||undef, # ref gene name
				ref_transcript_id => $attr->{reference_id}||undef, # ref trx id
				start=>$idx->{$attr->{gene_id}}{info}{start},
				end=>$idx->{$attr->{gene_id}}{info}{end}
			};
			delete $idx->{$attr->{gene_id}}{info}{start};
			delete $idx->{$attr->{gene_id}}{info}{end};
			$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{exon}[0]=0; # initial length
		}
		# add genomic coordinates for this exon, here, relate the index number in the A-ref
		$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{exon}[$attr->{exon_number}]=[$c[3], $c[4]]; # genome-start, genome-end, trx-start, trx-end
		$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{exon}[0] += $c[4] - $c[3] +1;
		# print Dumper $idx->{$attr->{gene_id}}{$attr->{transcript_id}}{exon};<>;
	last if $.==1000;
	}
	#need to reverse strand==2 exons
	foreach my $gid (keys %$idx) {
		if ($idx->{$gid}{info}{strand} eq '2') { # on reverse strand. using `eq` but not `==` to avoid conflicts in "sample"
			foreach my $tid (keys %{$idx->{$gid}}) {
				next if $tid eq 'info';
				# die Dumper $idx->{$gid}{$tid};
				my $e0=dclone $idx->{$gid}{$tid}{exon};
				print Dumper $e0;
				my $e1=[$e0->[0]];
				shift @$e0;
				push @$e1, (reverse @$e0);
				$idx->{$gid}{$tid}=dclone $e1;
				# print Dumper $e1;
				# print Dumper $idx->{$gid}{$tid};<>;
			}
		}
	}
	return $idx;
}
