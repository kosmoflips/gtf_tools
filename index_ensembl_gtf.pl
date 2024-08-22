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

# index standard ensembl gtf file
# output : a perl hash file


# test id:
# ENSG00000160072, ENST00000308647
# has 3utr over many exons

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
	my $idx=ensembl_gtf_indexer($file);
	# map in-exon/in-cDNA positions per transcript
	foreach my $gid (keys %{$idx}) {
		# my $strand=$idx->{$gid}{info}{strand};
		# next if $strand!=2;
		next if $gid eq '_sample_';
		foreach my $tid (keys %{$idx->{$gid}}) {
			next if $tid eq 'info';
			my $e_start=1;
			my $exon2=[$idx->{$gid}{$tid}{exon}[0]];
			foreach my $i (1..(scalar @{$idx->{$gid}{$tid}{exon}}-1)) { # [0] is total exon len
				# plus/minus strands are the same calculation. because input exon A-ref is ordered by exon number
				my $ex=$idx->{$gid}{$tid}{exon}[$i];
				my $p=$e_start;
				my $q=$ex->[1] - $ex->[0] + $p;
				$e_start=$q+1;
				push @$exon2, [ $ex->[0], $ex->[1], $p, $q];
			}
			# print Dumper $idx->{$gid}{$tid}{exon};
			delete $idx->{$gid}{$tid}{exon};
			$idx->{$gid}{$tid}{exon} = dclone $exon2;
			# die Dumper $idx->{$gid}{$tid}{exon};
		}
	}
	my $ofile=$file.'.index.hash';
	printf "  writing indexed coordinate info to %s . . .\n", $ofile;
	nstore($idx, $ofile);
}

print "\n\nall done.";



sub ensembl_gtf_indexer {
	my ($gtffile)=@_; # must be a stringtie-generated GTF
	open (my $fh, $gtffile);
	# with sample data structure
	my $idx={
		"_sample_" => {
			'info' => "{ Hash ref for this gene as in GTF, Plus strand, start, end }",
			'transcript' => {
				"info" => "{ Hash ref for this transcript as in GTF, Plus start, end }",
				"exon" => [ "total_length", ['start-exon1', 'end-exon1', 'cdna-start', 'cdna-end'], ['start-exon2', 'end-exon2', 'cdna-start2', 'cdna-end2'] ],
				"CDS" => [ "total_length", ['start-exon1', 'end-exon1'], ['start-exon2', 'end-exon2'] ],
				"start_codon" => [ "total_length", ['start', 'end'] ],
				"stop_codon" => [ "total_length", ['start', 'end'] ],
				"five_prime_utr" => [ "total_length", ['start', 'end'] ],
				"three_prime_utr" => [ "total_length", ['start', 'end'] ],
				"Selenocysteine" => "[untested. buggy]",
			}
		}
	};
	while (<$fh>) {
		# last if $.==2000;
		next if /^#/;
		next if !/\S/;
		chomp;
		my @c=split /\t/;
		my $attr=parse_gtf_attr($c[-1]);
		# $attr->{gene_id} is the key link
		# currently known types in both human/mouse are:
		# gene   transcript   CDS   exon   five_prime_utr   start_codon   stop_codon   three_prime_utr   Selenocysteine; unsure about "Selenocysteine" but for now, treat it as another feature at same level as exon, start_codon, etc.
		if ($c[2]=~/gene|transcript/) { # save gene meta only if the feature type is 'gene'
			$attr->{start}=$c[3];
			$attr->{end}=$c[4];
			if ($c[2] eq 'gene') {
				$attr->{chr}=$c[0];
				$attr->{strand}=$c[6] eq "+"?1:2;
				$idx->{$attr->{gene_id}}{info}=dclone $attr;
				# die Dumper \@c;
				# die Dumper $attr;
			}
			else {
				$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{info}=dclone $attr;
			}
		} else {
			# cds, exon, start_codon, stop_codon should have "exon_number" key.
			# note that start_codon is within CDS, while stop_codon is right after CDS
			my $inner_idx=-1;
			if ($attr->{exon_number} and $c[2]!~/_codon/) { # use exon_number is only assigned to exon/cds items
				$inner_idx=$attr->{exon_number};
			}
			if (!$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{$c[2]}[0]) { # no length count yet
				$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{$c[2]}[0]=$c[4]-$c[3]+1;
			} else { # has exising length count, get sum
				$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{$c[2]}[0]+=$c[4]-$c[3]+1;
			}
			# add start/end pos
			if ($inner_idx>0) { # for CDS/exon
				$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{$c[2]}[$inner_idx]=[$c[3], $c[4]];
			} else { # the rest feat
				if ((scalar @{$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{$c[2]}})<=1) { # new item, [0] should be occupied
					$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{$c[2]}[1]=[$c[3], $c[4]];
				} else { # existing item
					push @{$idx->{$attr->{gene_id}}{$attr->{transcript_id}}{$c[2]}}, [$c[3], $c[4]];
				}
			}
		}
		
		# if ($c[2]=~/three/) {
			# print Dumper $idx->{$attr->{gene_id}}{$attr->{transcript_id}};<>;
		# }
	}
	return $idx;
}
