use strict;
use warnings;

use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
use File::Spec;
use Getopt::Long;

use lib ".";
use GTFsupport;


# test id:
# ENSG00000160072, ENST00000308647
# has 3utr over many exons

my ($gtf, $type,$ofile);
my $help;
GetOptions(
	"inputfile|gtf=s{1}"=>\$gtf,
	"outfile=s{1}"=>\$ofile,
	"type=s{1}"=>\$type,
	"help"=>\$help
);

if ($help or !$gtf or !-e $gtf or !$type) {
die <<HELP;
--------------------------------------------------
** this script should be run from its storing directory **

index a Stringtie-generated GTF file for each exon.

----- required -----
[-g INPUT_GTF] # full path preferred
 - extension must be *.gtf
[-t TYPE_OF_GTF] # use one of the following: ensembl, stringtie
# ensembl: GTF downloaded from ensembl.org
# stringtie: GTF of de novo assembly by StringTie
# other types may be supported in the future

# ----- optional -----
[-o OUTPUT.hash] # where to save output 
# output file is a perl hash and can be read by Perl's Storable::retrieve,
access the hash by key `_sample_` to view a sample data structure.

--------------------------------------------------

HELP
}

# --- chk gtf ---
printf ">>%s . . .\n", $gtf;
if (!-e $gtf) {
	print "  file doesn't exist!\n";
}
if ($gtf!~/\.gtf$/i) {
	print "  file must have a *.gtf extension!\n";
}

# --- convert to hash based on GTF type---
my $idx;
if ($type eq 'ensembl') {
	$idx=index_ensembl_gtf($gtf);
} elsif ($type eq 'stringtie') {
	$idx=index_stringtie_gtf($gtf);
}

# --- map trx-related positions ---
foreach my $gid (keys %{$idx}) {
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
		delete $idx->{$gid}{$tid}{exon};
		$idx->{$gid}{$tid}{exon} = dclone $exon2;
	}
}

# --- save data ---
if (!$ofile) {
	$ofile=$gtf.'.index.hash';
}
printf "  writing indexed coordinate info to %s . . .\n", $ofile;
nstore($idx, $ofile);


# ---------- subs ------------

sub index_ensembl_gtf {
	my ($gtffile, $testmode)=@_; # must be a stringtie-generated GTF
	open (my $fh, $gtffile);
	# with sample data structure
	my $idx={};
	while (<$fh>) {
		if ($testmode) {
			last if $.==500;
		}
		next if /^#/;
		next if !/\S/;
		chomp;
		my @c=split /\t/;
		my $attr=parse_gtf_attr($c[-1]);
		my $gid=$attr->{gene_id};
		my $tid=$attr->{transcript_id};
		if ($c[2]=~/gene|transcript/) { # save gene meta only if the feature type is 'gene'
			$attr->{start}=$c[3];
			$attr->{end}=$c[4];
			if ($c[2] eq 'gene') {
				$attr->{chr}=$c[0];
				$attr->{strand}=$c[6] eq "+"?1:2;
				$idx->{$attr->{gene_id}}{info}=dclone $attr;
			}
			else {
				$idx->{$gid}{$tid}{info}=dclone $attr;
			}
		} else {
			# cds, exon, start_codon, stop_codon should have "exon_number" key.
			# note that start_codon is within CDS, while stop_codon is right after CDS
			my $inner_idx=-1;
			if ($attr->{exon_number} and $c[2]!~/_codon/) { # use exon_number is only assigned to exon/cds items
				$inner_idx=$attr->{exon_number};
			}
			if (!$idx->{$gid}{$tid}{$c[2]}[0]) { # no length count yet
				$idx->{$gid}{$tid}{$c[2]}[0]=$c[4]-$c[3]+1;
			} else { # has exising length count, get sum
				$idx->{$gid}{$tid}{$c[2]}[0]+=$c[4]-$c[3]+1;
			}
			# add start/end pos
			if ($inner_idx>0) { # for CDS/exon
				$idx->{$gid}{$tid}{$c[2]}[$inner_idx]=[$c[3], $c[4]];
			} else { # the rest feat
				if ((scalar @{$idx->{$gid}{$tid}{$c[2]}})<=1) { # new item, [0] should be occupied
					$idx->{$gid}{$tid}{$c[2]}[1]=[$c[3], $c[4]];
				} else { # existing item
					push @{$idx->{$gid}{$tid}{$c[2]}}, [$c[3], $c[4]];
				}
			}
		}
	}
	# add sample gene ID
	$idx->{_sample_}=$idx->{"ENSG00000160072"};
	return $idx;
}

sub index_stringtie_gtf {
# note: stringtie GTF doesn't treat transcripts on reverse strand in biological order. so this sub also adjust them here.
	my ($gtffile, $testmode)=@_; # must be a stringtie-generated GTF
	open (my $fh, $gtffile);
	# with sample data structure
	my $idx={};
	while (<$fh>) {
		next if /^#/;
		next if !/\S/;
		if ($testmode) {
			last if $.==500;
		}
		chomp;
		my @c=split /\t/;
		my $attr=parse_gtf_attr($c[-1]);
		my $gid=$attr->{gene_id};
		my $tid=$attr->{transcript_id};
		if ($c[2] eq 'transcript') { # first line in this trx, save gene level info
			$idx->{$gid}{info}={
				gene_id=>$gid,
				strand => $c[6] eq '+'? 1:2,
				chr => $c[0],
				start=>$c[3],
				end=>$c[4]
			};
			next;
		}
		if ($attr->{exon_number}==1) { # this row is a new transcript, add trx level info
			$idx->{$gid}{$tid}{info}={
				gene_id=>$gid,
				transcript_id => $tid,
				ref_gene_id => $attr->{ref_gene_id}||undef, # ref gene id
				ref_gene_name => $attr->{ref_gene_name}||undef, # ref gene name
				ref_transcript_id => $attr->{reference_id}||undef, # ref trx id
				start=>$idx->{$gid}{info}{start},
				end=>$idx->{$gid}{info}{end}
			};
			delete $idx->{$gid}{info}{start};
			delete $idx->{$gid}{info}{end};
			$idx->{$gid}{$tid}{exon}[0]=0; # initial length
		}
		# add genomic coordinates for this exon, here, relate the index number in the A-ref
		$idx->{$gid}{$tid}{exon}[$attr->{exon_number}]=[$c[3], $c[4]];
		$idx->{$gid}{$tid}{exon}[0] += $c[4] - $c[3] +1; # total length of all exons
	}
	close ($fh);

	# need to reverse strand==2 exons to biological order
	foreach my $gid (keys %$idx) {
		if ($idx->{$gid}{info}{strand} eq '2') { # on reverse strand. using `eq` but not `==` to avoid conflicts in "sample"
			foreach my $tid (keys %{$idx->{$gid}}) {
				next if $tid eq 'info';
				my $e0=dclone $idx->{$gid}{$tid}{exon};
				my $e1=[$e0->[0]];
				shift @$e0;
				push @$e1, (reverse @$e0);
				$idx->{$gid}{$tid}{exon}=dclone $e1;
			}
		}
	}
	# add ref
	$idx->{_sample_}=$idx->{"STRG.1"};
	return $idx;
}
