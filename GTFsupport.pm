package GTFsupport;

# unless otherwise indicated, coordinates dealt here are all 1-based!

use strict;
use warnings;
use Data::Dumper;
use Storable qw/:DEFAULT nstore dclone/;
# use DBI;

use base 'Exporter';
our @EXPORT = qw/
parse_gtf_attr
convert_trx2gen
convert_gen2trx
/;

our @EXPORT_OK = qw//;

our %EXPORT_TAGS = (
# translation=>[qw//]
);


sub parse_gtf_attr {
	# need gene/trx id and exon number
	my ($attrline)=@_; # $needed_attrs, A ref of names to be returned
	#gene_id "ENSG00000284662"; gene_version "1"; transcript_id "ENST00000332831"; transcript_version "4"; exon_number "1"; gene_name "OR4F16"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "OR4F16-201"; transcript_source "ensembl_havana"; transcript_biotype "protein_coding"; tag "CCDS"; ccds_id "CCDS41221"; exon_id "ENSE00002324228"; exon_version "3"; tag "basic"; transcript_support_level "NA (assigned to previous version 3)"
	my @lines=split /\s*;\s*/, $attrline;
	my $attr;
	# convert all attributions to a H ref
	foreach my $fx (@lines) {
		my (@x)=$fx=~/(\S+)\s+"(.+)"/;
		$attr->{$x[0]}=$x[1];
	}
	return $attr;
}


# convert transcript-coordinates to genome
# confirmed with ensembl.org
sub convert_trx2gen {
	my ($exon_ref, $trx_site, $strand)=@_; # exon_ref file must be an A-ref at the $gtf->{ENSG}{ENST}{exon} level, and must be produced by GTFsupport.pm's	"xxx_gtf_indexer"
	# $gtf->{ENSG}{ENST}{exon} = [ total_exon_length, [g1,g2,t1,t2], [g1,g2,t1,t2], [g1,g2,t1,t2], ...]
	$strand=1 if !$strand;
	foreach my $i (1..(scalar(@$exon_ref)-1)) {
		my $list=$exon_ref->[$i];
		my ($t1,$t2)=($list->[2]<$list->[3]?($list->[2],$list->[3]) : ($list->[3],$list->[2]));
		if ($trx_site >= $t1 and $trx_site <=$t2) { # this site in inside the exon, do convertion
			my ($g1,$g2)=($list->[0]<$list->[1]?($list->[0],$list->[1]) : ($list->[1],$list->[0]));
			if ($strand==1) { # forward
				return ($g1+($trx_site-$t1));
			} else { # reverse
				return ($g2-($trx_site-$t1));
			}
		}
	}
	return 0; # error. since all position must be >= 1
}


# convert genome-coordinates to transcript
sub convert_gen2trx {
	my ($exon_ref, $gen_site, $strand)=@_; # exon_ref, same as &convert_trx2gen
	$strand=1 if !$strand;
	foreach my $i (1..(scalar(@$exon_ref)-1)) {
		my $list=$exon_ref->[$i];
		my ($g1,$g2)=($list->[0]<$list->[1]?($list->[0],$list->[1]) : ($list->[1],$list->[0]));
		if ($gen_site >= $g1 and $gen_site <=$g2) { # this site in inside the exon, do convertion
			my ($t1,$t2)=($list->[2]<$list->[3]?($list->[2],$list->[3]) : ($list->[3],$list->[2]));
			if ($strand==1) { # forward
				return ($t1+($gen_site-$g1));
			} else { # reverse
				return ($t1-($gen_site-$g2));
			}
		}
	}
	return 0; # all positions must >= 1
}

1;
