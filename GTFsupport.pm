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
map_pos_in_feat
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
	# my $flen=$exon_ref->[0];
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



# find what elem a genomic position belongs to
sub map_pos_in_feat {
# note. because ensembl gtf annotates stop codon individually, this script will combine CDS and stop_codon as one feature.
	my ($anno, $pos_gen, $strand)=@_;
	my $pinfo;
	foreach my $elem (keys %{$anno}) {
		next if $elem eq 'info';
		next if !$anno->{$elem};
		# die Dumper $anno->{five_prime_utr} if $anno->{five_prime_utr};
		# my $pinfo1=map_pos_in_feat_one_elem($anno, $elem, $pos_gen, $strand);
		my $pinfo1=convert_gen2trx($anno->{$elem}, $pos_gen, $strand);
		if ($pinfo1) {
			$pinfo->{$elem}=[$anno->{$elem}[0], $pinfo1];
		}
	}
	return $pinfo;
}



# RETIRE
sub map_pos_in_feat_one_elem { # essentially the same as gen2trx, but combine stop codon into CDS
	my ($anno, $elem, $pos_gen, $strand, $merge_stop_to_cds)=@_;
	my $pinfo=[undef,undef,undef];
	if (!$anno->{$elem} or !$pos_gen or $pos_gen<=0) {
		return $pinfo;
	}
	my $curr_elem_len=0; # total length of this element, from block-1 to block-current
	foreach my $i (1.. ( (scalar @{$anno->{$elem}})-1 )  ) {
		next if !$anno->{$elem}[$i]; # CDS may start from a later exon
		if ($anno->{$elem}[$i][0]<=$pos_gen and $pos_gen<=$anno->{$elem}[$i][1]) { # in this block
			my $cds_len=0;
			# combine stop-codon into CDS
			if ($elem eq 'stop_codon' or $elem eq 'CDS') {
				$pinfo->[0]='CDS';
				$pinfo->[1]=$anno->{CDS}[0]+3;
				if ($elem eq 'stop_codon') {
					$cds_len=$anno->{CDS}[0];
				}
			} else {
				$pinfo->[0]=$elem;
				$pinfo->[1]=$anno->{$elem}[0];
			}

			# calc in-elem pos
			if ($strand eq '2') {
				# q-x+1
				$pinfo->[2]=$anno->{$elem}[$i][1]-$pos_gen+1+$cds_len+$curr_elem_len;
			} else {
				# x-p+1
				$pinfo->[2]=$pos_gen-$anno->{$elem}[$i][0]+1+$cds_len+$curr_elem_len;
			}
			# die Dumper $pinfo, $pos_gen, $anno->{$elem};
			# $ct++;
			last;
		}
		$curr_elem_len+=$anno->{$elem}[$i][1] - $anno->{$elem}[$i][0] +1; # do this here
	}
	return $pinfo;
}

1;
