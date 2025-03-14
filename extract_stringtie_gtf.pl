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
my ($infile);
GetOptions(
	"input=s{1}"=>\$infile,
	# "output=s{1}"=>\$ofile,
	"help"=>\$help,
);


if (!$infile or !-e $infile) {
	print "need input stringtie GTF [-i]\n";
	$help=1;
}
# if (!$ofile) {
	# print "need output parsed file save path [-o]\n";
	# $help=1;
# }


if ($help) {die <<USAGE;

-----------------------------------------
# extract stringtie GTF attribution info (FPKM, TPM, etc) to txt

-- input --
[-i stringtie_output.gtf]
-----------------------------------------

USAGE
}

my $ofile=$infile.'.txt';
open (my $fh2, ">", $ofile) or die "can't open output file, maybe the containing folder doesn't exist\n";

open (my $fh, $infile);
my $line1=0;

my @header=qw/gene_id transcript_id cov FPKM TPM /;
while (<$fh>) {
	chomp;
	next if /^#/;
	if (!$line1) {
		printf $fh2 "%s\n", join "\t", @header;
		$line1=1;
		next;
	}

	my @c=split "\t";
	next if $c[2] ne 'transcript';
	# die Dumper \@c;
	my $attr=parse_gtf_attr($c[-1]);
	next if ($attr->{FPKM}+$attr->{TPM}) ==0;
	my @line;
	foreach my $h1 (@header) {
		push @line, $attr->{$h1}||'';
	}
	# die Dumper \@line;
	printf $fh2 "%s\n", join "\t", @line;
	# exit;
}

