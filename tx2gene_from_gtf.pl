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

# extract gene id/trx id from gtf, output tsv

my ($infile, $ofile,$force);
my $help;
GetOptions(
	"inputfile|gtf=s{1}"=>\$infile,
	"outfile=s{1,}"=>\$ofile,
	"force"=>\$force,
	"help"=>\$help
);

if ($help or !$infile or !-e $infile or -z $infile) {
die <<HELP;
--------------------------------------------------
** extract GTF and write tsv table of trx/gene IDs **

----- required -----
[-g INPUT_GTF] # full path of input GTF file

----- optional -----
[-o OUTPUT_PATH] # where to save output
[-f] # if outfile exists, force to overwrite

--------------------------------------------------

HELP
}

printf "> %s . . .\n", $infile;

if ($ofile) {
	my @sdir=File::Spec->splitpath($ofile);
	pop @sdir;
	my $d1=File::Spec->catpath(@sdir);
	mkpath $d1 if !-d $d1;
} else {
	$ofile=$infile.'.tx2gene_id.txt';
}

if (-e $ofile and !$force) {
	printf "output file %s exists. if need to force overwrite, use flag [-f]\n", $ofile;
	die;
}

open (my $fh2, ">", $ofile);
printf $fh2 "%s\n", join "\t", qw/chr transcript_id transcript_version gene_id gene_version gene_name/;

my $saved;
open (my $fh, $infile);
while (<$fh>) {
	next if /^#/;
	chomp;
	my @c=split /\t/;

	my $attr=parse_gtf_attr($c[-1]);
	next if !$attr->{transcript_id};

	# save trx info if new
	if (!$saved->{$attr->{transcript_id}}) {
		# die Dumper $saved;
		$saved->{$attr->{transcript_id}}=1;
		printf $fh2 "%s\n", join "\t", ( $c[0], $attr->{transcript_id}||'', $attr->{transcript_version}||'', $attr->{gene_id}||'', $attr->{gene_version}||'', $attr->{gene_name}||'' );
	}
}

printf "  - done. file saved to %s.\n", $ofile;
