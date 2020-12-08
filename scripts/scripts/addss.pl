#!/usr/bin/perl
#
# addss.pl
# Add PSIPRED secondary structure prediction (and DSSP annotation) to an MSA or HMMER file.
# Output format is A3M (for input alignments) or HMMER (see User Guide).

#     HHsuite version 3.0.0 (15-03-2015)
#
#     Reference:
#     Remmert M., Biegert A., Hauser A., and Soding J.
#     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
#     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

#     (C) Johannes Soeding and Michael Remmert, 2012

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#     We are very grateful for bug reports! Please contact us at soeding@mpibpc.mpg.de

use lib  "/home/vanille/gelly/PROJECTS/HHBLITS/hh-suite/scripts";
my $hhscripts="/home/vanille/gelly/PROJECTS/HHBLITS/hh-suite/scripts";
#use HHPaths;    # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use Align;      # Needleman-Wunsch and Smith-Waterman alignment functions
use File::Temp qw/ tempfile tempdir /;
use File::Copy;
use strict;

my $ss_cit =
"PSIPRED: Jones DT. (1999) Protein secondary structure prediction based on position-specific scoring matrices. JMB 292:195-202.";

# Module needed for aligning DSSP-sequence

$| = 1;       # Activate autoflushing on STDOUT

# Default values:
our $v = 4;    # verbose mode

my $numres   = 0;     # number of residues per line for secondary structure
my $informat = "a3m"; # input format
my $neff     = 7;     # use alignment with this diversity for PSIPRED prediction
my $program  = $0;    # name of perl script
my $pdbfile;

my $help = "
addss.pl from HHsuite  
Add PSIPRED secondary structure prediction (and DSSP annotation) to a multiple sequence alignment (MSA) 
or HMMER (multi-)model file. 

If the input file is an MSA, the predicted secondary structure and confidence values are added as 
special annotation sequences with names >ss_pred, >ss_conf, and >ss_dssp to the top of the output 
A3M alignment. If no output file is given, the output file will have the same name as the input file, 
except for the extension being replaced by '.a3m'. Allowed input formats are A3M (default), 
A2M/FASTA (-fas, -a2m), CLUSTAL (-clu), STOCKHOLM (-sto), HMMER (-hmm).

If the input file contains HMMER models, records SSPRD and SSCON containing predicted secondary 
structure and confidence values are added to each model. In this case the output file name is 
obligatory and must be different from the input file name.

Usage: perl addss.pl <ali_file> [<outfile>] [-fas|-a3m|-clu|-sto]  
  or   perl addss.pl <hhm_file> <outfile> -hmm  
\n";

# Variable declarations
my $line;
my @seqs;    # sequences from infile (except >aa_ and >ss_pred sequences)
my $query_length;
my $header;    # header of MSA: everything before first '>'
my $name;      # query in fasta format: '>$name [^\n]*\n$qseq\n'
my $qseq;      # residues of query sequence
my $infile;
my $outfile;
my $ss_pred = "";    # psipred ss states
my $ss_conf = "";    # psipred confidence values
my $ss_dssp;         # dssp states as string
my $sa_dssp;    # relative solvent accessibility from dssp as string {A,B,C,D,E} A:absolutely buried, B:buried, E:exposed
my $aa_dssp;    # residues from dssp file as string
my $aa_astr;    # residues from infile as string
my $q_match;    # number of match states in query sequence
my $xseq;       # sequence x returned from Align.pm
my $yseq;       # sequence y returned from Align.pm
my $Sstr;       # match sequence returned from Align.pm

###############################################################################################
# Processing command line input
###############################################################################################

if ( @ARGV < 1 ) { die($help); }

my $options = "";
for ( my $i = 0 ; $i < @ARGV ; $i++ ) { $options .= " $ARGV[$i] "; }

#Input format fasta?
if    ( $options =~ s/ -fas\s/ /g ) { $informat = "fas"; }
elsif ( $options =~ s/ -a2m\s/ /g ) { $informat = "a2m"; }
elsif ( $options =~ s/ -a3m\s/ /g ) { $informat = "a3m"; }
elsif ( $options =~ s/ -clu\s/ /g ) { $informat = "clu"; }
elsif ( $options =~ s/ -sto\s/ /g ) { $informat = "sto"; }
elsif ( $options =~ s/ -hmm\s/ /g ) { $informat = "hmm"; }

if ( $options =~ s/ -v\s+(\d+) / /g ) { $v = $1; }

# Set input and output file
if ( $options =~ s/ -i\s+(\S+) // )   { $infile  = $1; }
if ( $options =~ s/ -o\s+(\S+) // )   { $outfile = $1; }
if ( $options =~ s/^\s*([^-]\S*) // ) { $infile  = $1; }
if ( $options =~ s/^\s*([^-]\S*) // ) { $outfile = $1; }

# Warn if unknown options found or no infile/outfile
if ( $options !~ /^\s*$/ ) {
	$options =~ s/^\s*(.*?)\s*$/$1/g;
	die("Error: unknown options '$options'\n");
}
if ( !$infile ) { print($help); exit(1); }

my $v2 = $v - 1;
if ( $v2 > 2 ) { $v2--; }
if ( $v2 < 0 ) { $v2 = 0; }

if ( $informat eq "hmm" && !$outfile ) {
	print(
"Error: no output file given. With the -hmm option an output file is obligatory\n"
	);
	exit(1);
}

###############################################################################################
# Reformat input alignment to a3m and psiblast-readable format and generate file with query sequence
###############################################################################################

my $inbase;    # $inbasename of infile: remove extension
my $inroot;    # $inbasename of infile: remove path and extension
if ( $infile =~ /(.*)\..*/ ) { $inbase = $1; }
else { $inbase = $infile; }    # remove extension
if ( $inbase =~ /.*\/(.*)/ ) { $inroot = $1; }
else { $inroot = $inbase; }    # remove path

# Create tmpfile
my $tmpdir;
if ( $v <= 3 ) { $tmpdir = tempdir( CLEANUP => 1 ); }
else { $tmpdir = tempdir( CLEANUP => 0 ); }
my ( $tmpf, $tmpfile ) = tempfile( DIR => $tmpdir );
my $tmpfile_no_dir;
if ( $tmpfile =~ /.*\/(.*)/ ) {
	$tmpfile_no_dir = $1;
}
else {
	$tmpfile_no_dir = $tmpfile;
}

my $tmp_outfile = "$tmpfile.out";

if ( $infile eq "stdin" ) {
	my @stdin = <STDIN>;
	$infile = "$tmpfile.stdin";
	open(OUT, ">$infile");
	foreach my $line (@stdin) {
  	print OUT $line;
	}
	close(OUT);
}

############################################################################################

if ( $informat ne "hmm" ) {
	if ( !$outfile ) { $outfile = "$inbase.a3m"; }

	# Use first sequence to define match states and reformat input file to a3m and psi
	if ( $informat ne "a3m" ) {
		system("$hhscripts/reformat.pl -v $v2 -M first $informat a3m $infile $tmpfile.1.in.a3m"
		);
	}
	else {
		system("cp $infile $tmpfile.1.in.a3m");
	}

	# Sanitise the input file - remove all '>' characters, except the one at the start of the string
	# Note that this will corrupt the file if the "header" is not separated by a newline from the
	# first line that starts with '>'
	open( INFILE, "<$tmpfile.1.in.a3m" );
	open( OUTFILE, ">$tmpfile.in.a3m" );
	while ( $line = <INFILE> ) {
		$line =~ s/([^>]+)>/$1/g;
		print OUTFILE $line;
	}
	close(INFILE);
	close(OUTFILE);

	# Read query sequence
	open( INFILE, "<$tmpfile.in.a3m" )
	  or die("ERROR: cannot open $tmpfile.in.a3m!\n");
	$/ = ">";    # set input field separator
	my $i = 0;
	$qseq   = "";
	$header = <INFILE>;
	$header =~ s />$//;
	while ( $line = <INFILE> ) {
		$line =~ s/>$//;
		if ( $line =~ /^ss_/ || $line =~ /^aa_/ ) { next; }
		$seqs[ $i++ ] = ">$line";
		if ( !$qseq ) {
			$line =~ s/^(.*)[^\n]*//;
			$name = $1;
			$qseq = $line;
			$qseq =~ s/\n//g;
		}
	}
	close(INFILE);
	$/ = "\n";    # set input field separator

	if ( $qseq =~ /\-/ ) {

		# First sequence contains gaps => calculate consensus sequence
		system("hhconsensus -i $tmpfile.in.a3m -s $tmpfile.sq -o $tmpfile.in.a3m > /dev/null");

	}
	else {
		$query_length = ( $qseq =~ tr/A-Z/A-Z/ );
		$qseq =~ tr/A-Z//cd;    # remove everything except capital letters

		# Write query sequence file in FASTA format
		open( QFILE, ">$tmpfile.sq" )
		  or die("ERROR: can't open $tmpfile.sq: $!\n");
		printf( QFILE ">%s\n%s\n", $name, $qseq );
		close(QFILE);
	}

	# Filter alignment to diversity $neff
	if ( $v >= 1 ) { printf(STDERR "Filtering alignment to diversity $neff ...\n"); }
	system("hhfilter -v $v2 -neff $neff -i $tmpfile.in.a3m -o $tmpfile.in.a3m");

	# Reformat into PSI-BLAST readable file for jumpstarting
	system("$hhscripts/reformat.pl -v $v2 -r -noss a3m psi $tmpfile.in.a3m $tmpfile.in.psi"
	);

	open( ALIFILE, ">$tmp_outfile" )
	  || die("ERROR: cannot open $tmp_outfile: $!\n");
	printf( ALIFILE "%s", $header );

#	# Add DSSP sequence (if available)
#	if ( $dssp ne "" ) {
#		if ( !&AppendDsspSequences("$tmpfile.sq") ) {
#			if ($numres) {
#
#				# insert a line break every $numres residues
#				$ss_dssp =~ s/(\S{$numres})/$1\n/g;
#			}
#			printf( ALIFILE ">ss_dssp\n%s\n", $ss_dssp );
#			if ( $v >= 1 ) { print(STDERR "\nAdding DSSP state sequence ...\n"); }
#		}
#	}
#
	# Secondary structure prediction with psipred
	if ( $v >= 2 ) {
		print(STDERR "Predicting secondary structure with PSIPRED ... ");
	}
	&RunPsipred("$tmpfile.sq");

	if ( open( PSIPREDFILE, "<$tmpfile.horiz" ) ) {
		$ss_conf = "";
		$ss_pred = "";

		# Read Psipred file
		while ( $line = <PSIPREDFILE> ) {
			if    ( $line =~ /^Conf:\s+(\S+)/ ) { $ss_conf .= $1; }
			elsif ( $line =~ /^Pred:\s+(\S+)/ ) { $ss_pred .= $1; }
		}
		close(PSIPREDFILE);
		$ss_conf =~ tr/0-9/0/c;    # replace all non-numerical symbols with a 0
		if ($numres) {
			$ss_pred =~ s/(\S{$numres})/$1\n/g
			  ;                    # insert a line break every $numres residues
			$ss_conf =~ s/(\S{$numres})/$1\n/g
			  ;                    # insert a line break every $numres residues
		}
		if (length($ss_pred) > 0 && length($ss_conf) > 0) {
			printf( ALIFILE ">ss_pred PSIPRED predicted secondary structure\n%s\n", $ss_pred);
			printf( ALIFILE ">ss_conf PSIPRED confidence values\n%s\n", $ss_conf );
		}
	}

	# Append alignment sequences to psipred sequences
	for ( $i = 0 ; $i < @seqs ; $i++ ) {
		printf( ALIFILE "%s", $seqs[$i] );
	}
	close(ALIFILE);

	if ( !$outfile ) {
		$outfile = "$inbase.a3m";
	}

	if($outfile ne "stdout") {
		move( $tmp_outfile, $outfile );
	}
	else {
		open(IN, "<$tmp_outfile");
		foreach $line(<IN>) {
			print $line;
		}
		close(IN)
	}

	if ( $v >= 2 ) { print(STDERR "done \n"); }
}

if ( $v <= 3 ) {
	unlink("$tmpfile.in.a3m");
	unlink("$tmpfile.in.psi");
	unlink("$tmpfile.horiz");
	unlink("$tmpfile.dssp");
}

exit;

##############################################################################################
# Run SS prediction starting from alignment in $tmpfile.in.psi (called by BuildAlignment)
##############################################################################################
sub RunPsipred() {

	# This is a simple script which will carry out all of the basic steps
	# required to make a PSIPRED V2 prediction. Note that it assumes that the
	# following programs are in the appropriate directories:
	# blastpgp - PSIBLAST executable (from NCBI toolkit)
	# makemat - IMPALA utility (from NCBI toolkit)
	# psipred - PSIPRED V2 program
	# psipass2 - PSIPRED V2 program

	my $infile = $_[0];
	my $basename;    #file name without extension
	my $rootname;    #basename without directory path
	if ( $infile =~ /^(.*)\..*?$/ ) { $basename = $1; }
	else { $basename = $infile; }
	if ( $basename =~ /^.*\/(.*?)$/ ) { $rootname = $1; }
	else { $rootname = $basename; }

	# Does dummy database exist?
	if ( !-e "$dummydb.phr" ) {
		if ( !-e "$dummydb" ) {
			die "Error in addss.pl: Could not find $dummydb\n";
		}

		system("cp $infile $dummydb");
		system("/home/vanille/souza/PROJECTS/NANORION/LAST_ORION/bin/PSI_BLAST/bin/formatdb -i $dummydb");
		if ( !-e "$dummydb.phr" ) {
			die "Error in addss.pl: Could not find nor create index files for $dummydb\n";
		}
	}

# Start Psiblast from checkpoint file tmp.chk that was generated to build the profile
	system(
		"/home/vanille/souza/PROJECTS/NANORION/LAST_ORION/bin/PSI_BLAST/bin/blastpgp -b 1 -j 1 -h 0.001 -d $dummydb -i $infile -B $tmpfile.in.psi -C $tmpfile.chk 1> $tmpfile.blalog 2> $tmpfile.blalog"
	);

	#print("Predicting secondary structure...\n");
	system( "echo " . "$tmpfile_no_dir" . ".chk > $tmpfile.pn\n" );
	system( "echo " . "$tmpfile_no_dir" . ".sq  > $tmpfile.sn\n" );
	system("/home/vanille/souza/PROJECTS/NANORION/LAST_ORION/bin/PSI_BLAST/bin/makemat -P $tmpfile");

#	# Start Psiblast from checkpoint file tmp.chk that was generated to build the profile
#	if ( -e "$datadir/weights.dat4" ) {    # Psipred version < 3.0
#		system(
#			"$execdir/psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 $datadir/weights.dat4 > $tmpfile.ss"
#		);
#	}
#	else {
#		&HHPaths::System(
#			"$execdir/psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $tmpfile.ss"
#		);
#	}
#
#	&HHPaths::System(
#		"$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $tmpfile.ss2 $tmpfile.ss > $tmpfile.horiz"
#	);
#
#	# Remove temporary files
#	if ( $v <= 3 ) {
#		unlink(
#			split ' ',
#				"$tmpfile.pn $tmpfile.sn $tmpfile.mn $tmpfile.chk $tmpfile.blalog $tmpfile.mtx $tmpfile.aux $tmpfile.ss $tmpfile.ss2 $tmpfile.sq"
#		);
#	}
	return;
}


