#!/usr/bin/perl -w
use strict;


my $i=shift;
my $window=15;

if ($#ARGV == 0)
{
        $window=shift;
}
#print STDERR "$window\n";

chomp $i;	
my $s1=`basename $i`;
chomp $s1;
$s1=~s/\.vector_aaindex//g;

#NAME OF FILE (HERE PDB CODE)
my $name=$s1;

#OPEN AAINDEX
open(F1,"$i") or die "Cannot open file \"$i\" : $!\n";	
my @tab_file1=<F1>;
close F1;


######## WITHOUT GAP ###########
#OPEN AAMTX
#CHECK
if (not -e "pdb_vector_aamtx/$name.vector_aamtx")
{
	print STDERR "File \"pdb_vector_aamtx/$name.vector_aamtx\" not exist!\n";
	exit 1;
}

open(F2,"pdb_vector_aamtx/$name.vector_aamtx") or die "Cannot open file \"pdb_vector_aamtx/$name.vector_aamtx\" : $!\n";	
my @tab_file2=<F2>;
close F2;


#CHECK SIZE
if ($#tab_file1 != $#tab_file2)
{
	print STDERR "Not the same file size : \"$i\" and \"pdb_vector_aamtx/$name.vector_aamtx\" ( $#tab_file1 != $#tab_file2 )!\n";
	exit 1;
}

#CREATE MERGE VECOTR
open(F3,">merge_vector_aaindex_aamtx/$name.merge_vector") or die "Cannot create file \"merge_vector_aaindex_aamtx/$name.merge_vector\" : $!\n";

#FOREACH POSITION
for (my $x=0 ; $x <= $#tab_file1 ; $x++)
{
	chomp $tab_file2[$x];
	my @tab_vector_aam=split(/\s+/,$tab_file2[$x]);
	chomp $tab_file1[$x];
	my @tab_vector_aaindex=split(/\s+/,$tab_file1[$x]);


	if ($#tab_vector_aam != (($window)*21)-1)
	{
		my $tmp=(($window)*21)-1;
		print STDERR "Not the correct number of value for aam (no gaps) : $#tab_vector_aam != $tmp \n";
		exit 1;
	}
	
	if ($#tab_vector_aaindex != (($window)*59)-1)
	{
		my $tmp=(($window)*59)-1;
		print STDERR "Not the correct number of value for aaindex : $#tab_vector_aaindex != $tmp \n";
		exit 1;
	}




	my $vector_out="";
        for(my $i=0 ; $i <=  $window-1 ; $i++)
        {
                        my $start_aam=$i*21;
                        my $end_aam=$start_aam+21-1;

			
                        my $start_aaindex=$i*59;
                        my $end_aaindex= $start_aaindex+59-1;

			#print "$i/$window SIZE:$#tab_vector_aaindex $start_aaindex : $end_aaindex \n";
                        #print "\n\$start_aam:$start_aam \$end_aam:$end_aam | \$start_aaindex:$start_aaindex \$end_aaindex:$end_aaindex\n";
                        $vector_out.=sprintf join(' ',@tab_vector_aam[$start_aam..$end_aam]);
                        $vector_out.=" ";
                        $vector_out.=sprintf join(' ',@tab_vector_aaindex[$start_aaindex..$end_aaindex]);
                        $vector_out.=" ";

        }
        chop $vector_out;
        print F3 "$vector_out\n";
}
close F3;	

####### WITH GAP #######
if (not -e "pdb_vector_aamtx/$name.vector_aamtx_gaps")
{
	print STDERR "File \"pdb_vector_aamtx/$name.vector_aamtx_gaps\" not exist!\n";
	exit 1;
#next FOREACH13;
}


open(F2,"pdb_vector_aamtx/$name.vector_aamtx_gaps") or die "Cannot open file \"pdb_vector_aamtx/$name.vector_aamtx_gaps\" : $!\n";	
@tab_file2=<F2>;
close F2;

if ($#tab_file1 != $#tab_file2)
{
	print STDERR "Not the same file size : \"$i\" and \"pdb_vector_aamtx/$name.vector_aamtx_gaps\" ( $#tab_file1 != $#tab_file2 ) !\n";
	exit 1;
#next FOREACH13;
}

#$correct++;
#print " correct=$correct\n";
open(F3,">merge_vector_aaindex_aamtx/$name.merge_vector_gaps") or die "Cannot create file \"merge_vector_aaindex_aamtx/$name.merge_vector_gaps\" : $!\n";
for (my $x=0 ; $x <= $#tab_file1 ; $x++)
{
	chomp $tab_file2[$x];
	my @tab_vector_aam=split(/\s+/,$tab_file2[$x]);
	chomp $tab_file1[$x];
	my @tab_vector_aaindex=split(/\s+/,$tab_file1[$x]);
	
	if ($#tab_vector_aam != (($window)*22)-1)
	{
		my $tmp=(($window)*22)-1;
		print STDERR "Not the correct number of value for aam (gaps) : $#tab_vector_aam != $tmp \n";
		exit 1;
	}
	
	if ($#tab_vector_aaindex != (($window)*59)-1)
	{
		my $tmp=(($window)*59)-1;
		print STDERR "Not the correct number of value for aaindex : $#tab_vector_aaindex != $tmp \n";
		exit 1;
	}



	my $vector_out="";
        for(my $i=0 ; $i <=  $window-1 ; $i++)
        {
                        my $start_aam=$i*22;
                        my $end_aam=$start_aam+22-1;
                        my $start_aaindex=$i*59;
                        my $end_aaindex= $start_aaindex+59-1;

                        #print "\n\$start_aam:$start_aam \$end_aam:$end_aam | \$start_aaindex:$start_aaindex \$end_aaindex:$end_aaindex\n";
                        $vector_out.=sprintf join(' ',@tab_vector_aam[$start_aam..$end_aam]);
                        $vector_out.=" ";
                        $vector_out.=sprintf join(' ',@tab_vector_aaindex[$start_aaindex..$end_aaindex]);
                        $vector_out.=" ";

        }
        chop $vector_out;
        print F3 "$vector_out\n";













}
close F3;


