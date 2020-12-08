#!/usr/bin/perl -w

use strict;


my $i=shift; #file aaindex_vector

my $s1=`basename $i`;
chomp $s1;
$s1=~s/\..*$//g;
my $name=$s1;

# $value can be any regex. be safe
open(F1,"$i") or die "Cannot open file \"$i\" : $!\n";	
my @tab_file1=<F1>;
close F1;

open(F2,"$name.vector_aamtx") or die "Cannot open file \"$name.vector_aamtx\" : $!\n";	
my @tab_file2=<F2>;
close F2;

if ($#tab_file1 != $#tab_file2)
{
    print STDERR "Not the same file size : \"$i\" and \"$name.vector_aamtx\" !\n";
}


open(F3,">$name.merge_vector") or die "Cannot create file \"$name.merge_vector\" : $!\n";
for (my $x=0 ; $x <= $#tab_file1 ; $x++)
{
    chomp $tab_file1[$x];
    chomp $tab_file2[$x];
    print F3 "$tab_file1[$x] $tab_file2[$x]\n";	
}
close F3;	


open(F2,"$name.vector_aamtx_gaps") or die "Cannot open file \"$name.vector_aamtx_gaps\" : $!\n";	
@tab_file2=<F2>;
close F2;

if ($#tab_file1 != $#tab_file2)
{
    print STDERR "Not the same file size : \"$i\" and \"$name.vector_aamtx_gaps\" !\n";
}


open(F3,">$name.merge_vector_gaps") or die "Cannot create file \"$name.merge_vector_gaps\" : $!\n";
for (my $x=0 ; $x <= $#tab_file1 ; $x++)
{
    chomp $tab_file1[$x];
    chomp $tab_file2[$x];
    print F3 "$tab_file2[$x] $tab_file1[$x]\n";	
}
close F3;	


