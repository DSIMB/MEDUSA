#!/usr/bin/perl -w

use strict;


my $i=shift; #file

my $s1=`basename $i`;
chomp $s1;
$s1=~s/\..*$//g;
my $name=$s1;

# $value can be any regex. be safe

    open(F1,"$name.merge_vector") or die "Cannot open file \"$name.merge_vector\" : $!\n";	
    my @tab_file1=<F1>;
    close F1;

    open(F11,"$name.merge_vector_gaps") or die "Cannot open file \"$name.merge_vector_gaps\" : $!\n";	
    my @tab_file11=<F11>;
    close F11;


    open(F2,"$name.PB.fasta") or die "Cannot open file \"$name.PB.fasta\" : $!\n";	
    my @tab_file2=<F2>;
    close F2;

    my $pb="";
FOREACH2:foreach my $line (@tab_file2)
         {
             next FOREACH2 if ($line =~/^>/);
             chomp $line;
             $pb.=$line;
         }
         my @tab_pb=split('',$pb);


         if ($#tab_file1 != $#tab_pb)
         {
             print STDERR "Not the same vector size : \"$i\":$#tab_file1 and PB \"$name.PB.fasta\":$#tab_file2 !\n";
         }

         open(F3,">$name.class_merge_vector") or die "Cannot create file \"$name.class_merge_vector\" : $!\n";
         for (my $x=0 ; $x <= $#tab_file1 ; $x++)
         {
             chomp $tab_file1[$x];
             chomp $tab_pb[$x];
             print F3 "$tab_pb[$x] $tab_file1[$x]\n";
         }
         close F3;

         if ($#tab_file11 != $#tab_pb)
         {
             print STDERR "Not the same vector size : \"$name.merge_vector_gaps\" and PB \"$name.PB.fasta\" !\n";
         }

         open(F3,">$name.class_merge_vector_gaps") or die "Cannot create file \"$name.class_merge_vector_gaps\" : $!\n";
         for (my $x=0 ; $x <= $#tab_file1 ; $x++)
         {
             chomp $tab_file11[$x];
             chomp $tab_pb[$x];
             print F3 "$tab_pb[$x] $tab_file11[$x]\n";
         }
         close F3;



