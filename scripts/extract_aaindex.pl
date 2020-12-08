#!/usr/bin/perl -w
use strict;
use List::Util qw[min max];

help_arg() if ($#ARGV != 1);


my $file_selected_aaindex=shift;
my $file_aaindex=shift;


# 1- OPEN AND PARSE SELECTED AAINDEX
#===================================
open(F,"$file_selected_aaindex") or die "Cannot open $file_selected_aaindex: $!\n";
my %hash_selected_aaindex;
WHILE0:while(my $line=<F>)
{
	chomp $line;
	my @tab_line=split(/\s+/,$line);
	$hash_selected_aaindex{$tab_line[1]}=1;
}
close F;

# 2- OPEN AND PARSE AAINDEX
#===========================
open(F,"$file_aaindex") or die "Cannot open $file_aaindex: $!\n";
my @tab_aaindex;
my @tab_index;
my $index=0;
WHILE1:while(my $line=<F>)
{
	chomp $line;
	if ($line=~/^#/)
	{
		print "$line\n";
		@tab_index=split(/\s+/,$line);
		shift @tab_index;
		next WHILE1;
	}
	my @tab_line=split(/\s+/,$line);
	my $id = shift @tab_line;
	next WHILE1 if (not exists $hash_selected_aaindex{$id});
	print "$line\n"; 	
}
close F;


# SUB PROGRAMS
#=============
sub help_arg
{
	print "Error !\n";
	print "Usage: $0 file_selected_aaindex file_aaindex\n";
	exit 1;
}



