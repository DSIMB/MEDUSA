#!/usr/bin/perl -w
use strict;

help_arg() if ($#ARGV != 1);


my $file_seq=shift;
my $type=shift;

my $num_classes=20;


# 0- Create hash table
my @tab_class=("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p");

my $verbose = 1;

my %hash_class;

if ($type eq 'PB')
{

	$num_classes=16;
	%hash_class=(
	"a" => 0,
	"b" => 1,
	"c" => 2,
	"d" => 3,
	"e" => 4,
	"f" => 5,
	"g" => 6,
	"h" => 7,
	"i" => 8,
	"j" => 9,
	"k" => 10,
	"l" => 11,
	"m" => 12,
	"n" => 13,
	"o" => 14,
	"p" => 15,
	"Z" => 16
);
}
elsif ($type eq 'AA')
{

$num_classes=20;


%hash_class=(
	"A" => 0,
	"C" => 1,
	"D" => 2,
	"E" => 3,
	"F" => 4,
	"G" => 5,
	"H" => 6,
	"I" => 7,
	"K" => 8,
	"L" => 9,
	"M" => 10,
	"N" => 11,
	"P" => 12,
	"Q" => 13,
	"R" => 14,
	"S" => 15,
	"T" => 16,
	"V" => 17,
	"W" => 18,
	"Y" => 19,
);
}
elsif  ($type eq 'ACC')
{

$num_classes=10;
	%hash_class=(
	"a" => 0,
	"b" => 1,
	"c" => 2,
	"d" => 3,
	"e" => 4,
	"f" => 5,
	"g" => 6,
	"h" => 7,
	"i" => 8,
	"j" => 9,
	"k" => 10,
	"l" => 11,
	"m" => 12,
	"n" => 13,
	"o" => 14,
	"p" => 15,
);



}

# 1- OPEN AND FILTERED PIR
#=========================
open(F,"$file_seq") or die "Cannot open $file_seq : $!\n";
my @tab_file=<F>;
close F;
my $current_name;
my $current_seq;

my @tab_name;
my @tab_seq;
my $number_seq=-1;
for (my $i=0; $i <= $#tab_file ; $i++)
{
	my $line=$tab_file[$i];
	next if ($line =~/^\s+/);
	chomp $line;

	if ($line =~/^>(.*)$/)
	{
		$number_seq++;
		$tab_name[$number_seq]="$1";
		$tab_seq[$number_seq]="";

	}
	#elsif ($line =~/^sequence|^structure/)
	#{
	#}
	else
	{		
		#if ($line =~/\*$/) {chop $line;}
		
		$tab_seq[$number_seq].=$line;
		if ($type eq 'AA') {$tab_seq[$number_seq]=uc($tab_seq[$number_seq]);};
	}
}
if ($number_seq==-1)
{
	print STDERR "No sequence in file!\n";
}

#############
#Parse position
my @tab_seq_tab_position;
my $number_pos=-1;
for (my $i=0 ; $i <= $number_seq ; $i++)
{
	my @tab_position=split('',$tab_seq[$i]);
	$number_pos=-1;
	print ">$tab_name[$i]\n";
	for (my $j=0 ; $j <= $#tab_position ; $j++)
	{
		$number_pos++;
		if ($type eq 'AA') {$tab_position[$j]=uc($tab_position[$j]);};
		my $id=$hash_class{$tab_position[$j]};
		my $onehot="";
		for (my $k=0 ; $k < $num_classes ; $k++)
		{
			if ($k == $id)
			{
				$onehot.=sprintf("%d ","1"); 
			}
			else
			{
				$onehot.=sprintf("%d ","0"); 
			}
		}
		chop $onehot;
		print "$onehot\n";
	}
}
