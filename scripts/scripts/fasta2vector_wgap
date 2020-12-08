#!/usr/bin/perl -w
use strict;
use List::Util qw[min max];

help_arg() if ($#ARGV != 1);


my $file_pir=shift;
my $file_aaindex=shift;


# 1- OPEN AND PARSE AAINDEX
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
		@tab_index=split(/\s+/,$line);
		shift @tab_index;
		next WHILE1;
	}
	my @tab_line=split(/\s+/,$line);
	shift @tab_line;
	#print "@tab_line\n";
	my @tab_values;
	my %hash_index;
	for(my $i=0; $i <= $#tab_line ; $i++)
	{
		push @tab_values,$tab_line[$i];
	}
	#my $ref_tab_values=zscore(\@tab_values);
	my $ref_tab_values=minmax(\@tab_values);
	for(my $i=0; $i <= $#{$ref_tab_values} ; $i++)
	{
		
		${$tab_aaindex[$index]}{"$tab_index[$i]"}=sprintf("%-6.3f",$$ref_tab_values[$i]);
		${$tab_aaindex[$index]}{"$tab_index[$i]"}=~s/\s+//g;
	}
	$index++;
}
close F;

#print_index(\@tab_aaindex);


# 2- OPEN AND FILTERED PIR
#=========================
open(F,"$file_pir") or die "Cannot open $file_pir : $!\n";
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
		$tab_seq[$number_seq]=uc($tab_seq[$number_seq]);
	}
}
$number_seq++;

if ($number_seq==0)
{
	print STDERR "No sequence in file!\n";
}

#############
#Parse position
my @tab_aa=split('',"ACDEFGHIKLMNPQRSTVWY-_");
my @tab_seq_tab_position;
my $number_pos=-1;
for (my $i=0 ; $i < $number_seq ; $i++)
{
	my @tab_position=split('',$tab_seq[$i]);
	$number_pos=-1;
	for (my $j=0 ; $j <= $#tab_position ; $j++)
	{
		$number_pos++;
		$tab_position[$j]=uc($tab_position[$j]);
		${$tab_seq_tab_position[$i]}[$j]=$tab_position[$j];
	}
}
$number_pos++;
###############
#my %hash_position_filtered;
for (my $j=0 ; $j < $number_pos ; $j++)
{
	#$hash_position_filtered{$j}=0;
	for (my $i=0 ; $i < $number_seq ; $i++)
	{
		#print "@{$tab_seq_tab_position[$i]}\n";
		my $aa="";
		if (defined ${$tab_seq_tab_position[$i]}[$j])
		{
			$aa=${$tab_seq_tab_position[$i]}[$j];
		}
		else
		{
			$aa="-";
		}
		$aa=uc($aa);
		if ($aa =~/[ACDEFGHIKLMNPQRSTVWY]/)
		{
		#	print "0"
		}
		#else
		#{
		#	$hash_position_filtered{$j}++;
		#	print "1";
		#}
	}
	#print "\n";
}

# 3- OPEN AND PARSE PIR 
#=====================
open(F,"$file_pir") or die "Cannot open $file_pir : $!\n";

my $id="";
my $second_pirline=0;
$index=0;
my $indice=1;
my $printer;


my @tab_line_f=<F>;
close F;
foreach my $line (@tab_line_f)
{
	chomp $line;
	if ($line =~/^\>/)
	{
		$id=$line;
		$index=0;
		$indice=1;
		print "$line\n";
		$printer="";
	}
	else 
	{
		#$printer.=sprintf("+1 ");
		my @tab_line=split('',$line);
		for (my $i=0 ; $i<= $#tab_line ; $i++)
		{
			my $aa=$tab_line[$i];
			FORJ:for (my $j=0 ; $j<= $#tab_aaindex ; $j++)
     				{
					#my $indice=$j+$index;
	     				if ($aa =~/[ACDEFGHIKLMNPQRSTVWY]{1}/)
	     				{
		     				my $value=${$tab_aaindex[$j]}{$aa};	
						$printer.=sprintf("$value ");
	     				}
					elsif ($aa =~/-/)
                                        {
                                                my $value="-1";
                                                $printer.=sprintf("$value ");
                                        }
					else
					{
						print STDERR "Unknow amino acids... skip (warning sequence file is reduced by one)\n";
					}
	     			$indice++;
				}
			$index++;
			#TMP     JCG RETURN LINE
			#$printer.=sprintf("\n");
			#END TMP JCG RETURN LINE
 			#if (length($printer) >0) {print "$printer\n";}
			print "$printer\n";
			$printer="";

		}
	}

}
close F;


# SUB PROGRAMS
#=============
sub help_arg
{
	print "Error !\n";
	print "Usage: $0 file_fasta file_aaindex\n";
	exit 1;
}


sub print_index
{
	my $ref_tab_aaindex=shift;
	foreach my $ref_hash (@{$ref_tab_aaindex})
	{
		foreach my $key (keys %{$ref_hash})
		{
			print "$key:${$ref_hash}{$key} "; 
		}
		print "\n";
	}

}

sub mean
{
	my $ref_tab=shift;
	my $sum=0;
	my $num=0;
	foreach my $value (@{$ref_tab})
	{
		$sum=$sum+$value;
		$num++;
	}

	return($sum/$num);
}

sub sd
{
	my $ref_tab=shift;
	my $mean=mean($ref_tab);
	my $sum=0;
	my $sd=0;
	my $num=0;
	foreach my $value (@{$ref_tab})
	{
		$sd=$sd+(($value-$mean)*($value-$mean));
		$num++;
	}

	return(sqrt($sd/$num));
}

sub zscore
{
	my $ref_tab=shift;

	my $mean=mean($ref_tab);
	my $sd  =sd($ref_tab);

	for (my $i=0 ; $i <=$#{$ref_tab} ; $i++)
	{ 
		${$ref_tab}[$i]=(${$ref_tab}[$i]-$mean)/$sd;
	}
	return($ref_tab);
}

sub minmax
{
        my $ref_tab=shift;

        my $mean=mean($ref_tab);
        my $sd  =sd($ref_tab);
        for (my $i=0 ; $i <=$#{$ref_tab} ; $i++)
        {
                ${$ref_tab}[$i]=(${$ref_tab}[$i]-$mean)/$sd;
        }
	my $min =min(@$ref_tab);
	my $max =max(@$ref_tab);
	#print "$max $min \n";
	for (my $i=0 ; $i <=$#{$ref_tab} ; $i++)
        {
		${$ref_tab}[$i]=(${$ref_tab}[$i]-$min)/($max-$min);

	}
        return($ref_tab);
}
