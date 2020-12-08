#!/usr/bin/perl -w
use strict;

#$pwd_init/psisearch_out/$name

my $f=shift;

my @tab_file=`ls -1 $f.it*`;

my $last=0;
foreach my $file (@tab_file)
{
	chomp $file;
	if ($file =~/^.*([0-5])$/)
	{
		#print "$1\n";
		if ($last < $1)
		{
			$last=$1;
		}	
	}	
}
print "$last";


