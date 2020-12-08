#!/usr/bin/perl -w
use strict;


my $file=shift;


open(F,"$file") or die "Cannot open file $file :$!\n";

my %hash_seq;
my @tab_id;
while(my $line=<F>)
{
    next if ($line =~/^#/);
    my $size = split(/\s+/,$line);
    next if ($size != 2);
    chomp $line;
    my ($id,$seq)=split(/\s+/,$line);
    if (exists $hash_seq{$id})
    {
        $hash_seq{$id}.=$seq;
    }
    else
    {
        $hash_seq{$id}=$seq;
        push @tab_id,$id;
    }
}
close F;

foreach my $id (@tab_id)
{
    print ">$id\n";
    my @tab_seq=split('',$hash_seq{$id});

    my $outseq="";
    my $num=1;
    for(my $i = 0 ; $i <= $#tab_seq ; $i++)
    {
        if ($num == 80) 
        {
            $outseq.=$tab_seq[$i];
            print "$outseq\n";
            $outseq="";
            $num=1;
        } 
        else 
        {
            $outseq.=$tab_seq[$i];
            $num++;
        }
    }
    if ($outseq ne '')
    {
        print "$outseq\n";
    }
}


