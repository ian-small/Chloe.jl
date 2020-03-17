#!/usr/bin/perl -w
use strict;

my $infile=shift;

open IN, "$infile" || die "$!\n";
my $head = <IN>;
my $id = $1 if($head =~ /^LOCUS\s+(\S+)/);
print ">$id\n";

while(<IN>){
	if (/^ORIGIN/){
		stop: while(<IN>){
			last stop if (/^\/\//);
			#print "test: $_\n";
			tr/atcgn/ATCGN/;
			my @info = split(/\s+/);
			shift @info;
			shift @info;
			print join "", @info, "\n";
		}

	}
}
close IN;
