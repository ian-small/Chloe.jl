#!/usr/bin/perl -w
use strict;

my $in_sff = shift;
my $in_intron = shift;

#my $output;
my %rps12_hash;
open IN, "$in_sff" || die "$!\n";
while(<IN>){
	chomp;
	if (/^(rps12\S+)\//){
		#print "$1\n";
		my @info = split;
		push @{$rps12_hash{$1}}, [@info];
	}else{
		print "$_\n";
		#$output .= $_;
	}

}
close IN;

open IN, "$in_intron" || die "$!\n";
while(<IN>){
	chomp;
	if (/^(rps12\S+)\//){
		#print "$1\n";
		my @info = split;
		push @{$rps12_hash{$1}}, [@info];
	}
}
close IN;

foreach my $gid(sort keys %rps12_hash){
	#print "$gid\n";
	@{$rps12_hash{$gid}} = sort {$a->[2]<=>$b->[2]} @{$rps12_hash{$gid}};
	#@{$gff{$gene_id}{CDS}} = sort {$a->[3]<=>$b->[3] or $a->[4]<=>$b->[4]} @{$gff{$gene_id}{CDS}};

	for(my $i=0; $i<@{$rps12_hash{$gid}}; $i++){
		my $num = $i + 1;
		print "$gid\/$num\t$rps12_hash{$gid}[$i][1]\t$rps12_hash{$gid}[$i][2]\t$rps12_hash{$gid}[$i][3]\t$rps12_hash{$gid}[$i][4]\t$rps12_hash{$gid}[$i][5]\n";
	}
}
