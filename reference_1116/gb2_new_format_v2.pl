#!/usr/bin/perl -w
use strict;

my $ingb=shift;

my %ghash;
my $genome_len;
my $accession;
open IN, "$ingb" || die "$!\n";
while(<IN>){
	chomp;

#	unless($genome_version){
#		$genome_version=$1 if(/VERSION\s+(\S+)/);
#	}

	unless($genome_len){
		if(/LOCUS\s+(\S+)\s+(\d+)\s+bp/){
			$accession = $1;
			$genome_len = $2;
			print "$accession\t$genome_len\n";
		}
	}

	if (/^\s+(CDS|tRNA|rRNA)\s+/){
		#print "$_\n";
		my $type = $1;
		#print STDERR "$type\n";

		my $s=(/complement/) ? ("-") : ("+");
		#print "$s\n";

		s/\<|\>//g;
		my @info=split(/\,/);
		my @start_positions; my @exon_lens;
		for(my $i=0; $i<@info; $i++){
			if ($info[$i]=~/(\d+)\.\.(\d+)/){
				#print "\t$1_$2";
				my $start = $1;
				my $stop = $2;
				my $exon_len = $stop - $start + 1;

				if($s eq "+"){
					#$start--;
				}else{
					#$start = $genome_len - $start;
					#$start = $genome_len - $stop;
					$start = $genome_len - $stop + 1; # 1-based
				}

				push @start_positions, $start;
				push @exon_lens, $exon_len;
			}
		}

		if($s eq '-'){
			@start_positions = reverse @start_positions;
			@exon_lens = reverse @exon_lens;
		}

		#@positions = sort {$a <=> $b } @positions;

		my $n_line=<IN>;
		if ($n_line=~/\/gene=\"(\S+)\"/){
			my $gid="$1";

			$ghash{$gid}++;
			#print STDERR "$ghash{$gid}\n";
			#my $gene_count = $ghash{$gid} - 1;
			my $gene_count = $ghash{$gid}; # 1-based

			my $exon_add = 0;
			for (my $i=0; $i<@start_positions; $i++){
				#my $phase = $exon_add % 3;
				my $phase;
				if ($exon_add % 3 == 0){
					$phase = 0;
				}elsif($exon_add % 3 == 2){
					$phase = 1;
				}elsif($exon_add % 3 == 1){
					$phase = 2;
				}
				#print STDERR "$phase";
				$exon_add += $exon_lens[$i];

				$phase = "0" if ($type ne "CDS");

				my $exon_count = $i + 1;
				print "$gid/$gene_count/$type/$exon_count\t$s\t$start_positions[$i]\t$exon_lens[$i]\t$phase\t$accession\n";

				if($i < @start_positions - 1){
					my $start = $start_positions[$i] + $exon_lens[$i];
					my $stop = $start_positions[$i+1] -1;
					my $intron_len = $stop - $start + 1;
					#print "$gid\tintron\t$start\t$stop\t$intron_len\n";
					print "$gid/$gene_count/intron/$exon_count\t$s\t$start\t$intron_len\t0\t$accession\n";
				}
			}
		}
	}
}
close IN;
#$/="\n"
