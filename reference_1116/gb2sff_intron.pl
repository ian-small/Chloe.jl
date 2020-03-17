#!/usr/bin/perl -w
use strict;

my $ingb=shift;

my %ghash;
my %ihash;
my $genome_len;
my $accession;
open IN, "$ingb" || die "$!\n";
while(<IN>){
	chomp;

#	unless($genome_version){
#		$genome_version=$1 if(/VERSION\s+(\S+)/);
#	}

	unless($genome_len){
		if(/LOCUS\s+(\S+).+?(\d+)\s+bp/){
			$accession = $1;
			$genome_len = $2;
			#print STDERR "$accession\t$genome_len\n";
		}
	}

	if (/^\s+(intron)\s+/){
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
					#$start--; # 0-based
				}else{
					#$start = $genome_len - $start;
					#$start = $genome_len - $stop; # 0-based
					#$start = $genome_len - $stop; # 0-based
					$start = $genome_len - $stop + 1; # 1-based
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
		if ($n_line=~/\/gene=\"(rps12\S+)\"/){
			my $gid="$1";
			#print "$gid\t$s\n";

			if (exists $ghash{"$gid"}){
				#print "$gid\t$s\n";
				if(exists $ghash{"$gid\_$s"}){

				}else{
					$ghash{"$gid"}++;
					$ghash{"$gid\_$s"} = 1;
				}
			}else{
				$ghash{"$gid"} = 1;
				$ghash{"$gid\_$s"} = 1;
			}

			if (exists $ihash{"$gid\_$s"}){
				#print "$gid\_$s\texists\n";
				$ihash{"$gid\_$s"} ++;
			}else{
				#print "$gid\_$s\tno\n";
				$ihash{"$gid\_$s"} = 1;
			}

			# print STDERR "$gid\_$s\t$ghash{$gid}\n";
			my $g_count = $ghash{"$gid"};
			my $i_count = $ihash{"$gid\_$s"};

			for (my $i=0; $i<@start_positions; $i++){
				print "$gid/$g_count/$type/$i_count\t$s\t$start_positions[$i]\t$exon_lens[$i]\t0\t$accession\n";

				# if($i < @start_positions - 1){
				# 	my $start = $start_positions[$i] + $exon_lens[$i];
				# 	my $stop = $start_positions[$i+1] -1;
				# 	my $intron_len = $stop - $start + 1;
				# 	#print "$gid\tintron\t$start\t$stop\t$intron_len\n";
				# 	print "$gid/$gene_count/intron/$i\t$s\t$start\t$intron_len\n";
				# }
			}
		}
	}
}
close IN;
#$/="\n";
