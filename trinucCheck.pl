#!/usr/bin/perl

use strict;
use warnings;
use Bio::DB::Fasta;

my $fasta = "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/";
my %refHash;

my $l;
my ($chr,$pos,$id,$ref,$alts,$qual,$filter,$info,$format,$g1,$g2, $start, $end, $seq);


while ($l = <STDIN>)
{
	chomp $l;

	unless ($l =~ /^#/)
	{
		($chr,$pos,$id,$ref,$alts,$qual,$filter,$info,$format,$g1,$g2) = split(/\t/, $l);

		for my $alt (split(/,/, $alts))
		{
			if (length($ref) == length($alt))
			{
				unless (exists $refHash{$chr})
				{
					$refHash{$chr} = Bio::DB::Fasta->new("$fasta/$chr.fa");
				}
				$start = $pos - 1;
				$end = $pos + 1;

				$seq = $refHash{$chr}->seq($chr, $start, $end);

				print "$chr\t$pos\t$ref\t$alt\t$seq\n";
			}
		}
	}
}


