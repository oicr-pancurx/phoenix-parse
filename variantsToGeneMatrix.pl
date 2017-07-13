#!/usr/bin/perl

use strict;
use warnings;

my $l;
my @f;

my %genes;
my %samples;

my @headers;
my %row;

my %voi = (
	"somatic snv,nonsynonymous" => 1,
	"somatic snv,splicing" => 1,
	"somatic snv,stopgain" => 1,
	"somatic snv,stoploss" => 1,
	"somatic indel,frameshift" => 1,
	"somatic indel,nonframeshift" => 1,
	"somatic indel,stopgain" => 1,
	"somatic indel,stoploss" => 1,
	"somatic indel,splicing" => 1,
	"somatic sv,deletion breakpoint" => 1,
	"somatic sv,duplication breakpoint" => 1,
	"somatic sv,inversion breakpoint" => 1,
	"somatic sv,translocation breakpoint" => 1,
	"somatic cnv,homozygous deletion" => 1,
);



while ($l = <STDIN>)
{
	chomp $l;
	if ($l =~ /^donor/)
	{
		@headers = split(/,/, $l);
	}
	else
	{
		@f = split(/,/, $l);
		%row = ();
		for (my $i = 0; $i < scalar(@f); $i++)
		{
			$row{$headers[$i]} = $f[$i];
		}
		if (exists $voi{"$row{mutation_class},$row{mutation_type}"})
		{
			$samples{"$row{donor},$row{tumour}"}{$row{gene}}++;
			$genes{$row{gene}}++;
		}
	}
}


my @geneOrder = sort keys %genes;
#my @geneOrder = qw/KRAS TP53 CDKN2A SMAD4 MAP2K4 ARID1A RNF43 TGFBR2 KDM6A/;

print "Donor,Sample";
for my $g (@geneOrder)
{
	print ",$g";
}
print "\n";

for my $s (sort keys %samples)
{
	print $s;
	for my $g (@geneOrder)
	{
		if (exists $samples{$s}{$g})
		{
#			print ",$samples{$s}{$g}";
			print ",1";
		}
		else
		{
			print ",0";
		}
	}
	print "\n";
}














