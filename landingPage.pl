#!/usr/bin/perl

use strict;
use warnings;

my $rootDir = "/.mounts/labs/PCSI/";
my $hpc = "https://www.hpc.oicr.on.ca/archive/projects/PCSI/";

my $resultsPath = $ARGV[0];		# /analysis/oicr/donor/sample/wgs/bwa/0.6.2/results
my $sample = $ARGV[1];		# PCSI_0218_Pa_P_526
my $normal = $ARGV[2];		# PCSI_0218_Ly_R


my $normalPath = $resultsPath;
$normalPath =~ s/$sample/$normal/;

warn "opening >$rootDir/$resultsPath/${sample}_landing.html\n";
open (HTML, ">$rootDir/$resultsPath/${sample}_landing.html") or die "Couldn't open $rootDir/$resultsPath/${sample}_landing.html\n";

print HTML "<a href=\"$hpc/$resultsPath/${sample}_summary.html\">Genomics Summary</a><br>\n";
print HTML "<a href=\"$hpc/$resultsPath/${sample}_slide.html\">Genomics Slide</a><br>\n";
print HTML "<a href=\"$hpc/$resultsPath/${sample}.onco.html\">Oncology Targets Summary</a><br>\n";
print HTML "<a href=\"$hpc/$resultsPath/${sample}.trials.html\">Clinical Trials Summary</a><br>\n";
print HTML "<br>\n";
print HTML "<a href=\"$hpc/$resultsPath/../final_strelka-mutect/${sample}.final.vcf.html\">Old Somatic VCF file Summary</a><br>\n";
print HTML "<a href=\"$hpc/$resultsPath/../celluloid/v11.2/\">Celluloid Directory</a><br>\n";
print HTML "<a href=\"$hpc/$resultsPath/../collapsed/${sample}.bam\">Tumour BAM File</a><br>\n";
print HTML "<a href=\"$hpc/$normalPath/../collapsed/${normal}.bam\">Normal BAM File</a><br>\n";

close HTML;


