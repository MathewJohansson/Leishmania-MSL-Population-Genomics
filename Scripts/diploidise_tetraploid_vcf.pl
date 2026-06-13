#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# Script: diploidise_tetraploid_vcf.pl
# Purpose: Convert tetraploid VCF genotypes to diploid by splitting each tetraploid sample into two diploid samples
# Author: Generated for L. infantum MSL project
# Date: 2025-06-29

# PSEUDOCODE:
# INITIALIZE:
# - Parse command line parameters for allele positions (default: a=0, b=1, c=2, d=3)
# - Set input from STDIN
#
# FOR EACH LINE:
#   IF line starts with '#':
#     IF line starts with '#CHROM':
#       - Print comment about script purpose
#       - Split line into fields
#       - Keep fields 1-9 as is
#       - For fields 10+ (sample names):
#         - Duplicate each sample name with .1 and .2 suffixes
#       - Print modified header line
#     ELSE:
#       - Print comment line as is
#   ELSE:
#     - Process data line
#     - Keep fields 1-9 (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
#     - For fields 10+ (sample genotypes):
#       CALL convert_tetraploid_to_diploid() subroutine
#     - Print modified data line
#
# SUBROUTINE convert_tetraploid_to_diploid(genotype_field):
#   - Split genotype_field on ':'
#   - Extract GT (first element)
#   - Extract remaining fields (DP, AD, etc.)
#   - Split GT on '/' to get 4 alleles
#   - Create GT1 = allele[a]/allele[b]  
#   - Create GT2 = allele[c]/allele[d]
#   - Reconstruct two complete genotype fields:
#     - field1 = GT1:remaining_data
#     - field2 = GT2:remaining_data
#   - Return field1 and field2 separated by tab

# Command line parameters
my $allele_a = 0;
my $allele_b = 1;
my $allele_c = 2;
my $allele_d = 3;
my $help = 0;

GetOptions(
    'a=i' => \$allele_a,
    'b=i' => \$allele_b,
    'c=i' => \$allele_c,
    'd=i' => \$allele_d,
    'help' => \$help
) or die "Error in command line arguments\n";

if ($help) {
    print STDERR "Usage: $0 [options] < input.vcf > output.vcf\n";
    print STDERR "Options:\n";
    print STDERR "  -a INT  Position of first allele for first diploid (default: 0)\n";
    print STDERR "  -b INT  Position of second allele for first diploid (default: 1)\n";
    print STDERR "  -c INT  Position of first allele for second diploid (default: 2)\n";
    print STDERR "  -d INT  Position of second allele for second diploid (default: 3)\n";
    print STDERR "  -help   Show this help message\n";
    print STDERR "\nExample: cat tetraploid.vcf | $0 > diploid.vcf\n";
    exit 0;
}

# Subroutine to convert tetraploid genotype to two diploid genotypes
sub convert_tetraploid_to_diploid {
    my ($genotype_field) = @_;
    
    # Split genotype field on ':'
    my @gt_components = split /:/, $genotype_field;
    
    # Extract GT (first element) and remaining data
    my $gt = shift @gt_components;
    my $remaining_data = join(':', @gt_components);
    
    # Handle missing genotypes
    if ($gt eq '.' || $gt eq './././.') {
        my $diploid_missing = './.';
        my $field1 = $diploid_missing . ($remaining_data ? ":$remaining_data" : "");
        my $field2 = $diploid_missing . ($remaining_data ? ":$remaining_data" : "");
        return ($field1, $field2);
    }
    
    # Split GT on '/' to get 4 alleles
    my @alleles = split /\//, $gt;
    
    # Check if we have exactly 4 alleles for tetraploid
    if (scalar @alleles != 4) {
        die "Error: Expected tetraploid genotype (4 alleles) but found " . scalar @alleles . " alleles in: $gt\n";
    }
    
    # Create two diploid genotypes using specified allele positions
    my $gt1 = $alleles[$allele_a] . '/' . $alleles[$allele_b];
    my $gt2 = $alleles[$allele_c] . '/' . $alleles[$allele_d];
    
    # Reconstruct complete genotype fields
    my $field1 = $gt1 . ($remaining_data ? ":$remaining_data" : "");
    my $field2 = $gt2 . ($remaining_data ? ":$remaining_data" : "");
    
    return ($field1, $field2);
}

# Process input line by line
while (my $line = <STDIN>) {
    chomp $line;
    
    # Handle comment lines
    if ($line =~ /^#/) {
        if ($line =~ /^#CHROM/) {
            # Print comment about script purpose before CHROM line
            print "##diploidise_tetraploid_vcf.pl: Converted tetraploid genotypes to diploid by splitting each sample into .1 and .2 copies\n";
            print "##diploidise_allele_positions: a=$allele_a, b=$allele_b, c=$allele_c, d=$allele_d\n";
            
            # Process CHROM line
            my @fields = split /\t/, $line;
            
            # Keep fields 1-9 as is
            my @output_fields = @fields[0..8];
            
            # For fields 10+ (sample names), duplicate with .1 and .2 suffixes
            for my $i (9..$#fields) {
                push @output_fields, $fields[$i] . '.1', $fields[$i] . '.2';
            }
            
            print join("\t", @output_fields) . "\n";
        } else {
            # Print other comment lines as is
            print $line . "\n";
        }
    } else {
        # Process data lines
        my @fields = split /\t/, $line;
        
        # Keep fields 1-9 (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
        my @output_fields = @fields[0..8];
        
        # For fields 10+ (sample genotypes), convert tetraploid to diploid
        for my $i (9..$#fields) {
            my ($field1, $field2) = convert_tetraploid_to_diploid($fields[$i]);
            push @output_fields, $field1, $field2;
        }
        
        print join("\t", @output_fields) . "\n";
    }
}
