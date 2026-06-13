#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(shuffle);

# Script: diploidise_tetraploid_vcf.pl
# Purpose: Convert tetraploid VCF genotypes to diploid by splitting each tetraploid sample into two diploid samples
# Author: Generated for L. infantum MSL project
# Date: 2025-06-29
# Version: 2.0: assigns allele positions at random for each sample
#
# Usage: cat tetraploid.vcf | $0 > diploid.vcf

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
#   - Randomly shuffle the 4 alleles for each sample
#   - Create GT1 = allele[0]/allele[1] (first two shuffled alleles)
#   - Create GT2 = allele[2]/allele[3] (last two shuffled alleles)
#   - Reconstruct two complete genotype fields:
#     - field1 = GT1:remaining_data
#     - field2 = GT2:remaining_data
#   - Return field1 and field2 separated by tab

# No command line parameters needed - using random assignment

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
    
    # Randomly shuffle the alleles for each sample
    @alleles = shuffle(@alleles);
    
    # Create two diploid genotypes using first two and last two shuffled alleles
    my $gt1 = $alleles[0] . '/' . $alleles[1];
    my $gt2 = $alleles[2] . '/' . $alleles[3];
    
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
            print "##diploidise_tetraploid_vcf_v2.pl: Converted tetraploid genotypes to diploid by randomly splitting alleles for each sample into .1 and .2 copies\n";
            print "##diploidise_method: Random allele assignment for each sample\n";
            
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
