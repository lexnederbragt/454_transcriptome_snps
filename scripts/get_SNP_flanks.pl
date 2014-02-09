#! /usr/bin/perl

use strict;
use warnings;

# first argument: 454AllContigs.fna file
# second argument: table to process

my %seqs;
my %lengths;
my %descriptions;

open FASTA, "<$ARGV[0]" or die $!;
open TABLE, "<$ARGV[1]" or die $!;
open ANNOT, "< path/to/Taeniopygia_guttata.taeGut3.2.4.62.cdna.biomart.fa" or die $!;


# load sequences
$/=">"; # set the record separator to the '>' symbol
<FASTA>;		# remove the empty first 'sequence'
while (<FASTA>){
	chomp;	# remove the trailing '>' symbol
	my @lines = split(/\n/,$_);	# split the entry into individual lines based on the newline character
	my $header = shift @lines;	# the header is the first line (now without the '>' symbol)
	# take the first 'word'
    ($header) = split /\s+/, $header;
	my $seq = join "", @lines;

	#build hashes of sequences and their lengths
	$seqs{$header} = $seq;
	$lengths{$header} = length($seq);
}
$/="\n"; # reset the record separator

#Load descriptions, e.g.
#>ENSTGUG00000000018|ENSTGUT00000000020|Vacuolar protein sorting-associated protein 11 homolog (hVPS11)(RING finger protein 108) [Taeniopygia guttata]
while (<ANNOT>){
	chomp;
	$descriptions{$1} = $3 || "None" if (/>(ENSTGUG\d+)\|(ENSTGUT\d+)\|?(.*)/)
}


# proces table
while (<TABLE>){
	chomp;
	my $left_flank;
	my $snp;
	my $right_flank;
	my $description;
	
	my @cols = split;
	# 0: ENSTGUG_reference from first mapping, e.g. 110621_DOM_vs_Tg_ml90%mi95
	# 1: ENSTGUG_map_start where contig starts relative to ENSTGUG
	# 2: ENSTGUG_map_end end
	# 3: DOM_mapping_contig_number from first mapping, e.g. 110621_DOM_vs_Tg_ml90%mi95
	# 4: DOM_contig_start start position of the variant in the contig
	# 5: DOM_contig_end end
	# 6: Ref_Nuc
	# 7: Var_Nuc

	# first line: prepare header
	if ($. == 1){
		$left_flank = "left_flank";
		$snp = "SNP";
		$right_flank = "right_flank";
		$description = "description";
	}
	else{
		# Initialize
		my $start = $cols[4] - 100;
		my $stop = $cols[5] + 100;
		
		# check for enough bases before  and after SNP
		$start = 1 if $start < 1;
		$stop = $lengths{$cols[3]} if $stop > $lengths{$cols[3]};
		$snp = "[$cols[6]/$cols[7]]";
		my $seq = $seqs{$cols[3]};
		$left_flank = substr($seq, $start-1, 100);
		$right_flank = substr($seq, $cols[5],$stop-$cols[5]);
		$description = $descriptions{$cols[0]};
	}

	print join "\t", (
		@cols[3..5],
		$left_flank,
		$snp,
		$right_flank,
		@cols[8..14],
		@cols[0..2],
		$description,
		);	
	print "\n";
}