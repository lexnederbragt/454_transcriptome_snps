#! /usr/bin/perl

use strict;
use warnings;

my %reads = ();	# keys: read ID (uaccno), value: sample name
my @samples;	# sorted list of unique sample IDs
my $accno;		# current contig etc

my @groups = sort qw/DOM HISP/;	# regular expressions to split samples into two
my %groups;					# will hold '0' or '1' depending on hit to @IDs

&reads;
&diffs;


sub reads{
	my %samples;	# unique sample IDs
	#print "Reading in read --> sample information\n";
	open READS, "< path/to/readID_sample.tsv" or die $!;
	while (<READS>){
		chomp;
		my ($read, $sample) = split;
		$reads{$read} = $sample;
		$samples{$sample} = 1;
	}
	close READS;
	@samples = sort keys %samples;
	# assign sample IDs to sample groups
	foreach my $sample (@samples){
		foreach my $group (@groups){
			$groups{$sample}=$group if $sample =~ /$group/
		}
# debug	print "$sample --> $groups{$sample}\n"
	}
}

sub diffs{
	open HCDIFF, "< 454HCDiffs.txt" or die $!;
	$/=">"; # set the record separator to the '>' symbol
			# this forces each sequence into $_
			# note however, that each sequence ENDS with the '>' symbol
			# and that the first 'entry' (record) is consisting of ONLY the '>' symbol
	while (<HCDIFF>){
		chomp;	# remove the trailing '>' symbol
		next if length($_) < 1;	# skip empty lines

		my @lines = split(/\n/,$_);	# split the entry into individual lines based on the newline character
		my $header = shift @lines;	# the header is the first line (now without the '>' symbol)
				
		if (/Reference/){
			print "$header";
			print "\tRef\tVar" foreach (@samples,1,2);
			print "\n";
			next
		}
		if (/Accno/){
			print "$header";
			print "\t$_\t$_" foreach @samples;
			print "\t$_"."_count\t$_"."_count" foreach @groups;
			print "\n";
			next
		}

		($accno) = split "\t", $header;
#		next if $accno ne "contig00033";	# debug
		print ">$header\t";

		my $read_type;
		my %counts;
		
		foreach my $line (@lines){
			last if $line =~ /Signal Distribution/;
			if ($line =~/Reads with Difference:/){
				$read_type="diff";
				next
			}
			if ($line =~/Other Reads:/){
				$read_type="other";
				next
			}
			if ($line =~/^(\w+)/){	# get read ID
				my $sampleID = $reads{$1};			# look up sample ID
				next if $1 eq $accno;
				next if $1 eq "reference";
				if (!$sampleID){print "died:$accno-->$1\n"}
				$counts{$read_type.$sampleID}++;	# add to count for this type and sample ID
#debug#				print "$read_type\t$1\t$reads{$1}\t$counts{$read_type.$sampleID}\n";
			}			
		}

#debug#		print "\n"; while (my ($k, $v) = each %counts) {print "$k=>$v\t"};print "\n"; next;	# debug

		# results for this variant
		my (%vars, %refs); # counts of number of samples showing reads for each group
		foreach my $sampleID (@samples){
			print join "\t", (
				$counts{"other".$sampleID}||0,
				$counts{"diff".$sampleID}||0
			);
			$vars{$groups{$sampleID}}++ if $counts{"diff".$sampleID};
			$refs{$groups{$sampleID}}++ if $counts{"other".$sampleID};
			print "\t";
		}
#		print "\t",$refs{$_}||0,"\t",$vars{$_}||0 foreach @groups;
		print join "\t", (
			$refs{$_}||0,
			$vars{$_}||0,
			""
		) foreach @groups;
		print "\n";
	}
}
