#!/usr/bin/perl
use warnings;
use strict;

## if dist = 0; F(d) = sum_x{N(x) * [N(x) - 1]/2};
## if dist = x > 0; F(d) = sum_x{N(x) * N(x+d)}; 

# Note by Peter: This is (mostly) not my code, the algorithm was given to me, I did however tweak it to increase
# efficiency and correctness, such as by using warnings and strict

my $inFile = $ARGV[0];
my $outFile = "../corr/" . $inFile;

open (my $INPUT1, "<$inFile");
open (my $OUTPUT1, ">$outFile");

my $max = -1;

my @all_data1;
my @all_count1;


while (<$INPUT1>) {
	chomp();

	my @data = split(/\t/);

	push(@all_data1, $data[0]);
	push(@all_count1, $data[1]);
	$max++;
}

$max = $max/2;


my %occurrence = ();

for (my $i = 0; $i <= $max; $i++) {
	$occurrence{$i} = 0;
}

for (my $j = 0; $j < @all_data1; $j++) {

	for (my $m = $j; $m < @all_data1; $m++) {
		my $dist = $all_data1[$m] - $all_data1[$j];
		if ($dist >= 0 && $dist <= $max) {
			$occurrence{$dist} = $occurrence{$dist} + $all_count1[$j] * $all_count1[$m];
		}
	}
}

foreach my $key (sort numeric keys %occurrence) {
	print $OUTPUT1 "$key\t$occurrence{$key}\n";
}

sub numeric { $a <=> $b }
