#!/usr/bin/perl

# Usage: extractData $locs $dataFolder $retWindow $outDir $offset
#  locs: the locations to be extracted. Tab-delimited, SORTED LOW TO HIGH, no header. Format:
# 	chr# \t <unused> \t bp position to extract \t line name (identifier)
# dataFolder: folder containing ALL files to process and NO OTHER files. 
# 	All files in the folder MUST be named: chr#_datasetName_<whatever>. Used to name output files as well.
# 	Must be tab-delimited. Format:
# 	position \t value
# retWindow: number of BP on EACH side of location to search for matches
# outDir: the directory into which to output the finished files

# Example usage:
# perl extractData.pl hg19positions.csv ./inputFiles/ 2000 ./outputFiles/



use warnings;
use strict;
#NEED TO REDIRECT TO SHARED PERL
use lib '/shared/fxcsbi/nucpipeline/perl_modules/Parallel-ForkManager-1.06/lib';
use Parallel::ForkManager;

#For nuclesome occupancy, 4 seems to be optimal amout of cores
#If major pipeline is being run, should increase the cores
my $pm = Parallel::ForkManager->new(4);



#key matrices to 2d array that is list of CSV lines
sub processCSV{
	my $file = shift;
	my $i = 0;
	my %rethash;
	while(<$file>){
		chomp;
# 		Skip the header
# 		if($i == 0){
# 			$i++;
# 			next;
# 		}
		my ($chr, $TSS, $RE, $name) = split(/\t/);
		my $info = [$RE, $name];
		push(@{$rethash{$chr}}, $info);
	}
	return(%rethash);
}

#Process a folder of text files with nucleosome data
sub processFolder{
	my $dir = shift;
	my %posData = %{shift @_};
	
	my $retWindow = shift;
	my $outDir = shift;
	
	opendir(DIR, $dir) or die $!;

	
	FILELOOP:
	while (my $file = readdir(DIR)) {
		my %retData;
		# Use a regular expression to ignore files beginning with a period
		next if ($file =~ m/^\./);
		
		my $pid = $pm->start and next FILELOOP;
		my $chr = $file;
		my @chr = split(/_/, $chr);
		my $outmarker = $chr[1];
		#print "$chr[0]\n";
		my $REs = $posData{$chr[0]};
		
		$file = $dir.$file;
		
		
		open my $datafile, "<", $file or die "could not open $file";
		
		my $lineNum = 1;
		
		while(<$datafile>){
			my $areMore;
			my @line = split(/\t/);
			
			#Handling for data without line numbers. Now deprecated.
			if(scalar @line == 1){
				unshift(@line, $lineNum);
			}
			
			my @REs = @$REs;
			
			foreach my $RE (@REs){
				# Get points within window
				if($RE->[0] - $retWindow > $line[0]){
					$areMore = 1;
					last;
				}
				if($RE->[0] + $retWindow >= $line[0]){
					my $outline = $line[0] . "\t" . $line[1];
					push @{$retData{$RE->[1]}}, $outline;
				}else{
					shift @$REs;
				}
				#Confirm that there is more data to retrieve from this file
				if($RE->[0] + $retWindow + 3000 >= $line[0]){
					$areMore = 1;
				}
			}
			#Next file if no more data to get
			if(!$areMore){
				last;
			}
			
			$lineNum++;
			
		}
		for my $RE ( keys %retData ) {
			my $filename = "./" . $outDir ."$RE=$outmarker";
			open(my $fh, '>', $filename);
			print $fh @{ $retData{$RE} };
			close $fh;
		}
		$pm->finish;
	}
	closedir(DIR);
	$pm->wait_all_children;
}

my ($locs, $dataFolder, $retWindow, $outDir) = @ARGV;

open my $fh, "<", $locs or die "could not open location file";

my %posData = processCSV($fh);

close $fh;



processFolder($dataFolder, \%posData, $retWindow, $outDir);

