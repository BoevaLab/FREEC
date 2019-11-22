#!/usr/bin/env perl -w
#prints lengths of sequences from a multifasta file

use strict;

my $Filename;

if ($ARGV[0]) {
	$Filename = $ARGV[0];
}
else {
	print "Not enough arguments\n\n";
	die "!\n\n";
}

my $result_file=$Filename;
$result_file =~ s/.*\\//;
$result_file =~ s/.*\///;
$result_file="res_".$result_file;
open (FILE, "<$Filename") or die "Cannot open file!!!!: $!";
open (OUT, ">$result_file") or die "Cannot open file!!!!: $!";

$_ = <FILE>;
if (/>(.*)/) {
	my @line_splits = split /\s+/, $1;	
	print OUT $line_splits[0],"\t";
}
my $length = 0;

while (<FILE>)  {
	chomp;
	if (/>(.*)/) {
		my @line_splits = split /\s+/, $1;	
		print OUT $length,"\n";
		$length = 0;
		print OUT $line_splits[0],"\t";
	}
	else {
		$length += length($_);
	}
}

print OUT $length,"\n";

close FILE;
close OUT;
