#!/usr/bin/env perl -w
#translate "*_ratio.txt" (output of FREEC) into a Circos track

use strict;

my $usage = qq{
    $0

    -----------------------------
    mandatory parameters:
    
    -f	 file            	file with ratio
    -p   ploidy                 ploidy (default 2)
    
    -----------------------------
    optional parameters:
    
    -v 				verbose mode
    
};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

##  arguments

my $filename = "";
my $verbose = 0;
my $ploidy = 2;

## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-f') {$filename = shift @ARGV;}
    elsif ( $this_arg eq '-p') {$ploidy = shift @ARGV;}   

    elsif ( $this_arg eq '-v') {$verbose = 1;}

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if( $filename eq ""){
    die "you should specify file with ratio\n";
}
   
 

#-----------read file with pairs-----
my $numberOfAllSites = 0;
my $totalCount=0;
my $lines = 0;


#open (OUT , ">$name") or die "Cannot open file $name!!!!: $!";

if ($filename eq "") {
	print "Please specify a file with ratio!\n" if ($verbose) ;
        exit();
}

  
if ($filename =~ /.gz$/) {
        open(FILE, "gunzip -c $filename |") or die "$0: can't open ".$filename.":$!\n";
} else {
        open (FILE, "<$filename") or die "Cannot open file $filename!!!!: $!";
}

#Chromosome	Start	Ratio	MedianRatio	CopyNumber
#1	1	-1	1.03966	2
#1	1001	-1	1.03966	2

# to

#hsY 0 999999 0.005786
#hsY 1000000 1999999 0.003958
#hsY 2000000 2999999 0.002578
#hsY 3000000 3999999 0.007451


my ($chr,$start,$ratio,$MedRatio,$CPN) = ("","","","","");

my $eventStart = "";
my $eventEnd = "";
my $eventValue = "";
my $eventChr = "";
my $count = 0;
my $eventRatio = "";

while (<FILE>) {		
    chomp;
    ($chr,$start,$ratio,$MedRatio,$CPN) = split /\t/;      
    next if ($chr =~ m/^Chromosome/);        
    $totalCount++;
    $count++;
    if ($eventStart eq "") {
        $eventStart = $start;
        $eventEnd = $start+999;
        $eventValue = $MedRatio;
	$eventRatio = $ratio;
        $eventChr =$chr;
	$count =0;
	print "hs$eventChr $eventStart $eventEnd ".$eventValue*$ploidy."\n" if ($eventRatio >= 0);

    } else {
        if ($eventChr ne $chr) {
            #new Chromosome, print the last event
            
            print "hs$eventChr $eventStart $eventEnd ".$eventValue*$ploidy."\n" if ($eventValue >= 0);

            $eventStart = $start;
            $eventEnd = $start+999;
            $eventValue = $MedRatio;
	    $eventRatio = $ratio;
            $eventChr =$chr;
	    $count =0;
	    print "hs$eventChr $eventStart $eventEnd ".$eventValue*$ploidy."\n" if ($eventRatio >= 0);

        } else {
             if ($eventValue ne $MedRatio) {
                #new event, print the last event
                print "hs$eventChr $eventStart $eventEnd ".$eventValue*$ploidy."\n" if ($eventValue >= 0);

                $eventStart = $start;
                $eventEnd = $start+999;
                $eventValue = $MedRatio;
		$eventRatio = $ratio;
                $eventChr =$chr;
		$count =0;
		print "hs$eventChr $eventStart $eventEnd ".$eventValue*$ploidy."\n" if ($eventRatio >= 0);

             } else {
                #move forward
                $eventEnd = $start+999;;                    
             }                         
        }
	if ($count>=10) {
	    $count =0;
	    print "hs$eventChr $eventStart $eventEnd ".$eventValue*$ploidy."\n" if ($ratio >= 0); #($eventValue >= 0); $eventRatio >0
	    $eventStart=$eventEnd;

	}
    }        
}
close FILE;

#print the last event
print "hs$eventChr $eventStart $eventEnd ".$eventValue*$ploidy."\n" if ($eventRatio >= 0);

print STDERR "Read: $filename\t$totalCount\n";


