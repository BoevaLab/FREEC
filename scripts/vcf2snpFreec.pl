#!/usr/bin/perl -w
#translate dbSNP vcf file into the format acceptable by FREEC (for dbSNP files)

use strict;

my $usage = qq{
    $0
    
    -f	 file            	dbSNP vcf file
    
####

    OUTPUT will be in the standard output!!! pipe it if you need with gzip!!!
    
};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory arguments

my $filename = "";



## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-f') {$filename = shift @ARGV;}
    

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if( $filename eq ""){
    die "you should specify VCF file from dbSNP\n";
}
   
 

#-----------read file with pairs-----

if ($filename =~ /.gz$/) {
        open(FILE, "gunzip -c $filename |") or die "$0: can't open ".$filename.":$!\n";
} else {
        open (FILE, "<$filename") or die "Cannot open file $filename!!!!: $!";
}

##INFO=<ID=MUT,Number=0,Type=Flag,Description="Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources">
#1	10389	rs373144384	AC	A	.	.	RS=373144384;RSPOS=10390;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005000002000200;WGT=1;VC=DIV;R5;ASP
#1	10433	rs56289060	A	AC	.	.	RS=56289060;RSPOS=10433;dbSNPBuildID=129;SSR=0;SAO=0;VP=0x050000020005000002000200;WGT=1;VC=DIV;R5;ASP
#1	10440	rs112155239	C	A	.	.	RS=112155239;RSPOS=10440;dbSNPBuildID=132;SSR=0;SAO=0;VP=0x050000020005000002000100;WGT=1;VC=SNV;R5;ASP

# to

#chr1	11017	G/T	+	G	rs78927064
#chr1	11019	G/T	+	T	rs79282076
#chr1	11022	A/G	+	G	rs28775022
#chr1	11022	A/G	+	G	rs62636509


my ($chr,$start,$possibilities,$strand,$refAllele,$ID) = ("","","","+","","");

my $numberOfSites = 0;
my $totalCount=0;
my $lines = 0;

my($dot,$dot2,$info,$AltAllele);

while (<FILE>) {	
    $lines++;
    next if (m/^\#\#/);
    chomp;
    ($chr,$start,$ID,$refAllele,$AltAllele,$dot,$dot2,$info) = split /\t/;      
    $totalCount++;    
    next unless (length($refAllele)==1 && length($AltAllele)==1);
    $numberOfSites++;
    unless($chr =~ m/^chr/) {
	$chr="chr".$chr;
    }
    
    $possibilities = $refAllele."/".$AltAllele;
    print "$chr\t$start\t$possibilities\t$strand\t$refAllele\t$ID\n";
}
close FILE;


print STDERR "Read: $filename\tlines: $lines;\ttotal sites: $totalCount\taccepted sites: $numberOfSites\n";


