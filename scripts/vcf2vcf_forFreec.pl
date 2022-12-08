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

my $minFreq = 0.05;

my $keepSingleOnly = 1;

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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

#or

#chr1	10177	rs367896724	A	T	.	.	RS=367896724;RSPOS=10177;dbSNPBuildID=138;SSR=0;SAO=0;VP=0x050000020005170026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;G5A;G5;KGPhase3;CAF=0.5747,0.4253;COMMON=1;TOPMED=0.76728147298674821,0.23271852701325178
#chr1	10352	rs555500075	T	G	.	.	RS=555500075;RSPOS=10352;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005170026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;G5A;G5;KGPhase3;CAF=0.5625,0.4375;COMMON=1;TOPMED=0.86356396534148827,0.13643603465851172
#chr1	10616	rs376342519	CCGCCGTTGCAAAGGCGCGCCG	C	.	.	RS=376342519;RSPOS=10617;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005040026000200;GENEINFO=DDX11L1:100287102;WGT=1;VC=DIV;R5;ASP;VLD;KGPhase3;CAF=0.006989,0.993;COMMON=1
#chr1	11012	rs544419019	C	G	.	.	RS=544419019;RSPOS=11012;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005150024000100;GENEINFO=DDX11L1:100287102;WGT=1;VC=SNV;R5;ASP;VLD;G5;KGPhase3;CAF=0.9119,0.08806;COMMON=1
#chr1	11063	rs561109771	T	G	.	.	RS=561109771;RSPOS=11063;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005040026000100;GENEINFO=DDX11L1:100287102;WGT=1;VC=SNV;R5;ASP;VLD;KGPhase3;CAF=0.997,0.002995;COMMON=1;TOPMED=0.99760289245667686,0.00239710754332313

# to (Freq >0.01 by default; $minFreq) and keep only single nucleotide variants ($keepSingleOnly)

##INFO=<ID=MUT,Number=0,Type=Flag,Description="Is mutation (journal citation, explicit fact): a low frequency variation that is cited in journal and other reputable sources">
#chr1	10177	rs367896724	A	T	.	.	RS=367896724;TOPMED=0.76728147298674821,0.23271852701325178
#chr1	10352	rs555500075	T	G	.	.	RS=555500075;TOPMED=0.86356396534148827,0.13643603465851172


my ($chr,$start,$possibilities,$refAllele,$ID, $freqString) = ("","","","+","","");

my $numberOfSites = 0;
my $totalCount=0;
my $lines = 0;

my($dot,$dot2,$info,$AltAllele);

while (<FILE>) {	
    $lines++;
    if (m/^\#/) {print $_ ; next;}
    chomp;
    ($chr,$start,$ID,$refAllele,$AltAllele,$dot,$dot2,$info) = split /\t/;      
    $totalCount++;    
    if ($keepSingleOnly) {next unless (length($refAllele)==1 && length($AltAllele)==1);}
    my @infoElements = split /TOPMED=/, $info;
    my $isPrint = 0;
    my $isAlleleFreqProvided = 0;
    if (scalar(@infoElements)>1) {
	$isAlleleFreqProvided = 1;
	my @infoElements2 = split (',' , $infoElements[1]);
#	print $infoElements2[0],"\n";
	if ($infoElements2[0]<(1-$minFreq)) {
		$isPrint=1;
	}

    }
    if (!$isPrint) {
	@infoElements = split /CAF=/, $info;
	if (scalar(@infoElements)>1) {
		$isAlleleFreqProvided = 1;
		my @infoElements2 = split (',' , $infoElements[1]);
		if ($infoElements2[0]<(1-$minFreq)) {
			$isPrint=1;
		}
	}
    }
    if ($isPrint || !$isPrint&&!$isAlleleFreqProvided) { 	#checked whether TOPMED or CAF were actually provided
    	unless($chr =~ m/^chr/) {
		$chr="chr".$chr;
    	}
    	$numberOfSites++;
    	print $_,"\n";
    }
}
close FILE;


print STDERR "Read: $filename\tlines: $lines;\ttotal sites: $totalCount\taccepted sites: $numberOfSites\n";
