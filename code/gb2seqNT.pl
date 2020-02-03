#! /usr/bin/env perl
# Auth: Jennifer Chang
# Date: 2018/05/14

use strict;
use warnings;

# ===== Check ARGS 
my $USAGE="USAGE: $0 <input.gb> > <tabular|fasta>\n";
$USAGE=$USAGE."    input.gb - genbank file\n";
$USAGE=$USAGE."    Edit the printEntry function to select tab or fasta output\n";

if(@ARGV != 1){
    die $USAGE;
    exit;
}
# ===== Variables
my $fn=$ARGV[0];

my $unknown="-";

my $gb="GenBank";
my $fasta="";
my $seq=-1;

# ===== Print line of GB data
sub printEntry(){

   # = Fasta output
    $fasta=~s/[0-9]//g;
    $fasta=~s/ //g;

    print ">$gb\n";
    print $fasta,"\n";

    # = Reset variables
    $gb=$unknown;
    $fasta="";
    $seq=-1;
}

# ===== Main
my $fh;
open($fh, "<:encoding(UTF-8)",$fn)
    or die "Could not open file '$fn'";

while(<$fh>){
    if(/^\/\//){
	printEntry;
    }
    if($seq>0){
	$fasta=$fasta.$_;
    }elsif(/LOCUS\s+(\S+)/){
	$gb=$1;
    }elsif(/^ORIGIN/){
	$seq=1;
    }
}
