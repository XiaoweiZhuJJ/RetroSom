#! /usr/bin/perl -w

use strict;

my $sample   = $ARGV[0];
my $retro    = $ARGV[1];
my $path     = $ARGV[2];
my $exo_file = $path.'/'.$retro.'/'.$sample.'.ALU.alignout';
my $ou1_file = $path.'/'.$retro.'/'.$sample.'.discover';
my $ou2_file = $path.'/'.$retro.'/'.$sample.'.sr.discover';
my $CUT = 98;
my ($line, @data);
my ($name, $map, $cord1, $cord2);

(open (EXO, "<$exo_file")) || die "cannot open the exo file\n";
(open (OU1, ">$ou1_file")) || die "cannot open the out file\n";
(open (OU2, ">$ou2_file")) || die "cannot open the out file\n";

while ($line = <EXO>)
   {
    chomp($line);
    next if ($line =~ /^\#/);
    @data = split("\t", $line);
    $name = $data[0];
    $map  = $data[2];
    $cord1= $data[8];
    $cord2= $data[9];
    
    if ($map > $CUT && ( ($cord1<245 && $cord2>260 ) || ($cord1>260 && $cord2<245) ) )
       {
        # print "$name\n";
        if ($name =~ /pe$/)
           {
            $name =~ s/\_/\t/g;
            print OU1 "$name\t$map\n";
           }
        elsif ($name =~ /sr$/)
           {
            $name =~ s/\_/\t/g;
            print OU2 "$name\t$map\n";
           }
       }
   }
close(OU1);
close(OU2);   
