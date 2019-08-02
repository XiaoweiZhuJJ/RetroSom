#! /usr/bin/perl -w

use strict;

my $sample = $ARGV[0];
my $workfolder = $ARGV[2];
my $retro = 'retro_v'.$ARGV[1];
my $m1_file  = $workfolder.'/'.$sample.'/'.$retro.'_0/'.$sample.'.pe.LINE.matrix';
my $m2_file  = $workfolder.'/'.$sample.'/'.$retro.'_1/'.$sample.'.pe.LINE.matrix';
my $g1_file  = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe.pred.G1.txt';
my $g2_file  = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe.pred.G2.txt';
my $g3_file  = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe.pred.G3.txt';
my $g4_file  = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe.pred.G4.txt';
my $g5_file  = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe.pred.G5.txt';
my $g6_file  = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe.pred.G6.txt';
my $g7_file  = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe.pred.G7.txt';
my $g8_file  = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe.pred.G8.txt';
my $ou1_file = $workfolder.'/'.$sample.'/'.$retro.'_0/LINE/'.$sample.'.pe2.pred.all.txt';
my $ou2_file = $workfolder.'/'.$sample.'/'.$retro.'_0/LINE/'.$sample.'.pe2.pred.summary';
my $ou3_file = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe2.pred.all.txt';
my $ou4_file = $workfolder.'/'.$sample.'/'.$retro.'_1/LINE/'.$sample.'.pe2.pred.summary';
my @order = (1,4,8,2,7,6,5,3); 
my ($line, @data);
my (%pred1, %pred2, %pred3, %pred4, %pred5, %pred6, %pred7, %pred8);
my ($i, $pred, $id, $pos, $neg, $g1, $g2, $g3, $g4, $g5, $g6, $g7, $g8);

(open (M1, "<$m1_file")) || die "cannot open the mat file\n";
(open (M2, "<$m2_file")) || die "cannot open the mat file\n";
(open (G1, "<$g1_file")) || die "cannot open the g1 file\n";
while ($line = <G1>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred1{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (G2, "<$g2_file")) || die "cannot open the g2 file\n";
while ($line = <G2>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred2{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (G3, "<$g3_file")) || die "cannot open the g3 file\n";
while ($line = <G3>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred3{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (G4, "<$g4_file")) || die "cannot open the g4 file\n";
while ($line = <G4>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred4{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (G5, "<$g5_file")) || die "cannot open the g5 file\n";
while ($line = <G5>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred5{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (G6, "<$g6_file")) || die "cannot open the g6 file\n";
while ($line = <G6>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred6{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (G7, "<$g7_file")) || die "cannot open the g7 file\n";
while ($line = <G7>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred7{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (G8, "<$g8_file")) || die "cannot open the g8 file\n";
while ($line = <G8>)
   {
    chomp($line);
    @data = split("\t",$line);
    $pred8{$data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3]} = $data[4];
   }

(open (OU1, ">$ou1_file")) || die "cannot open the out file\n";
(open (OU2, ">$ou2_file")) || die "cannot open the out file\n";
(open (OU3, ">$ou3_file")) || die "cannot open the out file\n";
(open (OU4, ">$ou4_file")) || die "cannot open the out file\n";

print OU1 "chr\tcord1\tcord2\tread\tpos\tneg\tG1\tG2\tG3\tG4\tG5\tG6\tG7\tG8\n";
<M1>;
while ($line = <M1>)
   {
    chomp($line);
    @data = split("\t", $line);
    if ( (($data[6] == 1) || ($data[7] == 1)) && ($data[24] == 0) )
       {
        $id = $data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3];
        $pos = $data[6];
        $neg = $data[7];
        $g1 = (exists($pred1{$id})) ? $pred1{$id} : 'NA';
        $g2 = (exists($pred2{$id})) ? $pred2{$id} : 'NA';
        $g3 = (exists($pred3{$id})) ? $pred3{$id} : 'NA';
        $g4 = (exists($pred4{$id})) ? $pred4{$id} : 'NA';
        $g5 = (exists($pred5{$id})) ? $pred5{$id} : 'NA';
        $g6 = (exists($pred6{$id})) ? $pred6{$id} : 'NA';
        $g7 = (exists($pred7{$id})) ? $pred7{$id} : 'NA';
        $g8 = (exists($pred8{$id})) ? $pred8{$id} : 'NA';
        print OU1 "$id\t$pos\t$neg\t$g1\t$g2\t$g3\t$g4\t$g5\t$g6\t$g7\t$g8\n";
       }
   }

print OU3 "chr\tcord1\tcord2\tread\tpos\tneg\tG1\tG2\tG3\tG4\tG5\tG6\tG7\tG8\n";
<M2>;
while ($line = <M2>)
   {
    chomp($line);
    @data = split("\t", $line);
    if ( (($data[6] == 1) || ($data[7] == 1)) && ($data[24] == 0) )
       {
        $id = $data[0]."\t".$data[1]."\t".$data[2]."\t".$data[3];
        $pos = $data[6];
        $neg = $data[7];
        $g1 = (exists($pred1{$id})) ? $pred1{$id} : 'NA';
        $g2 = (exists($pred2{$id})) ? $pred2{$id} : 'NA';
        $g3 = (exists($pred3{$id})) ? $pred3{$id} : 'NA';
        $g4 = (exists($pred4{$id})) ? $pred4{$id} : 'NA';
        $g5 = (exists($pred5{$id})) ? $pred5{$id} : 'NA';
        $g6 = (exists($pred6{$id})) ? $pred6{$id} : 'NA';
        $g7 = (exists($pred7{$id})) ? $pred7{$id} : 'NA';
        $g8 = (exists($pred8{$id})) ? $pred8{$id} : 'NA';
        print OU3 "$id\t$pos\t$neg\t$g1\t$g2\t$g3\t$g4\t$g5\t$g6\t$g7\t$g8\n";
       }
   }

close(OU1);
close(OU3);
(open (OU1, "<$ou1_file")) || die "cannot open the out file\n";
(open (OU3, "<$ou3_file")) || die "cannot open the out file\n";

my ($pos_0, $pos_1, $pos_na, $pos_total, $neg_0, $neg_1, $neg_na, $neg_total, $recall, $precision);
my ($basep_0, $basep_1, $basep_na, $basen_0, $basen_1, $basen_na, $base_recall, $base_precision);
my ($rnd_recall, $rnd_precision);
my %usage = (
   '1' => 0,
   '2' => 0,
   '3' => 0,
   '4' => 0,
   '5' => 0,
   '6' => 0,
   '7' => 0,
   '8' => 0,
  ); 

<OU1>;
while ($line = <OU1>)
   {
    chomp($line);
    @data = split("\t", $line);
    $pred = 'NA';
    for ($i=0; $i<=$#order; $i++)
       {
        if ($data[5+$order[$i]] ne 'NA')
           {
            $usage{$order[$i]} ++;
            $pred = ($data[5+$order[$i]] > 0.5) ? 1 : 0;
            last; 
           }
       }
    print OU2 "$line\t$pred\n";
   }

%usage = (
   '1' => 0,
   '2' => 0,
   '3' => 0,
   '4' => 0,
   '5' => 0,
   '6' => 0,
   '7' => 0,
   '8' => 0,
  );

<OU3>;
while ($line = <OU3>)
   {
    chomp($line);
    @data = split("\t", $line);
    $pred = 'NA';
    for ($i=0; $i<=$#order; $i++)
       {
        if ($data[5+$order[$i]] ne 'NA')
           {
            $usage{$order[$i]} ++;
            $pred = ($data[5+$order[$i]] > 0.5) ? 1 : 0;
            last;
           }
       }
    print OU4 "$line\t$pred\n";
   }

