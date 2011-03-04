#!/usr/bin/perl

#Tiffany Morris June 2010

use strict;
use warnings;

my $filesDir = $ARGV[0];
my $platform = $ARGV[1];
my $studyName = $ARGV[2];
my $patientSamples = $ARGV[3];
my $cellSamples = $ARGV[4];
           
my @allColumns;
my $cellDelete = 6;
my $cellAmp = 5;
my $patientDelete = 4;
my $patientAmp = 3;

my $lossBegin;
my $lossEnd;
my $lossChr;
my $gainBegin;
my $gainEnd;
my $gainChr;
my $delBegin;
my $delEnd;
my $delChr;
my $ampBegin;
my $ampEnd;
my $ampChr;

my $gainLoop = 0;
my $lossLoop = 0;
my $ampLoop = 0;
my $delLoop = 0;

my @gainRegions;
my @lossRegions;
my @ampRegions;
my @delRegions;

my @allGainRegions;
my @allLossRegions;
my @allAmpRegions;
my @allDelRegions;
my @freq;

my $chr;
my @chromosomes = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22');

my $cellLoss = $patientSamples+$cellSamples+3+3;
my $cellGain = $patientSamples+$cellSamples+3+2;
my $patientLoss = $patientSamples+$cellSamples+3+1;
my $patientGain = $patientSamples+$cellSamples+3;

foreach $chr (@chromosomes) 
{
print "this chrom $chr\n";
open(FREQ, "$filesDir/AllFreq.txt") || die "Sorry I couldn't open $filesDir/AllFreq.txt file";
while(<FREQ>)
{
  chomp;
  @allColumns=split;

if($allColumns[0] eq "probeID")
{
  @gainRegions = "";
  @lossRegions = "";
  next;
}
else
{
  if($chr eq $allColumns[1])
  {
    if(($allColumns[$patientGain]>(.2*$patientSamples)) or ($allColumns[$cellGain]>(.2*$cellSamples)))
    {
      if($gainLoop == 0)
      {
        $gainBegin = $allColumns[2];                                                                      
      }
      $gainChr = $allColumns[1];
      $gainEnd = $allColumns[2];
      $gainLoop++; 
    }
    else
    {
      if($gainLoop > 0)
      {
        push(@gainRegions, $gainChr, "\t", $gainBegin, "\t",$gainEnd, "\t",$gainLoop, "\n");
        push(@allGainRegions, $gainChr, "\t", $gainBegin, "\t",$gainEnd, "\t",$gainLoop, "\n");
        $gainLoop = 0;
      } 
    }
    if(($allColumns[$patientLoss]>(.2*$patientSamples)) or ($allColumns[$cellLoss]>(.2*$cellSamples)))
    {
      if($lossLoop == 0)
      {
        $lossBegin = $allColumns[2];
      }
      $lossChr = $allColumns[1];
      $lossEnd = $allColumns[2];
      $lossLoop++; 
    }
    else
    {
      if($lossLoop > 0)
      {
        push(@lossRegions, $lossChr, "\t", $lossBegin, "\t",$lossEnd, "\t",$lossLoop, "\n");
        push(@allLossRegions, $lossChr, "\t", $lossBegin, "\t",$lossEnd, "\t",$lossLoop, "\n");
        $lossLoop = 0;
      }
    }
    next;
  }
  else
  {
    next;
  }  
}
}
if($gainLoop > 0)
{
  push(@gainRegions, $gainChr, "\t", $gainBegin, "\t",$gainEnd, "\t",$gainLoop, "\n");
  push(@allGainRegions, $gainChr, "\t", $gainBegin, "\t",$gainEnd, "\t",$gainLoop, "\n");
  $gainLoop = 0;
}
if($lossLoop > 0)
{
  push(@lossRegions, $lossChr, "\t", $lossBegin, "\t",$lossEnd, "\t",$lossLoop, "\n");
  push(@allLossRegions, $lossChr, "\t", $lossBegin, "\t",$lossEnd, "\t",$lossLoop, "\n");
  $lossLoop = 0;
}
close(FREQ);

open(GAIN, ">$filesDir/Regions/computedGainRegions_$chr.txt") || die "Sorry cannot open computedGainRegions.txt file";
print GAIN "CHR\tSTARTPOS\tENDPOS\tConsecReg\n@gainRegions";
close(GAIN);                
print "Data for gains on chromosome has been saved.\n";

open(LOSS, ">$filesDir//Regions/computedLossRegions_$chr.txt") || die "Sorry cannot open computedLossRegions.txt file";
print LOSS "CHR\tSTARTPOS\tENDPOS\tConsecReg\n@lossRegions";
close(LOSS);               
print "Data for losses on chromosome has been saved.\n";


open(AMPFREQ, "$filesDir/AmpDelFreq.txt") || die "Sorry I couldn't open AmpDelFreq.txt file";
while(<AMPFREQ>)
{
  chomp;
  @allColumns=split;
if($allColumns[0] eq "probeID")
{
  @ampRegions = "";
  @delRegions = "";
  next;
}
else
{
  if($chr eq $allColumns[1])
  {
    if(($allColumns[$patientAmp]>0) or ($allColumns[$cellAmp]>0))
    {
      if($ampLoop == 0)
      {
        $ampBegin = $allColumns[2];
      }
      $ampChr = $allColumns[1];
      $ampEnd = $allColumns[2];
      $ampLoop++; 
    }
    else
    {
      if($ampLoop > 0)
      {
        push(@ampRegions, $ampChr, "\t", $ampBegin, "\t",$ampEnd, "\t",$ampLoop, "\n");
        push(@allAmpRegions, $ampChr, "\t", $ampBegin, "\t",$ampEnd, "\t",$ampLoop, "\n");
        $ampLoop = 0;
      }
    } 
    if(($allColumns[$patientDelete]>0) or ($allColumns[$cellDelete]>0))
    {
      if($delLoop == 0)
      {
        $delBegin = $allColumns[2];
      }
      $delChr = $allColumns[1];
      $delEnd = $allColumns[2];
      $delLoop++; 
    }
    else
    {
      if($delLoop > 0)
      {
        push(@delRegions, $delChr, "\t", $delBegin, "\t",$delEnd, "\t",$delLoop, "\n");
        push(@allDelRegions, $delChr, "\t", $delBegin, "\t",$delEnd, "\t",$delLoop, "\n");
        $delLoop = 0;
      }
    }
    next;
  }
  next;  
}
}
if($ampLoop > 0)
{
  push(@ampRegions, $ampChr, "\t", $ampBegin, "\t",$ampEnd, "\t",$ampLoop, "\n");
  push(@allAmpRegions, $ampChr, "\t", $ampBegin, "\t",$ampEnd, "\t",$ampLoop, "\n");
  $ampLoop = 0;
}
if($delLoop > 0)
{
  push(@delRegions, $delChr, "\t", $delBegin, "\t",$delEnd, "\t",$delLoop, "\n");
  push(@allDelRegions, $delChr, "\t", $delBegin, "\t",$delEnd, "\t",$delLoop, "\n");
  $delLoop = 0;
}
close(AMPFREQ);

open(AMP, ">$filesDir/Regions/computedAmpRegions_$chr.txt") || die "Sorry cannot open computedAmpRegions.txt file";
print AMP "CHR\tSTARTPOS\tENDPOS\tConsecReg\n@ampRegions";
close(AMP);               
print "Data for amps on chromosome has been saved.\n";

open(DEL, ">$filesDir/Regions/computedDelRegions_$chr.txt") || die "Sorry cannot open computedDelRegions.txt file";
print DEL "CHR\tSTARTPOS\tENDPOS\tConsecReg\n@delRegions";               
close(DEL);
print "Data for deletions on chromosome has been saved.\n";
} 


open(ALLGAIN, ">$filesDir/Regions/computedGainRegions_ALL.txt") || die "Sorry cannot open computedGainRegions.txt file";
print ALLGAIN "CHR\tSTARTPOS\tENDPOS\tConsecReg\n@allGainRegions";
close(ALLGAIN);                
print "Data for all gains has been saved.\n";

open(ALLLOSS, ">$filesDir/Regions/computedLossRegions_ALL.txt") || die "Sorry cannot open computedLossRegions.txt file";
print ALLLOSS "CHR\tSTARTPOS\tENDPOS\tConsecReg\n@allLossRegions";
close(ALLLOSS);               
print "Data for all losses has been saved.\n"; 

open(ALLAMP, ">$filesDir/Regions/computedAmpRegions_ALL.txt") || die "Sorry cannot open computedAmpRegions.txt file";
print ALLAMP "CHR\tSTARTPOS\tENDPOS\tConsecReg\n@allAmpRegions";
close(ALLAMP);               
print "Data for all amps has been saved.\n";

open(ALLDEL, ">$filesDir/Regions/computedDelRegions_ALL.txt") || die "Sorry cannot open computedDelRegions.txt file";
print ALLDEL "CHR\tSTARTPOS\tENDPOS\tConsecReg\n@allDelRegions";               
close(ALLDEL);
print "Data for all deletions has been saved.\n";     



if($platform eq "SNP100")
{
  my @filteredFreq;
  open(FREQ, "$filesDir/AllFreq.txt") or die "Sorry cannot open file $filesDir/$platform/$studyName.AllFreq!";
  while(<FREQ>)
  {
    chomp $_;
    @freq  = split/\t/;
    
    if($freq[0] eq "probeID")
    {
      foreach (@freq)
      { 
        push(@filteredFreq, $_,"\t");
      }
      pop(@filteredFreq);
      push(@filteredFreq, "\n");   
    }
    else
    {
      open(GAINREGIONS, "$filesDir/Regions/computedGainRegions_ALL.txt") or die "Sorry cannot open file $filesDir/$platform/$studyName/Regions/computedGainRegions_ALL!";
      while(<GAINREGIONS>)
      {
        chomp $_;
        my @gregion  = split/\t/;
        if($gregion[0] eq "CHR")
        {
          next;
        }
        elsif(($freq[1] == $gregion[0]) && ($gregion[3]>20))
        {
          if(($freq[2] == $gregion[1]) || (($freq[2]>$gregion[1]) && ($freq[2]<$gregion[2])))
          {
            #save data to an array
            foreach (@freq)
            { 
              push(@filteredFreq, $_,"\t");
            }
            pop(@filteredFreq);
            push(@filteredFreq, "\n");            
          }
          else
          {
          next;
          }   
        }
        else
        {
        next;
        }
      }
      open(LOSSREGIONS, "$filesDir/Regions/computedLossRegions_ALL.txt") or die "Sorry cannot open file $filesDir/$platform/$studyName/Regions/computedLossRegions_ALL!";
      while(<LOSSREGIONS>)
      {
        chomp $_;
        my @lregion  = split/\t/;
        if($lregion[0] eq "CHR")
        {
          next;
        }
        elsif(($freq[1] == $lregion[0]) && ($lregion[3]>20))
        {
          if(($freq[2] == $lregion[1]) || (($freq[2]>$lregion[1]) && ($freq[2]<$lregion[2])))
          {
            #save data to an array
            foreach (@freq)
            { 
              push(@filteredFreq, $_,"\t");
            }
            pop(@filteredFreq);
            push(@filteredFreq, "\n");            
          }
          else
          {
          next;
          }   
        }
        else
        {
        next;
        }
      }
      
    }
  }
  open(PLOT, ">$filesDir/Regions/filteredRegionsFreqPlot.txt") || die "Sorry cannot open filteredRegionsFreqPlot.txt file";
  #print PLOT "probeID\tchr\tPos\tpatient_Freq_Gain\tpatient_Freq_Loss\tCell_Freq_Gain\tCell_Freq_Loss\n@filteredFreq";
  print PLOT @filteredFreq;
  close(PLOT);               
  print "Data for freq plot has been saved.\n";
}                                                                                                                                              
