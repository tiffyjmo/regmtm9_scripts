#!/usr/bin/perl
#Claude Chelala March 2008
#Adapted by Tiffany Morris December 2009

use strict;
use warnings;
use Cwd;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $filesDir = $ARGV[0];
my $resultDir = $ARGV[1];
my $parseDir = $ARGV[2];
my $platform = $ARGV[3];
my $studyName = $ARGV[4];
my $posthreshold = $ARGV[5];
my $negthreshold = $ARGV[6];
my $ampthreshold = $ARGV[7];
my $delthreshold = $ARGV[8];

#######
my $summaryColumn = -1; #last column of file
my $file;
my $summaryFile ="Table_of_aCGH_smoothed_profiles.txt.gz"; 
my $unzippedSummary = "Table_of_aCGH_smoothed_profiles.txt";
my $fileName;
my $out;
my $lnes;


opendir(DIR, "$resultDir") or die "Sorry cannot open dir $resultDir: $!";
  my @ls=grep !/^\./, readdir(DIR); 
closedir(DIR);


#important to sort based on chr and chr-position before launching R
my %sets; my @names; my @R; my $j = 3 ; my $i = 0 ; 
my $initial_loop =0;
#keep the 4 columns with the probeID, chr, position and log2ratios
foreach my $patientDir (@ls)
{
  opendir(DIR, "$resultDir/$patientDir")or die "Sorry cannot open dir $resultDir/$patientDir: $!";
  gunzip "$resultDir/$patientDir/$summaryFile"=>"$resultDir/$patientDir/$unzippedSummary" or die "gunzip failed: $GunzipError\n";
  closedir(DIR);
  $fileName = $patientDir;
  $fileName =~ s/\_Output//;

  open(CLINICAL, "$filesDir/clinicalData.txt") or die "Sorry cannot open dir $filesDir/clinicalData: $!";
  while(<CLINICAL>)
  {
    chomp $_;
    my @sample  = split/\t/;
   
    if($sample[0] eq "fileID")
    {
      next;
    }
    elsif($sample[0] eq $fileName)
    {
      $fileName = $sample[1];
      if($sample[2] eq "Cell")
      {
        $out ="CellLines/results_".$fileName.".txt";
      }
      elsif($sample[2] eq "Patient")
      {
        $out ="Patient/results_".$fileName.".txt";
      }
      else
      {
        print "No Patient/CellLine information!!"; 
        $out ="results_".$fileName.".txt" ;
      }
    }
    else{next;}
    $file = $resultDir."/".$patientDir."/".$unzippedSummary;
    print "\nProcessing file : $file\n\n";
    print "Writing matrix for sort interval in $out\n\n\n";
  
    open(SUMMARY, $file) || die "Sorry cannot open $file: $!";
    open(OUT, ">$parseDir/$out") || die "File named '$parseDir/$out' not found ! \n";
    $i=0;
    while (<SUMMARY>)
    {
      chomp $_;
      my @line =  split/\t/;
      if($line[0] eq "ProbeID") 
      {
        next;
      }
      elsif($line[1] eq 'X')
      {
        next;
      }
      elsif($line[1] eq 'Y')
      {
        next;
      }
      elsif($line[1] eq '23')
      {
        next;
      }
      elsif($line[1] eq '24')
      {
        next;
      }
      else
      {
        #first study
        if($initial_loop == 0)
        {
          $R[$i][0] = $line[0];
          $R[$i][1] = $line[1];
          $R[$i][2] = $line[2]; #check clones are the same
          if(!$line[$summaryColumn])
          {
            $line[$summaryColumn]=" ";
          }
          $R[$i][$j] = $line[$summaryColumn];
        }
        #subsequent studies
        else
        {
          if($R[$i][0] = $line[0])
          {
            $R[$i][1] = $line[1];
            $R[$i][2] = $line[2];
            if(!$line[$summaryColumn])
            {
              $line[$summaryColumn]=" ";
            }
            $R[$i][$j]=$line[$summaryColumn];
          }
          else
          {
            print "ERROR: Clone doesn't match!!!"; next;
          }
        } 
        print OUT "$line[0]\t$line[1]\t$line[2]\t$line[$summaryColumn]\n";
      }
      $i++;
    }   
    $j++;
    close(SUMMARY);
    close(OUT);
    $initial_loop ++;
    $lnes = $i ; 
  }
  close(CLINICAL)
}

open(FREQ, ">$filesDir/freq_patient.txt");
print FREQ "probeID\tchr\tPos\t";
open(P_BINARY, ">$filesDir/patient_binary.txt");
print P_BINARY "probeID\tchr\tPos\t";

opendir(DIR2, "$parseDir/Patient") or die "Sorry cannot open dir $parseDir/Patient: $!";

my @pt=grep !/^\./, readdir(DIR2); 

my $patientcol =3;
foreach my $patientDir (@pt)
{
  $fileName = $patientDir;
  $fileName =~ s/results\_//;
  $fileName =~ s/\.txt//;
  print FREQ "$fileName\t";
  print P_BINARY "$fileName\t";
  $patientcol++;
}
closedir(DIR2);
print FREQ "patient_Freq_Amp\tpatient_Freq_Gain\tpatient_Freq_Loss\tpatient_Freq_Del\n";
print P_BINARY "\n";

for ($i=0;$i<$lnes;$i++) #changed from lnes
{
  #per_row calculate frequencies
  my $gain = 0; 
  my $loss = 0;
  my $amp = 0;
  my $del = 0;
  
for ($j=0;$j<$patientcol;$j++)
{
  if(!$R[$i][$j] or $R[$i][$j] eq " ")
  {
    print FREQ " ";
    print FREQ "\t";
    print P_BINARY " ";
    print P_BINARY "\t";   
  }
  else
  {
    print FREQ "$R[$i][$j]\t";
    if($j<=2)
    {
      print P_BINARY "$R[$i][$j]\t";
    }
    elsif($j >2 and $R[$i][$j] > $ampthreshold) 
    {
      $amp++;
      $gain++;
      print P_BINARY "2\t";
    }
    elsif($j >2 and $R[$i][$j] < $delthreshold) 
    {
      $del++;
      $loss++;
      print P_BINARY "-2\t";
    }
    elsif($j > 2 and $R[$i][$j] > $posthreshold) 
    {
      $gain++;
      print P_BINARY "1\t";
    }
    elsif($j >2 and $R[$i][$j] < $negthreshold) 
    {
      $loss++;
      print P_BINARY "-1\t";
    }
    else
    {
      print P_BINARY "0\t";
    }
  }
}
print FREQ "$amp\t$gain\t$loss\t$del\n";
print P_BINARY "\n";
}
close(FREQ);
close(P_BINARY);

open(FREQCELL, ">$filesDir/freq_cell.txt");
print FREQCELL "probeID\tchr\tPos\t";
open(C_BINARY, ">$filesDir/cell_binary.txt");
print C_BINARY "probeID\tchr\tPos\t";

opendir(DIR3, "$parseDir/CellLines") or die "Sorry cannot open dir $parseDir/CellLines: $!";
my @cl=grep !/^\./, readdir(DIR3); 

my $cellcol =3;

foreach my $patientDir (@cl)
{
#  opendir(DIR3, "$filesDir/parse_CGHweb_RESULTS/CellLines/$patientDir")or die "Sorry cannot open dir $filesDir/parse_CGHweb_RESULTS/CellLines/$patientDir: $!";
    $fileName = $patientDir;
    $fileName =~ s/results\_//;
    $fileName =~ s/\.txt//;
    print FREQCELL "$fileName\t";
    print C_BINARY "$fileName\t";
  $cellcol++;
}
closedir(DIR3);
print FREQCELL "cell_Freq_Amp\tcell_Freq_Gain\tcell_Freq_Loss\tcell_Freq_Del\n";
print C_BINARY "\n";

for ($i=0;$i<$lnes;$i++) #changed from lnes
{
  #per_row calculate frequencies
  my $gain = 0 ; 
  my $loss = 0 ;
  my $amp = 0;
  my $del = 0;
  
  for ($j=0;$j<$cellcol;$j++)
  {
  if(!$R[$i][$j] or $R[$i][$j] eq " ")
  {
    print FREQCELL " ";
    print FREQCELL "\t";
    print C_BINARY " ";
    print C_BINARY "\t";   
  }
  else
  {
    print FREQCELL "$R[$i][$j]\t";
    if($j<=2)
    {
    print C_BINARY "$R[$i][$j]\t";
    }
    elsif($j >2 and $R[$i][$j] > $ampthreshold) 
    {
    $amp++;
    $gain++;
    print C_BINARY "2\t";
    }#pdac
    elsif($j >2 and $R[$i][$j] < $delthreshold) 
    {
    $del++;
    $loss++;
    print C_BINARY "-2\t";
    }#pdac
    elsif($j > 2 and $R[$i][$j] > $posthreshold) 
    {
    $gain++;
    print C_BINARY "1\t";
    }#pdac
    elsif($j >2 and $R[$i][$j] < $negthreshold) 
    {
    $loss++;
    print C_BINARY "-1\t";
    }#pdac
    else
    {
    print C_BINARY "0\t";
    }
  }
  }
  print FREQCELL "$amp\t$gain\t$loss\t$del\n";
  print C_BINARY "\n";
}
close(FREQCELL);
close(C_BINARY);
