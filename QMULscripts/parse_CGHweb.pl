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


my $summaryColumn = -1; #last column of file
my $file;
my $summaryFile ="Table_of_aCGH_smoothed_profiles.txt.gz"; 
my $unzippedSummary = "Table_of_aCGH_smoothed_profiles.txt";
my $fileName;
my $out;
my $lnes;

opendir(DIR, $resultDir) or die "Sorry cannot open dir $resultDir: $!";
  my @ls=grep !/^\./, readdir(DIR); 
closedir(DIR);


#important to sort based on chr and chr-position before launching R
my %sets; my @names; my @R; my $j = 3 ; my $i = 0 ; 
my $initial_loop =0;
#keep the 4 columns with the probeID, chr, position and log2ratios
foreach my $patientDir (@ls)
{
  #unzips each patient file from CGHweb and obtains summary column
  opendir(DIR, "$resultDir/$patientDir")or die "Sorry cannot open dir $resultDir/$patientDir: $!";
  gunzip "$resultDir/$patientDir/$summaryFile"=>"$resultDir/$patientDir/$unzippedSummary" or die "gunzip failed: $GunzipError\n";
  closedir(DIR);
  $fileName = $patientDir;
  $fileName =~ s/\_Output//;

  #looks in clinical file called clinicalInfo.txt to determine if sample is patient or cell line to save in the appropriate folder
  open(CLINICAL, "$filesDir/clinicalInfo.txt") or die "Sorry cannot open dir $filesDir/clinicalInfo.txt: $!";
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
  
    #goes through each line of the summary file and adds it to a main file
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
        #initial loop is for first patient
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
        else
        {
          #need to make sure clones match for subsequent patients
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
      $i++; #next row of current patient
    }   
    $j++; #new column for next patient
    close(SUMMARY);
    close(OUT);
    $initial_loop ++; #for subsequent patients
    $lnes = $i ; #records total number of lines
  }
  close(CLINICAL)
}

###  WRITE DATA  ##################################################

open(ALLP, ">$filesDir/allPatientData.txt");
print ALLP "probeID\tchr\tPos\t";

opendir(DIR2, "$parseDir/Patient") or die "Sorry cannot open dir $parseDir/Patient: $!";

my @pt=grep !/^\./, readdir(DIR2); 
#closedir(DIR2);

my $patientcol =3;
foreach my $patientDir (@pt)
{
  $fileName = $patientDir;
  $fileName =~ s/results\_//;
  $fileName =~ s/\.txt//;
  print ALLP "$fileName\t";
  $patientcol++;
}
closedir(DIR2);
print ALLP "\n";


for ($i=0;$i<$lnes;$i++) #changed from lnes
{ 
  for ($j=0;$j<$patientcol;$j++)
  {
    if(!$R[$i][$j] or $R[$i][$j] eq " ")
    {
      print ALLP " ";
      print ALLP "\t";   
    }
    else
    {
      print ALLP "$R[$i][$j]\t";
    }
  }
  print ALLP "\n";
}
close(ALLP);

open(ALLC, ">$filesDir/allCellData.txt");
print ALLC "probeID\tchr\tPos\t";
opendir(DIR3, "$parseDir/CellLines") or die "Sorry cannot open dir $parseDir/CellLines: $!";
my @cl=grep !/^\./, readdir(DIR3); 

my $cellcol =3;

foreach my $patientDir (@cl)
{
    $fileName = $patientDir;
    $fileName =~ s/results\_//;
    $fileName =~ s/\.txt//;
    print ALLC "$fileName\t";
   $cellcol++;
}
closedir(DIR3);
print ALLC "\n";

for ($i=0;$i<$lnes;$i++)
{  
  for ($j=0;$j<$cellcol;$j++)
  {
    if(!$R[$i][$j] or $R[$i][$j] eq " ")
    {
      print ALLC " "; 
    }
    else
    {
      print ALLC "$R[$i][$j]\t";
    }
  }
  print ALLC "\n";
}
close(ALLC);
