#! /usr/bin/perl
#####File SNP_Rloop_CGHweb.pl

#use strict;
#use warnings;
#use Cwd;

#loops through all files in a directory
my $studyDir = $ARGV[0];
my $inputDir = $ARGV[1];
my $initial_outputDir = $ARGV[2];
my $platform = $ARGV[3];
my $studyName = $ARGV[4];

opendir DIR, $inputDir or die "cannot open dir $inputDir: $!";
my @ls=grep !/^\./, readdir(DIR);  # to ensure that . and .. directories aren't included: 
closedir DIR;

#my $_=~ s/\n//;

#open batch file to run in R
open(out, ">temp_cgh2.R") or die "cannot open out $!";
$counter = 0;

  foreach(@ls)
  {
    $_=~ s/\n//;     #removes line return after each file name      
    #R commands
    if($counter == 0) {print out "library(CGHweb)\n";}
    print out  "setwd(\"$inputDir/\")\n";
    print out "tab<-read.delim(file=\"$_\", header=TRUE, sep=\"\\t\")\n";
    print out "setwd(\"$studyDir/\")\n";
    $_=~ s/\_Input.txt/_Output/; #change filename
    my $outputDir = $initial_outputDir.$_;                                            
    print out "x<- runCGHAnalysis(tab,tempDir=getwd(),resultDir=\"$outputDir\")\n\n";
    $counter = $counter +1;
  }
  
  print "\nThere are $counter files to run in CGHweb for the $studyName study";
  close(out); #close batch file
      
exit; #!/usr/bin/perl
