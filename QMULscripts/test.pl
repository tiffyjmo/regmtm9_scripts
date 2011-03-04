#!/usr/bin/perl

#Tiffany Morris October 2010
#This is the main script to control the Online Tool analysis of CGH array data

#these input values come from the interface
my $workingDir = "/home/analysis/onlineTool/";
#my $workingDir = "C:/Users/Tiff/Desktop/onlineTool/";

my $platform = "100KSNParray";
my $study = "user1";
#   $study = ~ s/\@/\_/;
#   $study = ~ s/\./\_/;
my $totalNormals = 1;
my $patients = 2;
my $cells = 0;
my $totalSamples = $patients + $cells;
my $posthreshold = 0.2;
my $negthreshold = -0.2;
my $ampthreshold = 0.4;
my $delthreshold = -0.4;

my $studyDir = $workingDir."cghResults/".$study;
$directory = $studyDir;
mkdir $directory;
my $analysisDir = $workingDir."aromaAnalysis/";
my $scriptsDir = $workingDir."scripts/";

#Create directories for CEL files
#this will vary for other platforms
#system('mkdir -p /home/analysis/OnlineTool/$platform/$study/rawdata/$platform/{Mapping50K_Hind240,Mapping50K_Xba240}');

#save normal and test CEL files to these directories
#put raw data in $directory = $workingDir."/Analysis/".$platform."/rawData/".$study;


#mkdir for all necessary directories for pipeline
$inputDir = $studyDir."/CGHweb_InputFiles";
mkdir $inputDir;
$resultDir = $studyDir."/CGHweb_Results";
mkdir $resultDir;
$directory = $parseDir."/parse_CGHweb_RESULTS";
mkdir $parseDir;
$directory = $studyDir."/parse_CGHweb_RESULTS/Patient";
mkdir $directory;
$directory = $studyDir."/parse_CGHweb_RESULTS/CellLines";
mkdir $directory;
$directory = $studyDir."/Regions";
mkdir $directory;

#need to create clinicalData.txt form user input
#currently saved in Harada

 
#creates an R script with to send parameters to 100SNP_Analysis.R 
#100SNP_Analysis.R uses aroma.affymetrix to analyse SNP100K data
#This works!!
  open(out, ">setVariables1.R") or die "cannot open out $!";
  print out "platform = '"; print out $platform; print out "'\n";
  print out "dataSetName = '"; print out $study; print out "'\n";
  print out "totalSamples = $totalSamples\n";
  print out "totalNormals = $totalNormals\n";
  print out "analysisDir = '"; print out $analysisDir; print out "'\n";
  print out "resultsDir = '"; print out $studyDir; print out "'\n"; 
  print out "scriptsDir = '"; print out $scriptsDir; print out "'\n";
  print out "setwd(scriptsDir)\n";
  print out "source('100SNP_Analysis.R')\n";
  system ('R --vanilla < setVariables1.R');