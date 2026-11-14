#!/bin/bash 
#Export all enviroment variables 
#$ -V 
#Use current working directory 
#$ -cwd 
#Join stdout and stderr 
#$ -j y 
#Email address 
#$ -m ae 
#$ -M chending@mail.ustc.edu.cn
#Memory 
#$ -l virtual_free=120G 
##$ -l virtual_free=10G 
#Time 
#$ -q short-centos79  
#$ -l h_rt=6:00:00 
#Parallel environment 
#$ -pe smp 8 
##$ -pe smp 10
#$ -N coGAPb3_s8c 
##$ -t 1-10
source activate dimsum
experimentDesignPath="Abundance2_experiment_design.txt"
outputPath="./"
projectName="Abundance_block2_Q20_rbg_filter"
wildtypeSequence="GAGTACAGTGCAATGAGGGACCAGTACATGAGGACTGGGGAGGGCTTTCTTTGTGTATTTGCCATAAATAATACTAAATCATTTGAAGATATTCACCATTATAGAGAACAAATTAAAAGAGTTAAGGACTCTGAAGATGTACCTATGGTCCTAGTAGGAAATAAATGTGATTTGCCTTCTAGAACAGTAGAC"
startStage=4
stopStage=5
numCores=8
retainItermediateFiles=T
countPath="./variantCounts_Abundance_block2.tsv"
fitnessMinInputCountAll="1:20,2:10,3:5,4:5,5:5,6:5"
#sequenceType="coding"
#mutagenesisType="codon"
#barcodeIdentityPath="barcode_bind_concatenate_20210903.txt"
#DiMSum -i $fastqFileDir -e $experimentDesignPath --vsearchMinQual $vsearchMinQual --retainIntermediateFiles $retainItermediateFiles --vsearchMinovlen $vsearchMinovlen --numCores $numCores --cutadapt5First $cutadapt5First --gzipped $gzipped --cutadapt5Second $cutadapt5Second -a $cutadaptErrorRate -o $outputPath -p $projectName -w $wildtypeSequence --stopStage $stopStage --fitnessMinInputCountAll $fitnessMinInputCountAll --startStage $startStage --barcodeIdentityPath $barcodeIdentityPath --sequenceType $sequenceType --mutagenesisType $mutagenesisType

DiMSum  -e $experimentDesignPath --retainIntermediateFiles $retainItermediateFiles --numCores $numCores -o $outputPath -p $projectName -w $wildtypeSequence --stopStage $stopStage --startStage $startStage --countPath $countPath --fitnessMinInputCountAll $fitnessMinInputCountAll
