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
#$ -N coGAPb2_s8c 
##$ -t 1-10
source activate dimsum
fastqFileDir="/data1/member/mayujing/DiMSum/FASTQ"
experimentDesignPath="Abundance2_experiment_design.txt"
gzipped=T
cutadapt5First="CTCGACACAGCAGGTCAAGAG"
cutadapt5Second="TAAGTCCTGAGCCTGTTTTGT"
cutadaptErrorRate="0.2"
outputPath="./"
projectName="Abundance_block2_Q20_initial"
wildtypeSequence="GAGTACAGTGCAATGAGGGACCAGTACATGAGGACTGGGGAGGGCTTTCTTTGTGTATTTGCCATAAATAATACTAAATCATTTGAAGATATTCACCATTATAGAGAACAAATTAAAAGAGTTAAGGACTCTGAAGATGTACCTATGGTCCTAGTAGGAAATAAATGTGATTTGCCTTCTAGAACAGTAGAC"
startStage=0
stopStage=5
vsearchMinovlen=10
numCores=8
retainItermediateFiles=T
vsearchMinQual=20

DiMSum -i $fastqFileDir -e $experimentDesignPath --vsearchMinQual $vsearchMinQual --retainIntermediateFiles $retainItermediateFiles --vsearchMinovlen $vsearchMinovlen --numCores $numCores --cutadapt5First $cutadapt5First --gzipped $gzipped --cutadapt5Second $cutadapt5Second -a $cutadaptErrorRate -o $outputPath -p $projectName -w $wildtypeSequence --stopStage $stopStage
