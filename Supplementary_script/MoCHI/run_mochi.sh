#!/bin/bash
#Export all enviroment variables
#$ -V
#Use current working directory
#$ -cwd
#Join stdout and stderr
#$ -j y
#Email address
##$ -m ae
#$ -M chending@mail.ustc.edu.cn
#Memory
##$ -l virtual_free=40G
##$ -l virtual_free=20G
#Time
##$ -q long-sl7
##$ -q gpu
##$ -l h_rt=24:00:00
##$ -l gpu=1
#Parallel environment
##$ -pe smp 10
##$ -pe smp 4

source /home/user/miniconda3/etc/profile.d/conda.sh  # This will properly initialize conda
conda activate mochi-gpu  # Activate the correct environment

OUTPUT_DIR="./"

PROJECT_NAME="MoCHI_result"
L2_REG="0.000001"
seed=901
LEARN_RATE=0.025  # Set the learning rate to 0.025

# Call the Python script with the new learning rate
python3 MoCHI_setting.py --output_directory $OUTPUT_DIR --project_name $PROJECT_NAME --l2_regularization_factor $L2_REG --seed $seed --learn_rate $LEARN_RATE

