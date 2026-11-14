#!/usr/bin/env python

#######################################################################
## IMPORTS ##
#######################################################################

import os
import pandas as pd
import argparse

import pymochi
from pymochi.project import MochiProject

#######################################################################
## Functions ##
#######################################################################

def init_argparse() -> argparse.ArgumentParser:
    """
    Command line argument parser.

    :returns: ArgumentParser.
    """

    parser = argparse.ArgumentParser(
        description=""
    )
    parser.add_argument('--seed', type = int, default = 1, help = "random seed (default: 1)")
    parser.add_argument('--learn_rate', type=float, default=0.05, help="Learning rate for the optimizer")
    parser.add_argument('--features', default = [], help = "features file")
    parser.add_argument('--output_directory', help = "output directory")
    parser.add_argument('--project_name', help = "project name")
    parser.add_argument('--ensemble', type = bool, help = "ensemble")
    parser.add_argument('--max_interaction_order', type = int, help = "max interaction order")
    parser.add_argument('--downsample_observations', type = float, default = None, help = "downsample observations")
    parser.add_argument('--downsample_interactions', default = None, type = str, help = "downsample interactions")
    parser.add_argument('--l1_regularization_factor', default = "0", help = "l1 regularization factor")
    parser.add_argument('--l2_regularization_factor', default = "0", help = "l2 regularization factor")
    return parser

#######################################################################
## COMMANDLINE ARGUMENTS ##
#######################################################################

model_design = pd.DataFrame({
    'phenotype': [
        'Abundance1', 'Abundance2_1', 'Abundance2_2','Abundance3_1','Abundance3_2',
        'Binding1_K13', 'Binding2_K13','Binding3_K13',
        'Binding1_K19', 'Binding2_K19','Binding3_K19',
        'Binding1_K27', 'Binding2_1_K27','Binding2_2_K27', 'Binding3_1_K27','Binding3_2_K27',
        'Binding1_K55', 'Binding2_1_K55','Binding2_2_K55', 'Binding3_1_K55','Binding3_2_K55',
        'Binding1_PI3', 'Binding2_PI3','Binding3_PI3',
        'Binding1_RAF', 'Binding2_1_RAF', 'Binding2_2_RAF','Binding3_1_RAF','Binding3_2_RAF',
        'Binding1_RAL', 'Binding2_RAL','Binding3_RAL',
        'Binding1_SOS', 'Binding2_SOS','Binding3_SOS',
    ],
    'transformation': ['TwoStateFractionFolded', 'TwoStateFractionFolded', 'TwoStateFractionFolded','TwoStateFractionFolded','TwoStateFractionFolded',
                       'ThreeStateFractionBound', 'ThreeStateFractionBound', 'ThreeStateFractionBound',
                       'ThreeStateFractionBound', 'ThreeStateFractionBound', 'ThreeStateFractionBound',
                       'ThreeStateFractionBound', 'ThreeStateFractionBound', 'ThreeStateFractionBound','ThreeStateFractionBound','ThreeStateFractionBound',
                       'ThreeStateFractionBound', 'ThreeStateFractionBound', 'ThreeStateFractionBound','ThreeStateFractionBound','ThreeStateFractionBound',
                       'ThreeStateFractionBound', 'ThreeStateFractionBound', 'ThreeStateFractionBound',
                       'ThreeStateFractionBound', 'ThreeStateFractionBound', 'ThreeStateFractionBound','ThreeStateFractionBound','ThreeStateFractionBound',
                       'ThreeStateFractionBound', 'ThreeStateFractionBound', 'ThreeStateFractionBound',
                       'ThreeStateFractionBound', 'ThreeStateFractionBound', 'ThreeStateFractionBound',
    ],
    'trait': [['Folding'], ['Folding'], ['Folding'], ['Folding'],['Folding'],
              ['Folding', 'Binding_K13'], ['Folding', 'Binding_K13'], ['Folding', 'Binding_K13'],
              ['Folding', 'Binding_K19'], ['Folding', 'Binding_K19'], ['Folding', 'Binding_K19'],
              ['Folding', 'Binding_K27'], ['Folding', 'Binding_K27'], ['Folding', 'Binding_K27'],['Folding', 'Binding_K27'],['Folding', 'Binding_K27'],
              ['Folding', 'Binding_K55'], ['Folding', 'Binding_K55'], ['Folding', 'Binding_K55'],['Folding', 'Binding_K55'],['Folding', 'Binding_K55'],
              ['Folding', 'Binding_PI3'], ['Folding', 'Binding_PI3'], ['Folding', 'Binding_PI3'],
              ['Folding', 'Binding_RAF'], ['Folding', 'Binding_RAF'], ['Folding', 'Binding_RAF'],['Folding', 'Binding_RAF'],['Folding', 'Binding_RAF'],
              ['Folding', 'Binding_RAL'], ['Folding', 'Binding_RAL'], ['Folding', 'Binding_RAL'],
              ['Folding', 'Binding_SOS'], ['Folding', 'Binding_SOS'], ['Folding', 'Binding_SOS'],
    ],
    'file': [
        "CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
        "CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
        "MA_RAS_abundance_2_fitness_replicates_fullseq.RData",
        "CW_RAS_abundance_3_fitness_replicates_fullseq.RData",
        "MA_RAS_abundance_3_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K13_1_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K13_2_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K13_3_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K19_1_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K19_2_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K19_3_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K27_3_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_K55_3_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_PI3_1_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_PI3_2_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_PI3_3_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData",
        "MA_RAS_binding_RAF_3_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_RAL_1_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_RAL_2_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_RAL_3_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_SOS_1_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_SOS_2_fitness_replicates_fullseq.RData",
        "CW_RAS_binding_SOS_3_fitness_replicates_fullseq.RData",
    ],
})

holdout_orders = [2]
batch_size = "2048"

#######################################################################
## MAIN ##
#######################################################################

#Get command line arguments
parser = init_argparse()
args = parser.parse_args()

#MoCHI project
mochi_project = MochiProject(
    directory = os.path.join(args.output_directory, args.project_name),
    seed = args.seed,
    RT = 303*0.001987,
    model_design = model_design,
    downsample_observations = args.downsample_observations,
    downsample_interactions = args.downsample_interactions,
    holdout_orders = holdout_orders, 
    features = args.features,
    ensemble = args.ensemble,
    batch_size = batch_size,
    learn_rate = args.learn_rate,
    l1_regularization_factor = args.l1_regularization_factor,
    l2_regularization_factor = args.l2_regularization_factor)
