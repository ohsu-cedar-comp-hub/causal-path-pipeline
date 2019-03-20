__author__ = "Joey Estabrook"
__email__ = "estabroj@ohsu.edu"
__license__ = "MIT"

"""CausalPath pipeline"""


import datetime
import sys
import os
import pandas as pd
import json
from itertools import combinations
from pylab import *

configfile:"omic_config.yaml"
meta_data = pd.read_csv(config['meta'],sep='\t',index_col=0)
condition = list(map(str,meta_data.loc[:,config['condition']].unique().tolist()))
combs = list(combinations(condition,2))
comb_list = ['_'.join(x) for x in combs][:4]

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

with open('cluster.json') as json_file:
    json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


value_transformation = ['significant-change-of-mean','difference-of-means','fold-change-of-mean']
data_type = ['protein_only','protein_rna','rna_only']

rule all:
    input:
        expand("results/{transform}/{type}/{cond}/ProteomicsData.txt",transform=value_transformation,type=data_type,cond=comb_list),
        expand("results/{transform}/{type}/{cond}/parameters.txt",transform=value_transformation,type=data_type,cond=comb_list),
        expand("results/correlation/{condition}/parameters.txt",condition=condition),
        expand("results/correlation/{condition}/ProteomicsData.txt",condition=condition),
        expand("results/{transform}/{type}/{cond}/causative.sif",transform=value_transformation,type=data_type,cond=comb_list),
        expand("results/correlation/{condition}/causative.sif",condition=condition)

include: "rules/causal.smk"
