from utils import ensure_dir, generate_data_files, generate_proteomics_data, generate_parameter_file
import pandas as pd
import os
from itertools import combinations

meta_file = snakemake.params.meta
meta = pd.read_csv(meta_file,sep='\t',index_col=0)

phospho_prot_file = snakemake.params.phospho_prot
phospho_prot = pd.read_csv(phospho_prot_file,sep='\t')

condition_id = snakemake.params.condition
permutations = snakemake.params.perm
fdr = snakemake.params.fdr
site_match = snakemake.params.site_match
site_effect = snakemake.params.site_effect
ds_thresh = snakemake.params.ds_thresh


transform, type_,cond = snakemake.output[0].split('/')[1:-1]
relnm = os.path.join(*[os.getcwd(),'results',transform, type_, cond])
ensure_dir(relnm)
kwargs = {condition_id:list(map(int,cond.split('_')))}
sub_data, baseline, contrast = generate_data_files(phospho_prot, meta, condition_id, **kwargs)
generate_proteomics_data(sub_data, relnm)
generate_parameter_file(relnm, contrast, baseline, transform, fdr, site_match, site_effect, permutations)

