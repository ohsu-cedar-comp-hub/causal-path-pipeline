from utils import ensure_dir, generate_data_files, generate_proteomics_data, generate_parameter_file, generate_rna_data
import pandas as pd
import os
from itertools import combinations

meta_file = snakemake.params.meta
meta = pd.read_csv(meta_file,sep='\t',index_col=0)
meta = meta.astype(str)


phospho_prot_file = snakemake.params.phospho_prot
phospho_prot = pd.read_csv(phospho_prot_file,sep='\t').drop_duplicates()
phospho_prot.ID = phospho_prot.ID.str.upper()

condition_id = snakemake.params.condition
permutations = snakemake.params.perm
fdr = snakemake.params.fdr
site_match = snakemake.params.site_match
site_effect = snakemake.params.site_effect
ds_thresh = snakemake.params.ds_thresh
rna_file = snakemake.params.rna_file

transform, ctype, cond = snakemake.output[0].split('/')[1:-1]
relnm = os.path.join(*[os.getcwd(),'results',transform, ctype, cond])
ensure_dir(relnm)
kwargs = {condition_id:list(map(str,cond.split('_')))}
sub_data, baseline, contrast = generate_data_files(phospho_prot, meta, condition_id, **kwargs)

generate_proteomics_data(sub_data, relnm)

if rna_file != None:
    print('Incorporating RNAseq into causal relations')

    rna_frame = pd.read_csv(rna_file,sep='\t',index_col=0)
    print('total RNAseq expression matrix of shape {},{}'.format(rna_frame.shape[0],rna_frame.shape[1]))
    print(rna_frame.head())
    sub_rna = rna_frame.reindex(sub_data.columns,axis=1).iloc[:,3:]
    print(sub_rna.head())
    generate_rna_data(sub_rna, relnm)

generate_parameter_file(relnm, contrast, baseline, ctype, transform, fdr, site_match, site_effect, permutations)

