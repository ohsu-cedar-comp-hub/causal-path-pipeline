from utils import ensure_dir, generate_data_files, generate_data_files_causal, generate_proteomics_data, generate_parameter_file
import pandas as pd
import os
from itertools import combinations

meta_file = snakemake.params.meta
meta = pd.read_csv(meta_file,sep='\t',index_col=0)

condition_id = snakemake.params.condition
permutations = snakemake.params.perm
fdr = snakemake.params.fdr
site_match = snakemake.params.site_match
site_effect = snakemake.params.site_effect

phospho_prot_file = snakemake.params.phospho_prot
phospho_prot = pd.read_csv(phospho_prot_file,sep='\t')

correlation, cond = snakemake.output[0].split('/')[1:-1]
causal_relnm = os.path.join(*[os.getcwd(),'results', 'correlation', cond])
ensure_dir(causal_relnm)
kwargs = {condition_id:list(map(str,[cond]))}
print(kwargs)

sub_data, baseline, contrast = generate_data_files_causal(phospho_prot, meta, condition_id, **kwargs)
generate_proteomics_data(sub_data, causal_relnm)
generate_parameter_file(relnm=causal_relnm, test_samps=contrast, control_samps=baseline, value_transformation='correlation', fdr_threshold=fdr, site_match=site_match, site_effect=site_effect, permutations=permutations)

