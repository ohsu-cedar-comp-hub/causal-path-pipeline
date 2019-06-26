import pandas as pd
import os
import sys

def generate_proteomics_data(sub_data, relnm):
    """ Write out phospho-proteomices datafile to designated output directory

    Args:
        sub_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
        relnm (str): out path for analysis files

    """

    sub_data.iloc[:, 4:] = sub_data.iloc[:, 4:].fillna('NaN')
    sub_data.to_csv(path_or_buf=os.path.join(relnm, 'ProteomicsData.txt'),sep='\t')


def generate_rna_data(rna_data, relnm):
    """ Write out RNAseq data to accompany phosphoproteomic data

    Args:
        rna_data (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
        relnm (str): out path for analysis files
    """
    rna_data.to_csv(path_or_buf=os.path.join(relnm, 'rna_file.txt'), sep='\t')
    

def generate_data_files(phospho_prot, merged_meta, condition, **kwargs):
    """ Subset and index phosphoproteomic data to generate a comparison using two factor levels in the condition
         column of the meta data file

    Args:
        phospho_prot (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
        merged_meta (:obj: `pandas DataFrame`) : pandas DataFrame containing experiment specific meta data of
            shape [n_samples, ]
    Returns:
        sub_data (:obj: `pandas DataFrame`) : pandas DataFrame containing subset of processed phosphoproteomic data
        pre_rx (:obj: `list`): list of boolean values to index control/baseline samples
        post_rx (:obj: `list`): list of boolean values to index treatment/contrast samples

    """
    print(kwargs)
    print('Generating contrast based analysis files')
    subset_meta = merged_meta[merged_meta[condition].isin(kwargs[condition])]
    samps_idx = ['ID','Symbols','Sites','Effect'] + subset_meta.index.tolist()
    pre_rx = subset_meta[(subset_meta[condition] == kwargs[condition][0])].index
    post_rx = subset_meta[(subset_meta[condition] == kwargs[condition][1])].index
    sub_data = phospho_prot.loc[:, samps_idx]
    sub_data.set_index('ID',inplace=True)

    return sub_data, pre_rx, post_rx


def generate_data_files_causal(phospho_prot, merged_meta, condition, **kwargs):
    """
    Args:
        phospho_prot (:obj: `pandas DataFrame`) : pandas DataFrame containing processed phospho data
        merged_meta (:obj: `pandas DataFrame`) : pandas DataFrame containing experiment specific meta data of
            shape [n_samples, ]
    Returns:
        sub_data (:obj: `pandas DataFrame`) : pandas DataFrame containing subset of processed phosphoproteomic data
        pre_rx (:obj: `list`): list of boolean values to index control/baseline samples
        post_rx (:obj: `list`): list of boolean values to index treatment/contrast samples
 
   """
    print(kwargs)
    print('Generating correlation based analysis files')
    subset_meta = merged_meta[merged_meta[condition].isin(kwargs[condition])]
    samps_idx = ['ID','Symbols','Sites','Effect'] + subset_meta.index.tolist()
    pre_rx = subset_meta[(subset_meta[condition] == kwargs[condition][0])].index
    post_rx = None 
    sub_data = phospho_prot.loc[:, samps_idx]
    sub_data.set_index('ID',inplace=True)
    
    return sub_data, pre_rx, post_rx


def ensure_dir(relnm):
    """ Accept relative filepath string, create it if it doesnt already exist
        return filepath string

    Args:
        relnm (str) : Relative name/path

    Returns:
        relnm (str)

    """

    d = os.path.join(os.getcwd(), relnm)
    if not os.path.exists(d):
        os.makedirs(d)


def generate_parameter_file(relnm, test_samps, control_samps, ctype, value_transformation = 'significant-change-of-mean',
                            fdr_threshold = '0.1', ds_thresh='1.0', site_match = '5', site_effect = '5', permutations=1000):
    """ Generate the required CausalPath input parameters file that associates the ProteomicsData file with the
        parameters for the analysis

    Args:
        relnm (str): out path for analysis files
        test_samps (list): list of samples to be used as the test group or contrast in the analysis
        control_samps (list): list of samples to be used as the control group or baseline in the analysis
        value_transformation (str): means of detecting changes between the test and control group
        fdr_threshold (str): False discovery rate threshold
        ds_thresh (str): Threshold for data significance
        site_match (str): Associate sites with known site if site falls within n=site_match of known site
        site_effect (str): Associate sites with known activity of site if site falls within n=site_match of known site

    """
    panda_out = os.path.join(os.getcwd(),'resources/panda')
    initial_lines = ["proteomics-values-file = ProteomicsData.txt", "id-column = ID", "symbols-column = Symbols",
                     "sites-column = Sites", "effect-column = Effect", "do-log-transform = false"]
    out_f = open(os.path.join(relnm, 'parameters.txt'), 'w')
    for v in initial_lines:
        out_f.write(v + '\n')

    if value_transformation == 'significant-change-of-mean':
        sig_dif_mean_params(out_f, test_samps, control_samps, panda_out, value_transformation, fdr_threshold = fdr_threshold, site_match = site_match, site_effect = site_effect, permutations = permutations, ctype = ctype)

    if value_transformation == 'difference-of-means':
        dif_mean_params(out_f, test_samps, control_samps, panda_out, value_transformation, fdr_threshold = fdr_threshold, ds_thresh = ds_thresh, site_match = site_match, site_effect = site_effect, permutations = permutations, ctype = ctype)        
    
    if value_transformation == 'fold-change-of-mean':
        fc_mean_params(out_f, test_samps, control_samps, panda_out, value_transformation, fdr_threshold = fdr_threshold, ds_thresh = ds_thresh, site_match = site_match, site_effect = site_effect, permutations = permutations, ctype = ctype)        

    if value_transformation == 'correlation':
        correlation_params(out_f, control_samps, panda_out, value_transformation, fdr_threshold = fdr_threshold, site_match = site_match, site_effect = site_effect, permutations = permutations)
    else:
        print('value_transformation : {} not supported'.format(value_transformation))
        sys.exit()


def sig_dif_mean_params(out_f, test_samps, control_samps, panda_out, value_transformation = 'significant-change-of-mean',
                        fdr_threshold = '0.1', site_match = '5', site_effect = '5', permutations=1000, ctype=None):
    """ Generate the required CausalPath input parameters file that associates the ProteomicsData file with the
        parameters for the analysis

    Args:
        out_f (_io.TextIOWrapper): out file handle for analysis files
        test_samps (list): list of samples to be used as the test group or contrast in the analysis
        control_samps (list): list of samples to be used as the control group or baseline in the analysis
        panda_out (str): path to causalpath stored resources
        value_transformation (str): means of detecting changes between the test and control group
        fdr_threshold (str): False discovery rate threshold
        site_match (str): Associate sites with known site if site falls within n=site_match of known site
        site_effect (str): Associate sites with known activity of site if site falls within n=site_match of known site
        permutations (int): number of random data permutations to see if the result network is large, or any protein's
            downstream is enriched

    """
    out_f.write('fdr-threshold-for-data-significance = {} protein'.format(fdr_threshold) + '\n')
    out_f.write('fdr-threshold-for-data-significance = {} phosphoprotein'.format(fdr_threshold) + '\n')
    out_f.write('value-transformation = {}'.format(value_transformation) + '\n')
    out_f.write('minimum-sample-size = 3' + '\n')
    out_f.write('calculate-network-significance = true' + '\n')
    out_f.write('permutations-for-significance = {}'.format(permutations) + '\n')
    out_f.write('color-saturation-value = 2' + '\n')
    out_f.write('site-match-proximity-threshold = {}'.format(site_match) + '\n')
    out_f.write('site-effect-proximity-threshold = {}'.format(site_effect) + '\n')
    out_f.write('show-insignificant-data = false' + '\n')
    out_f.write('use-network-significance-for-causal-reasoning = true' + '\n')
    out_f.write('custom-resource-directory = {}'.format(panda_out) + '\n')
    if ctype == 'protein_rna':
        out_f.write('data-type-for-expressional-targets = rna' + '\n')
        out_f.write('data-type-for-expressional-targets = protein' + '\n')
        out_f.write('rna-expression-file = rna_file.txt' + '\n')
    if ctype == 'rna_only':
        out_f.write('data-type-for-expressional-targets = rna' + '\n')
        out_f.write('rna-expression-file = rna_file.txt' + '\n') 
    for c in control_samps:
        out_f.write('control-value-column = ' + c + '\n')
    for t in test_samps:
        out_f.write('test-value-column = ' + t + '\n')

    out_f.close()


def dif_mean_params(out_f, test_samps, control_samps, panda_out, value_transformation = 'difference-of-means',
                        fdr_threshold = '0.1', ds_thresh='1.0', site_match = '5', site_effect = '5', permutations=1000):
    """ Generate the required CausalPath input parameters file that associates the proteomic data file with the
        parameters for the analysis

    Args:
        out_f (_io.TextIOWrapper): out file handle for analysis files
        test_samps (list): list of samples to be used as the test group or contrast in the analysis
        control_samps (list): list of samples to be used as the control group or baseline in the analysis
        panda_out (str): path to causalpath stored resources
        value_transformation (str): means of detecting changes between the test and control group
        fdr_threshold (str): False discovery rate threshold
        ds_thresh (str): Threshold for data significance
        site_match (str): Associate sites with known site if site falls within n=site_match of known site
        site_effect (str): Associate sites with known activity of site if site falls within n=site_match of known site
                permutations (int): number of random data permutations to see if the result network is large, or any protein's
            downstream is enriched

    """
    out_f.write('fdr-threshold-for-data-significance = {} protein'.format(fdr_threshold) + '\n')
    out_f.write('fdr-threshold-for-data-significance = {} phosphoprotein'.format(fdr_threshold) + '\n')
    out_f.write('threshold-for-data-significance = {} protein'.format(ds_thresh) + '\n')
    out_f.write('threshold-for-data-significance = {} phosphoprotein'.format(ds_thresh) + '\n')
    out_f.write('value-transformation = {}'.format(value_transformation) + '\n')
    out_f.write('minimum-sample-size = 3' + '\n')
    out_f.write('calculate-network-significance = true' + '\n')
    out_f.write('permutations-for-significance = {}'.format(permutations) + '\n')
    out_f.write('color-saturation-value = 2' + '\n')
    out_f.write('site-match-proximity-threshold = {}'.format(site_match) + '\n')
    out_f.write('site-effect-proximity-threshold = {}'.format(site_effect) + '\n')
    out_f.write('show-insignificant-data = false' + '\n')
    out_f.write('use-network-significance-for-causal-reasoning = true' + '\n')
    out_f.write('custom-resource-directory = {}'.format(panda_out) + '\n')
    if ctype == 'protein_rna':
        out_f.write('data-type-for-expressional-targets = rna' + '\n')
        out_f.write('data-type-for-expressional-targets = protein' + '\n')
        out_f.write('rna-expression-file = rna_file.txt' + '\n')
    if ctype == 'rna_only':
        out_f.write('data-type-for-expressional-targets = rna' + '\n')
        out_f.write('rna-expression-file = rna_file.txt' + '\n')
    for c in control_samps:
        out_f.write('control-value-column = ' + c + '\n')
    for t in test_samps:
        out_f.write('test-value-column = ' + t + '\n')

    out_f.close()


def fc_mean_params(out_f, test_samps, control_samps, panda_out, value_transformation = 'fold-change-of-mean',
                        fdr_threshold = '0.1', ds_thresh='1.0', site_match = '5', site_effect = '5', permutations=1000):
    """ Generate the required CausalPath input parameters file that associates the proteomic data file with the
        parameters for the analysis

    Args:
        out_f (_io.TextIOWrapper): out file handle for analysis files
        test_samps (list): list of samples to be used as the test group or contrast in the analysis
        control_samps (list): list of samples to be used as the control group or baseline in the analysis
        panda_out (str): path to causalpath stored resources
        value_transformation (str): means of detecting changes between the test and control group
        fdr_threshold (str): False discovery rate threshold
        ds_thresh (str): Threshold for data significance
        site_match (str): Associate sites with known site if site falls within n=site_match of known site
        site_effect (str): Associate sites with known activity of site if site falls within n=site_match of known site
                permutations (int): number of random data permutations to see if the result network is large, or any protein's
            downstream is enriched

    """
    out_f.write('fdr-threshold-for-data-significance = {} protein'.format(fdr_threshold) + '\n')
    out_f.write('fdr-threshold-for-data-significance = {} phosphoprotein'.format(fdr_threshold) + '\n')
    out_f.write('threshold-for-data-significance = {} protein'.format(ds_thresh) + '\n')
    out_f.write('threshold-for-data-significance = {} phosphoprotein'.format(ds_thresh) + '\n')
    out_f.write('value-transformation = {}'.format(value_transformation) + '\n')
    out_f.write('minimum-sample-size = 3' + '\n')
    out_f.write('calculate-network-significance = true' + '\n')
    out_f.write('permutations-for-significance = {}'.format(permutations) + '\n')
    out_f.write('color-saturation-value = 2' + '\n')
    out_f.write('site-match-proximity-threshold = {}'.format(site_match) + '\n')
    out_f.write('site-effect-proximity-threshold = {}'.format(site_effect) + '\n')
    out_f.write('show-insignificant-data = false' + '\n')
    out_f.write('use-network-significance-for-causal-reasoning = true' + '\n')
    out_f.write('custom-resource-directory = {}'.format(panda_out) + '\n')
    if ctype == 'protein_rna':
        out_f.write('data-type-for-expressional-targets = rna' + '\n')
        out_f.write('data-type-for-expressional-targets = protein' + '\n')
        out_f.write('rna-expression-file = rna_file.txt' + '\n')
    if ctype == 'rna_only':
        out_f.write('data-type-for-expressional-targets = rna' + '\n')
        out_f.write('rna-expression-file = rna_file.txt' + '\n')
    for c in control_samps:
        out_f.write('control-value-column = ' + c + '\n')
    for t in test_samps:
        out_f.write('test-value-column = ' + t + '\n')

    out_f.close()


def correlation_params(out_f, control_samps, panda_out, value_transformation = 'fold-change-of-mean',
                        fdr_threshold = '0.1', site_match = '5', site_effect = '5', permutations=1000):
    """ Generate the required CausalPath input parameters file that associates the proteomic data file with the
        parameters for the analysis

    Args:
        out_f (_io.TextIOWrapper): out file handle for analysis files
        control_samps (list): list of samples to be used as the control group or baseline in the analysis
        panda_out (str): path to causalpath stored resources
        value_transformation (str): means of detecting changes between the test and control group
        fdr_threshold (str): False discovery rate threshold
        site_match (str): Associate sites with known site if site falls within n=site_match of known site
        site_effect (str): Associate sites with known activity of site if site falls within n=site_match of known site
                permutations (int): number of random data permutations to see if the result network is large, or any protein's
            downstream is enriched

    """

    out_f.write('fdr-threshold-for-correlation = {}'.format(fdr_threshold) + '\n')
    out_f.write('generate-data-centric-graph = true' + '\n')
    out_f.write('show-insignificant-data = false' + '\n')
    out_f.write('custom-resource-directory = {}'.format(panda_out) + '\n')
    out_f.write('value-transformation = {}'.format(value_transformation) + '\n')
    out_f.write('site-match-proximity-threshold = {}'.format(site_match) + '\n')
    out_f.write('site-effect-proximity-threshold = {}'.format(site_effect) + '\n')
    out_f.write('permutations-for-significance = {}'.format(permutations) + '\n')
    for c in control_samps:
        out_f.write('value-column = ' + c + '\n')
    out_f.close()

