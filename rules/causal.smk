
rule generate_contrast_data_files:
    params:
        meta = config['meta'],
        phospho_prot = config['phos_prot'],
        condition = config['condition'],
        perm = config['perm'],
        fdr = config['fdr'],
        site_match = config['site_match'],
        site_effect = config['site_effect']
    output:
        "results/{transform}/{type}/{cond}/ProteomicsData.txt",
        "results/{transform}/{type}/{cond}/parameters.txt"
    script:
        "../scripts/partition_data.py"        


rule generate_correlation_data_files:
    params:
        meta = config['meta'],
        phospho_prot = config['phos_prot'],
        condition = config['condition'],
        perm = config['perm'],
        fdr = config['fdr'],
        site_match = config['site_match'],
        site_effect = config['site_effect']
    output:
        "results/correlation/{condition}/parameters.txt",
        "results/correlation/{condition}/ProteomicsData.txt"
    script:
        "../scripts/partition_data_causal.py"


rule run_causal_path:
    input:
        "results/{transform}/{type}/{cond}/ProteomicsData.txt",
        "results/{transform}/{type}/{cond}/parameters.txt"
    output:
        "results/{transform}/{type}/{cond}/causative.sif"
    shell:
        "java -jar resources/causalpath/target/causalpath.jar results/{wildcards.transform}/{wildcards.type}/{wildcards.cond}"


rule run_causal_path_correlation:
    input:
        "results/correlation/{condition}/ProteomicsData.txt",
        "results/correlation/{condition}/parameters.txt"
    output:
        "results/correlation/{condition}/causative.sif"
    shell:
        "java -jar resources/causalpath/target/causalpath.jar results/correlation/{wildcards.condition}"
            
