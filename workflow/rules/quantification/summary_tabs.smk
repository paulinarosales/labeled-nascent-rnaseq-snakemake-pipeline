        
rule summary_ERCC:
    input:
        'results/quality_control/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.idxstats'
    output:
        idx_percentTSV = 'results/quality_control/summary_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.idxstats.percentage.tsv',
        ercc_summaryTSV = 'results/quality_control/summary_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.ERCC_content.tsv'
    log: 
        'logs/summary_tabs/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_ercc_content.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    script:
        '../../scripts/summary_tabs/ercc_content.R'
