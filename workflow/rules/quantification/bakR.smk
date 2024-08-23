
rule bakR_mut_call:
    input:
        filterBAM = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.bam',
        # filterBAI = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.bam.bai',
        snpTXT = 'results/snps/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.snp.txt'
    output:
        tracks_dir = directory('results/track_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        mutCSV = 'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionCounts.tsv',
        cUCSV = 'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_cU.tsv'
    log:
        'logs/bakR/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_mut_call.log'
    params:
        mutType = config['HISAT3N']['BASE_CHANGE'], # Type of mutation to record (default: TC)
        minDist = config['BAKR']['MIN_BASE_DIST'], # Base distance from read-end filter (default: 5)
        minQual = config['BAKR']['MIN_BASE_Q'], # Base minimal quality filter (default: 40)
        strandedness = config['BAKR']['STRANDEDNESS'] # default='F',choices=['F','R'] Is first read forward or reverse orientation? F = forward is the default
    threads: 32
    conda:
        '../../envs/downstream/bakR.yaml'
    resources:
        mem = '240G'
    script:
        '../../scripts/quantification/mut_call_bakR.py'
        


rule bakR_merge:
    input:
        featureCountsTSV = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.bam.featureCounts',
        conversionCountsTSV = 'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionCounts.tsv'
    output:
        metaTSV = 'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionCounts.metadata.tsv',
        collapsedTSV = 'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionCounts.transcripts.tsv'
    log:
        'logs/bakR/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_merge_features.log'
    threads: 24
    resources:
        mem = '64G'
    conda:
        '../../envs/downstream/r-basic.yaml'
    script:
        '../../scripts/quantification/merge_features_and_muts.R'
