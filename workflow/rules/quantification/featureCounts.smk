def _input_refGTF(wildcards):
    return expand('resources/external/merged_ERCC_gencode_{release}/merged_ERCC_{genome}.annotation.gtf', release=GENCODE_RELEASE, genome=GENOME) # de-compressed

# rule featurecounts_matrix:
#     input:
#         bam = TARGETS['samtools_sort'],
#         genome_gtf = _input_refGTF
#     output:
#         matrix = 'results/read_counts/featureCounts_allSamples_counts.transcripts.tsv',
#         summary = 'results/read_counts/featureCounts_allSamples_counts.transcripts.tsv.summary'
#     log: 
#         'logs/featureCounts/count_matrix.log'
#     conda:
#         '../../envs/alignment/hisat3n.yaml'
#     threads: 24
#     resources:
#         mem = '64G'
#     params:
#         # count_mode = config['FEATURE_COUNTS']['COUNT_MODE'],
#         ftr_type = config['FEATURE_COUNTS']['FEATURE_TYPE'],
#         atr_type = config['FEATURE_COUNTS']['ATTRIBUTE_TYPE'],
#         extra =config['FEATURE_COUNTS']['EXTRA']
#     shell:
#         """
#             featureCounts -T {threads} -p --countReadPairs\
#             -t {params.ftr_type} -g {params.atr_type} {params.extra}\
#             -a {input.genome_gtf} -o {output.matrix} {input.bam} 2> {log} 
#         """

rule featurecounts_tx:
    input:
        filterBAM = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.bam',
        genome_gtf = _input_refGTF
    output:
        report = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.bam.featureCounts',
        matrix = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_counts.transcripts.tsv',
        summary = 'results/counts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_counts.transcripts.tsv.summary'
    log: 
        'logs/featureCounts/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_counts.transcripts.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '64G'
    params:
        # count_mode = config['FEATURE_COUNTS']['COUNT_MODE'],
        ftr_type = config['FEATURE_COUNTS']['FEATURE_TYPE'],
        atr_type = config['FEATURE_COUNTS']['ATTRIBUTE_TYPE'],
        strandness = config['FEATURE_COUNTS']['STRANDNESS'],
        extra =config['FEATURE_COUNTS']['EXTRA']
    shell:
        """
            featureCounts -T {threads} -p --countReadPairs -R CORE\
            -t {params.ftr_type} -g {params.atr_type} -s {params.strandness} {params.extra}\
            -a {input.genome_gtf} -o {output.matrix} {input.filterBAM} 2> {log} 
        """


