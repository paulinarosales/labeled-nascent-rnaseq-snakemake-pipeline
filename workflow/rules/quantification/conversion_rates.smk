def _input_tx_infoTSV(wildcards):
    return expand('resources/external/merged_ERCC_gencode_{release}/merged_ERCC_{genome}.transcripts_info.tsv', release=GENCODE_RELEASE, genome=GENOME)

rule conv_rates:
    input:
        collapsedTSV = 'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionCounts.transcripts.tsv',
        transcriptBED = _input_tx_infoTSV
    output:
        globalTSV = 'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates.global.tsv',
        txTSV = 'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates.transcripts.tsv'
    log:
        'logs/rates/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    script:
        '../../scripts/quantification/conversion_rates.R'