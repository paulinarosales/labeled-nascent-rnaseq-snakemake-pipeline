rule conv_rates_plot:
    input:
        'results/conversion_tables/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates.transcripts.tsv'
    output:
        'results/figures/conversions/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates.pdf'
    params:
        base_change = config['HISAT3N']['BASE_CHANGE']
    log:
        'logs/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    threads: 24
    resources:
        mem = '24G'
    script:
        '../../scripts/figures/conversion_rates_boxplot.R'

rule global_conv_rates_plot:
    input:
        sample_manifestTSV = config['SAMPLE_MANIFEST'],
        convFiles = TARGETS['global_conv_rates']
    output:
        all_globalRatesTSV = 'results/conversion_tables/all_samples/global_conversionRates.tsv',
        barplotPDF = 'results/figures/conversions/all_samples/global_conversionRates.pdf'
    log:
        'logs/figures/all_samples/global_conversionRates.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    script:
        '../../scripts/figures/global_conversion_rates_barplot.R'

rule ercc_conv_rates_plot:
    input:
        TARGETS['conv_rates']
    output:
        'results/figures/conversions/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates.ERCCs.pdf'
    params:
        base_change = config['HISAT3N']['BASE_CHANGE'],
        use_genes = config['CONVERSION_RATES']['ERCC_PLOT']
    log:
        'logs/figures/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_conversionRates_ERCCs.log'
    conda:
        '../../envs/downstream/r-basic.yaml'
    threads: 24
    resources:
        mem = '24G'
    script:
        '../../scripts/figures/conversion_rates_barplot_gois.R'
