# --------- Parsing functions ----------
def _params_trim_galore(wildcards):
    # add --2colour option depending on sequencing technology
    sample_type = wildcards['sample_type']
    treatment = wildcards['treatment']
    bio_rep = wildcards['bio_rep']

    _sequencer = SAMPLES.loc[( 
                          sample_type, 
                          treatment,
                          bio_rep), 
                        'Sequencer'][0]

    qual_th = config['TRIM_GALORE']['QUALITY']

    if _sequencer == 'HiSeq2500':
        quality = f'--quality {qual_th}'

    elif _sequencer == 'HiSeq4000':
        quality = f'--quality {qual_th}'
        
    elif _sequencer == 'NovaSeq':
        quality = f'--2colour {qual_th}'

    elif _sequencer == 'NextSeq500':
        quality = f'--2colour {qual_th}'
    
    elif _sequencer == 'NovaSeqX':
        quality = f'--2colour {qual_th}'

    return quality

# --------- Rules ----------

rule trim_galore:
    input:
        fq_1 = 'resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz',
        fq_2 = 'resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz'
    output:
        out_dir = directory('results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        fastqc_dir = directory('results/quality_control/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        report1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz_trimming_report.txt',
        report2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz_trimming_report.txt',
        fq1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_2.fq.gz'
    threads: 20
    resources:
        mem_mb = 7000
    conda:
        '../../envs/pre_processing/trim-galore.yaml'
    log: 
        'logs/trim_galore/{sample_type}_{treatment}_Bio-rep_{bio_rep}.log'
    params:
        basename = '{sample_type}_{treatment}_Bio-rep_{bio_rep}',
        quality = _params_trim_galore,
        clip_5 = config['TRIM_GALORE']['CLIP_5'],
        clip_3 = config['TRIM_GALORE']['CLIP_3'],
        extra = config['TRIM_GALORE']['EXTRA']
    shell:
        """
            mkdir -p {output.fastqc_dir} && \
            trim_galore {params.quality} {params.extra} -j {threads} \
            --basename {params.basename} \
            --clip_R1 {params.clip_5} --clip_R2 {params.clip_5} \
            --gzip --paired --fastqc_args "--outdir {output.fastqc_dir}"\
            -o {output.out_dir} {input} 2> {log}
        """

rule unzip_trim:
    input:
        fq1_gz = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq2_gz = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_2.fq.gz'
    output:
        fq1 = temp('results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_1.fq'),
        fq2 = temp('results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_2.fq')
    shell:
        """
             gzip -dc {input.fq1_gz} > {output.fq1} &&
             gzip -dc {input.fq2_gz} > {output.fq2}
        """