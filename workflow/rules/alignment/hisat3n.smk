# --------- Parsing functions ----------
# -- Build --
def _input_refGenome(wildcards):
    return expand('resources/external/merged_ERCC_gencode_{release}/merged_ERCC_{genome}.genome.fa', release=GENCODE_RELEASE, genome=GENOME) # de-compressed


def _params_for_hisat3n(wildcards):
    # add --repeat-index if requested
    _repeat_index = config['HISAT3N']['REPEAT_INDEX']
    if _repeat_index == 'Yes':
        repeat_index = '--repeat-index'
    elif _repeat_index == 'No':
        repeat_index = ''

    # add options for low memory consumption if requested
    _resources_mode = config['HISAT3N']['RESOURCES_MODE']
    if _resources_mode == 'Default':
        resources = ''

    elif _resources_mode == 'Low':
        resources = '--noauto --bmaxdivn 8 --dcv 512'

    return f'{repeat_index} {resources}'

# -- Align --

def _input_for_hisat3n_align(wildcards):
    # make sure that hisat3n_build is executed first
    # INDEX_DIR = os.path.dirname(checkpoints.hisat3n_build.get(release=GENCODE_RELEASE, genome=GENOME).output[0])
    INDEX_DIR = os.path.dirname(f'resources/external/index/hisat3n/gencode_{GENCODE_RELEASE}/{GENOME}.genome')
    INDEX_FILES = []
    # print(INDEX_DIR)
    for root, directories, files in os.walk(INDEX_DIR):
        for file in files:
            # list files under hisat3n_build output directory
            INDEX_FILES.append(file)
    return expand('{index_dir}/{index_file}', index_dir=INDEX_DIR, index_file=INDEX_FILES)


def _index_basename_for_hisat3n(wildcards):
    return expand('resources/external/index/hisat3n/merged_ERCC_gencode_{release}/merged_ERCC_{genome}.genome', release=GENCODE_RELEASE, genome=GENOME)

def _params_for_hisat3n_align(wildcards):
    # add --repeat if requested
    _repeat_index = config['HISAT3N']['REPEAT_INDEX']
    if _repeat_index == 'Yes':
        repeat_index = '--repeat'

    elif _repeat_index == 'No':
        repeat_index = ''
    return repeat_index

# --------- Rules ----------

rule hisat3n_align:
    input:
        index = _input_for_hisat3n_align,
        fq1 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_1.fq.gz',
        fq2 = 'results/fastq/trimmed/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_val_2.fq.gz'
    output:
        alignSAM = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sam',
        unalignSAM_1 = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_unalign.1.sam',
        unalignSAM_2 = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_unalign.2.sam'
    log:
        'logs/hisat3n/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_align.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 32
    resources:
        mem = '300G',
        time = '1-00:00:00'
    params:
        command = f'{HISAT3N_BIN_PATH}/hisat-3n', # run HISAT-3N script from git clone folder
        index_basename = _index_basename_for_hisat3n,
        repeat = _params_for_hisat3n_align,
        base_change = config['HISAT3N']['BASE_CHANGE'],
        unalign_basename = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_unalign.%.sam',
        # strandness = config['HISAT3N']['STRANDNESS'],
        extra = config['HISAT3N']['ALIGN_EXTRA']
    shell:
        """
            {params.command} -p {threads} -x {params.index_basename}\
            -q -1 {input.fq1} -2 {input.fq2}\
            -S {output.alignSAM} --base-change {params.base_change}\
            {params.repeat} {params.extra}\
            --un-conc {params.unalign_basename}\
            --new-summary --summary-file {log} 
        """