def _input_merge_fq(wildcards):
    sample_type = wildcards.sample_type
    treatment = wildcards.treatment
    bio_rep = wildcards.bio_rep
    
    identifier = SAMPLES.loc[(sample_type, treatment, bio_rep), 'Identifier'][0]
    handle = SAMPLES.loc[(sample_type, treatment, bio_rep), 'Fastq_handle'][0]
    lane_1 = int(SAMPLES.loc[(sample_type, treatment, bio_rep), 'Lanes_ID'][0])
    lane_2 = lane_1+1

    _n_lanes = SAMPLES.loc[(sample_type, treatment, bio_rep), 'Fastq_lanes'][0]
    

    if _n_lanes == 1:
        l_r1 = f'resources/fastq_seq/raw/{identifier}_S{handle}_L00{lane_1}_R1_001.fastq.gz'
        l_r2 = f'resources/fastq_seq/raw/{identifier}_S{handle}_L00{lane_1}_R2_001.fastq.gz'

    elif _n_lanes == 2:
        l_r1 = [f'resources/fastq_seq/raw/{identifier}_S{handle}_L00{lane_1}_R1_001.fastq.gz',  f'resources/fastq_seq/raw/{identifier}_S{handle}_L00{lane_2}_R1_001.fastq.gz']
        l_r2 = [f'resources/fastq_seq/raw/{identifier}_S{handle}_L00{lane_1}_R2_001.fastq.gz',  f'resources/fastq_seq/raw/{identifier}_S{handle}_L00{lane_2}_R2_001.fastq.gz']

    else: 
        print('Fastq_lanes argument on sample_manifest.tsv non-compatible.\n')

    return {'lanes_r1': l_r1, 'lanes_r2': l_r2}



rule merge_fq_lanes:
    input:
        unpack(_input_merge_fq)
    output:
        outdir = directory('resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}'),
        merged_fq1 = 'resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R1.fq.gz',
        merged_fq2 = 'resources/fastq_seq/merged/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_R2.fq.gz'
    threads: 6
    shell:
        """
            mkdir -p {output.outdir} && \
            cat {input.lanes_r1} > {output.merged_fq1} && \
            cat {input.lanes_r2} > {output.merged_fq2}
        """

