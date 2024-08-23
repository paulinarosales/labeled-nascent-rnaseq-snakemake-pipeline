rule samtools_sort_coord:
    input:
        'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sam'
    output:
        'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam'
    log:
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_sort_coord.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '120G',
        time = '1-00:00:00'
    shell:
        """
            samtools sort --threads {threads} -O bam -o {output} {input} 2> {log}
        """



rule samtools_index:
    input:
        'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam'
    output:
        'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam.bai'
    log: 
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_index.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '64G',
        time = '1-00:00:00',
        tmp_dir = 'tmp'
    shell:
        """
            samtools index {input} {output} 2> {log}
        """


rule samtools_sort_name:
    input:
        'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sam'
    output:
        temp('results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.sam')
    log:
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_sort_name.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '120G',
        time = '1-00:00:00',
        tmp_dir = 'tmp'
    shell:
        """
            samtools sort -n --threads {threads} {input} |\
            samtools fixmate --threads {threads} - {output} 2> {log}
        """


rule samtools_filter:
    input:
        'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.sam'
    output:
        temp('results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.sam')
    log:
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filter.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '20G',
        time = '1-00:00:00',
        tmp_dir = 'tmp'
    shell:
        """
            samtools view --threads {threads} -q 2 -h {input} |\
            awk '$1 ~ /^@/ {{print}}
                  ($2 == 147 || $2 == 99) || ($2 == 83 || $2 == 163) {{print}}' > {output} 2> {log}
        """

rule samtools_bam:
    input:
        'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.sam',
    output:
        'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sortn.filter.bam'
    log: 
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_filter_bam.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '64G',
        time = '1-00:00:00',
        tmp_dir = 'tmp'
    shell:
        """
            samtools sort --threads {threads} -n -O bam -o {output} {input} 2> {log}
        """