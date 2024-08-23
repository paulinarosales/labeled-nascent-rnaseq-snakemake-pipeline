rule samtools_stats:
    input:
        filteredBAM = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam',
        filteredBAI = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam.bai'
    output:
        'results/quality_control/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.stats'
    log: 
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_stats.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '20G',
        time = '1-00:00:00'
    shell:
        """
            samtools stats {input.filteredBAM} > {output} 2> {log}
        """
        
rule samtools_idxstats:
    input:
        filteredBAM = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam',
        filteredBAI = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam.bai'
    output:
        'results/quality_control/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.idxstats'
    log: 
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_idxstats.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '20G',
        time = '1-00:00:00'
    shell:
        """
            samtools idxstats {input.filteredBAM} > {output} 2> {log}
        """


rule samtools_flagstat:
    input:
        filteredBAM = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam',
        filteredBAI = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam.bai'
    output:
        'results/quality_control/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.flagstat'
    log: 
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_flagstat.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '20G',
        time = '1-00:00:00'
    shell:
        """
            samtools flagstat {input.filteredBAM} > {output} 2> {log}
        """

rule samtools_coverage:
    input:
        filteredBAM = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam',
        filteredBAI = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam.bai'
    output:
        'results/quality_control/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.coverage'
    log: 
        'logs/samtools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_coverage.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    threads: 24
    resources:
        mem = '20G',
        time = '1-00:00:00'
    shell:
        """
            samtools coverage -o {output} {input.filteredBAM}  2> {log}
        """