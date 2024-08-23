# --------- Parsing functions ----------
def _input_refGenome(wildcards):
    return expand('resources/external/merged_ERCC_gencode_{release}/merged_ERCC_{genome}.genome.fa', release=GENCODE_RELEASE, genome=GENOME) # de-compressed

def _input_control_sortedSAMs(wildcards):
    sample_type = wildcards.sample_type
    treatment = wildcards.treatment
    bio_rep = wildcards.bio_rep

    sortedBAM = []
    sortedBAI = []

    if SAMPLES.loc[(sample_type, treatment, bio_rep), "Control"] == 1: 
        sortedBAM.append(f'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam')
        sortedBAIappend(f'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam.bai')
    return [] # de-compressed



# --------- Rules ----------
rule bcftools_snp:
    input:
        refGenome = _input_refGenome,
        sortBAM = 'results/sam_files/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.hisat3n_align.sort.bam'
    output:
        'results/snps/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.snp.vcf'
    log:
        'logs/bcftools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_snp_call.log'
    params:
        ploidy = config['BCFTOOLS']['PLOIDY'],
        min_qual = config['BCFTOOLS']['MIN_QUAL'],
        min_dp = config['BCFTOOLS']['MIN_DEPTH']
    threads: 24
    conda:
        '../../envs/alignment/hisat3n.yaml'
    shell:
        """
            bcftools mpileup --threads {threads} -f {input.refGenome} {input.sortBAM} |\
            bcftools call --threads {threads} --ploidy {params.ploidy} -mv |\
            bcftools filter --threads {threads} -i 'QUAL>{params.min_qual} && DP>{params.min_dp}' > {output} 2> {log}
        """

rule parse_snp_txt:
    input:
        'results/snps/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.snp.vcf'
    output:
        'results/snps/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}.snp.txt'
    log:
        'logs/bcftools/{sample_type}_{treatment}_Bio-rep_{bio_rep}/{sample_type}_{treatment}_Bio-rep_{bio_rep}_snp_parse.log'
    conda:
        '../../envs/alignment/hisat3n.yaml'
    shell:
        """
            awk '$1 !~ /^#/ && length($4) == 1 {{if (length($5) == 1) {{print $4":"$5":"$1":"$2}}
                                            else if (length($5) == 3 && $5 ~ /,/) {{split($5, mut, ",")
                                                                                   print $4":"mut[1]":"$1":"$2
                                                                                   print $4":"mut[2]":"$1":"$2}}
                                            }}' {input} |\
            sort |\
            uniq > {output} 2> {log}
        """