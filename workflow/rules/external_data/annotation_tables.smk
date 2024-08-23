def _params_get_ensembl(wildcards):
  
    if GENCODE_RELEASE.startswith('M'):
        organism = 'mouse'
    else:
        organism = 'human'
    return organism


rule gft_transcriptBED:
    input:
        'resources/external/merged_ERCC_gencode_{release}/merged_ERCC_{genome}.annotation.gtf.gz'
    output:
        'resources/external/merged_ERCC_gencode_{release}/merged_ERCC_{genome}.transcripts_info.tsv'
    log:
        'logs/external_data/merged_ERCC_gencode_{release}_{genome}transcripts_info.log'
    shell:
        #   select colums : seqname start end strand attribute
        #   separate attributes
        #   add attribute as columns and select attributes to print on bed: id gene_name gene_source gene_biotype
        """
            zgrep -P "\ttranscript\t" {input} | cut -f1,4,5,7,9 | \
            sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
            awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $1,$2,$3,$8,$4,$16,$6,$12,$14}}' | \
            sort -k1,1 -k2,2n  > {output} && \
            sed  -i '1i\chromosome\tstart\tend\ttranscript_id\tstrand\ttranscript_name\tgene_id\tgene_name\tgene_type' {output} 
        """

# chr1    HAVANA  transcript      134487909       134563023       .       +       .       
# gene_id "ENSMUSG00000042207.18"; # 5 -6
# transcript_id "ENSMUST00000047714.14"; 
# gene_type "protein_coding"; 
# gene_name "Kdm5b"; # 11-12
# transcript_type "protein_coding"; 
# transcript_name "Kdm5b-201"; 
# level 2; 
# protein_id "ENSMUSP00000038138.8"; 
# transcript_support_level "1"; 
# mgi_id "MGI:1922855"; 
# tag "basic"; 
# tag "Ensembl_canonical"; 
# tag "appris_principal_1"; 
# tag "CCDS"; 
# ccdsid "CCDS35716.1"; 
# havana_gene "OTTMUSG00000017092.2"; 
# havana_transcript "OTTMUST00000041390.2";

rule gft_geneBED:
    input:
        'resources/external/merged_ERCC_gencode_{release}/merged_ERCC_{genome}.annotation.gtf.gz'
    output:
        'resources/external/merged_ERCC_gencode_{release}/merged_ERCC_{genome}.geneset.tsv'
    shell:
        #   select colums : seqname start end strand attribute
        #   separate attributes
        #   add attribute as columns and select attributes to print on bed: id gene_name gene_source gene_biotype
        """
            zgrep -P "\tgene\t" {input} | cut -f1,4,5,7,9 | \
            sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | \
            awk -F $'\t' 'BEGIN {{ OFS=FS }} {{ print $1,$2,$3,$6,$4,$3-$2+1,$10,$8 }}' | \
            sort -k1,1 -k2,2n > {output} && \
            sed  -i '1i\chromosome\tstart\tend\tgene_id\tstrand\tgene_length\tgene_name\tgene_type' {output} 
        """