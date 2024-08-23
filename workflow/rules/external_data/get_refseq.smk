# def _params_get_seq(wildcards):
  
#     if GENCODE_RELEASE.startswith('M'):
#         organism = 'mouse'
#     else:
#         organism = 'human'
#     return organism
#     # return dict(organism=organism, link_handle=link_handle) # use params[0][organism] o shell



# rule get_genome:
#     output:
#         'resources/external/gencode_{release}/{genome}.genome.fa.gz'
#     params:
#         organism = _params_get_seq
#     log:
#         'logs/external_data/get_gencode_{release}_{genome}_genome.log'
#     threads: 1
#     resources:
#         http = 1
#     shell:
#         """
#             wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params.organism}/release_{wildcards.release}/{wildcards.genome}.primary_assembly.genome.fa.gz -O  {output} 2> {log}
#         """


# rule unzip_genome:
#     input:
#         rules.get_genome.output
#     output:
#         'resources/external/gencode_{release}/{genome}.genome.fa'
#     shell:
#         """
#              gzip -dc {input} > {output}
#         """

# rule get_annotation:
#     output:
#         'resources/external/gencode_{release}/{genome}.annotation.gtf.gz'
#     params:
#         organism = _params_get_seq
#     log:
#         'logs/external_data/get_gencode_{release}_{genome}_annotation.log'
#     threads: 1
#     resources:
#         http = 1
#     shell:
#         """
#             wget --quiet http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{params.organism}/release_{wildcards.release}/gencode.v{wildcards.release}.primary_assembly.annotation.gtf.gz -O  {output} 2> {log}
        # """


rule unzip_gft:
    input:
        rules.get_annotation.output
    output:
        'resources/external/merged_ERCC_gencode_{release}/{genome}.annotation.gtf'
    shell:
        """
             gzip -dc {input} > {output}
        """