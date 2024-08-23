# Pipeline for single nucelotide conversion detection on RNA-seq data

## Description
This pipeline uses the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow manager to process RNA-seq data from nucleotide converted datasets such as [SLAM-Seq](http://doi.org/10.1186/s12859-019-2849-7), [TimeLapse-seq](https://www.nature.com/articles/nmeth.4582) or [DART-seq](https://www.nature.com/articles/s41592-019-0570-0), based on the **detection of single nucleotide conversions** for targeted RNA identification. The code is based on the [HISAT-3N](https://daehwankimlab.github.io/hisat2/hisat-3n/) alignment tool within the Snakemake framework for automated and parallel running, and combined with different stategies for **nucelotide-conversion oriented analysis**.

## Settings requierments

### HISAT-3N 
Since the HISAT-3N tool is not yet available as a conda package, the user has to __download HISAT-3N scripts manually__ following the [HISAT-3N installation steps](https://daehwankimlab.github.io/hisat2/hisat-3n/) and indicate the location of the `hisat-3n` folder as a variable within the `config.yaml` (see below). 

As soon as a compatible installation is implemented we will add it to the pipeline.

### Snakemake 
To use this pipeline you would first need to create a Snakemake environment and and activate it. [See more.](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Getting started
1. Use this template to create your own project repo.

2. Create the `config/sample_manifest.tsv` using the following structure:

    1. **Identifier:** Unique sample ID used on the raw `.fastq.gz` sequencing files
    2. **Sample_type:** Tissue, cell type, etc.
    3. **Treatment**: Condition (treated/untreated, fed/starved, etc.)
    4. **Bio_rep:** Biological replicate number
    5. **Labeling_method:** Metabolic labeling method (SLAMseq or TimeLapse)
    6. **Labeling_time:** Incubation time for methabolic labeling
    7. **Target_genome:** ENCODE target genome for sequencing (*supported options:* M25, M32, 19, 38)
    8. **Sequencer**: Sequencer to define `—2colour` parameter for the trimming (*supported options:* HiSeq4000, NovaSeq, NextSeq500, NovaSeqX)
    9. **Seq_mode**: Library preparation strategy (*supported options:* mRNA, totalRNA)
    10. **Seq_length**: Sequencing length
    11. **Fastq_handle:** Particular handle useful for raw `.fastq.gz` files selection using the name from the sequencing facility (number, extension, etc.)
    12. **Fastq_lanes:** Number of `.fastq.gz` files lanes from the sequencing facility (used for `merge_fq_lanes` rule)
 
    An example of the sample manifest is available [here](./config/sample_manifest_example.tsv).

3. Copy your raw sequencing data to `resources/fastq_seq/raw`, the names of the files must correspond to the values indicated on the sample manifest columns, using the following structure:

    ```
        # Read 1 - Lane 1
        {Identifier}_S{Fastq_handle}_L001_R1_001.fastq.gz

        # Read 1c
        {Identifier}_S{Fastq_handle}_L002_R1_001.fastq.gz

        # Read 2 - Lane 1
        {Identifier}_S{Fastq_handle}_L001_R2_001.fastq.gz

        # Read 2 - Lane 1
        {Identifier}_S{Fastq_handle}_L002_R2_001.fastq.gz
    ```

    The names inside the `{}` correspond to the column name on the `config/sample_manifest.tsv`. The pipeline autoamtically merges lane files (if applicable) and renames the files to:

    ```
        # Read 1
        {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}_R1.fastq.gz

        # Read 2
        {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}_R2.fastq.gz
    ```
    Again, with the corresponding values on the `config/sample_manifest.tsv`. If yopu wish to change the name of the raw files yo can do so changing the values of ` TARGETS['raw_seq']` on the `Snakefile` and `merge_lanes.smk` rule.

<blockquote class="callout warning">
  <h4>⚠️ Warning</h4>
  <p>Changing the structure `{Sample_type}_{Treatment}_Bio-rep_{Bio_rep}` of filenames is not recommended and we don't ensure the proper functioning of the pipeline once it is done.</p>
</blockquote>

4. Check and change processing parameters on your `config/config.yaml` according to your analysis set-up.

5. Run the pipeline. 

    Your pipeline should look like this:

    ```    
    <YOUR PROJECT FOLDER>
        ├── config
        │   ├── config.yaml
        │   └── sample_manifest.tsv
        ├── profiles
        |   └── profile
        ├── resources
        |   └── fastq_seq
        |       └── raw
        |            ├── {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}_R1.fastq.gz
        |            └── {Sample_type}_{Treatment}_Bio-rep_{Bio_rep}_R2.fastq.gz
        └── workflow
    ```

    Use the following command on your project folder:

    ```
        snakemake --profile profiles/<YOUR_PROFILE> -j<N_JOBS>
    ```

    > For more options check the [--help](https://snakemake.readthedocs.io/en/stable/executing/cli.html)

    
## Important notes

- This pipeline only works with _paired-end_ sequencing reads.


- The default configurations for this workflow are suitable to run snakemake on the `Helmholtz-Munich HPC` and using `conda` to deal with software dependencies. If one wishes to run the pipeline on a different computing platform, the profiles need to be adapted accordingly.


- The pipeline only supports GENCODE genomes as reference, choose the right `GENCODE Realease` as value for the `Target_genome` entry in your `config/sample_manifest.tsv`. These are the equivalent realeases for GENCODE and UCSC databases:


    | Organism | GENCODE Realease | GENCODE Genome | UCSC Genome |
    | ----- | ---- | ----- | ---- |
    | Mouse | M25 | GRCm38 | mm10 |
    | Mouse | M32 (Latest) | GRCm38 | mm39 |
    | Human | 19 | GRCh37 | h19 |
    | Human | 38 | GRCh38 | h38 |


## Authors and acknowledgment

Paulina Rosales-Becerra, Kevin Brokers & Robert Schneider

## Contact
paulina.rosales@helmholtz-munich.de
