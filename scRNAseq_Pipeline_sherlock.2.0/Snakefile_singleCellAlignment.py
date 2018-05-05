###############################################
# Snakemake rules associated with alignment
# of single cell RNAseq data and marking duplicates
# this file must be included into another 
# Snakemake file
# This file also deals with organizing aligned data
# In order to generate counts and QC data
# Update: added shell.prefix("set -euo pipefail;") to solve problems involved in running multi-line command in shell 
#### Update for Sherlock
# htseq-counting
# shell("source activate py2env")
# shell("""date                                                                                                                                                                              
# htseq-count -m intersection-strict  -s no -r name -f bam -o {wildcards.subsample}/htseq.sam.out {input_on_scratch[0]} {params.gtf_file} > {output_on_scratch}                            
# date""")
###############################################

# rules
shell.prefix("set -euo pipefail;")

rule STAR_MarkDuplicate:
  input:
    "{subsample}/P1.{subsample}.fastq",
    "{subsample}/P2.{subsample}.fastq"
  output:
    "{subsample}/STAR.{subsample}.sortedByCoord.bam",
    "{subsample}/STAR.{subsample}.sortedByCoord.bam.bai",
    "{subsample}/dedup.{subsample}.sortedByCoord.bam",
    "{subsample}/dedup.{subsample}.sortedByCoord.bam.bai",
    "{subsample}/picard.{subsample}.metric_file.txt",
    "{subsample}/dedup.{subsample}.sortedByName.bam"
  params:
    name="STAR_MarkDuplicate",
    partition=parameters.ix['STAR_partition','entry'],
    mem=parameters.ix['STAR_memory','entry'],
    qos="normal",
    time="05:00:00",
    STAR_genome_dir=parameters.ix['STAR_genome_directory','entry']
  threads: int(parameters.ix['STAR_thread','entry'])
  version: "1.0"
  run:
    # Managing files and obtain scratch location
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    # perform file size check
    file_size_check = 1
    for f in input:
      if os.path.getsize(f) == 0:
        file_size_check = 0
    # if file size check passes then proceed
    if file_size_check:
      # check fastq files have the same ids
      assert(check_fastq_ids(input[0],input[1])),"Input fastq files have different ids."
      cp_to_scratch(input, scratch)
      
      # STAR
      shell("bash {code_dir}/snakehelper_STARalign.sh {input_on_scratch} {params.STAR_genome_dir} \
        {threads} STAR_output_{wildcards.subsample} {scratch}")
      
      print(scratch+'/STAR_output_'+wildcards.subsample+'/Aligned.sortedByCoord.out.bam')
      shell("rsync -avrP {scratch}/STAR_output_{wildcards.subsample} {wildcards.subsample}/ &&\
        mv {scratch}/STAR_output_{wildcards.subsample}/Aligned.sortedByCoord.out.bam {output_on_scratch[0]} &&\
        mv {wildcards.subsample}/STAR_output_{wildcards.subsample}/Aligned.sortedByCoord.out.bam {output[0]} &&\
        mv {wildcards.subsample}/STAR_output_{wildcards.subsample}/Aligned.sortedByCoord.out.bam.bai {output[1]}")
      # call Piccard tools mark duplicated
      sortedByName_prefix = output_on_scratch[5][0:-4]  # Gets rid of .bam
      print('Picard Mark Duplicates Prefix = '+sortedByName_prefix)
      # print(output[0])
      shell("bash {code_dir}/snakehelper_markduplicates.sh {output_on_scratch[0]} {output_on_scratch[2]} \
        {output_on_scratch[5]} {tool_dir} {output_on_scratch[4]} {scratch}")
      cp_from_scratch(output[2:], scratch)
    else:
      print('Input files have size 0')

rule mpileup_alignment_vcf:
    """ Generates VCF from sorted STAR bam file """
    input:
        input_bam= "{subsample}/STAR.{subsample}.sortedByCoord.bam"
    output: '{subsample}/var.{subsample}.vcf'
    threads: 2
    params:
        name="mpile_vcf",
        reference_fasta=parameters.ix['reference_fasta','entry'],     
        partition="quake,owners",
        mem="10600"
    run:
        scratch = os.environ["L_SCRATCH_JOB"]
        input_on_scratch = names_on_scratch(input, scratch)
        output_on_scratch = names_on_scratch(output, scratch)
        cp_to_scratch(input, scratch)
        shell("echo {input_on_scratch} &&"
              "samtools mpileup --uncompressed "
              "--fasta-ref {params.reference_fasta} {input_on_scratch}| "
              "bcftools call -m --variants-only --output-type v -o {output_on_scratch}")
        cp_from_scratch(output, scratch)

rule htseq_rarefaction:
  # It is important that the input here is sorted by Name
  input: "{subsample}/dedup.{subsample}.sortedByName.bam"
  output: "{subsample}/htseq_downsample.{subsample}.{downsample}.out"
  params:
    name="htseq_downsample_and_count",
    partition="quake,owners",
    qos="normal",
    time="02:00:00",
    mem="10000",
    max_reads=parameters.ix['Down_Sample_Read_Number','entry'],
    gtf_file=parameters.ix['GTF_reference_annotation','entry']
  threads: 2
  version: "1.0"
  run:
    # manage input
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Compute subsample ratio a number between 0 and 1
    subsample_ratio = "%0.3f" % (float(wildcards.downsample)/float(params.max_reads))
    # run process
    shell("module load use.singlecell &&\
      module load python/2.7.9 &&\
      module load samtools/1.1 &&\
      date &&\
      echo {subsample_ratio} &&\
      samtools view -s {subsample_ratio} {input_on_scratch[0]} | htseq-count -m intersection-strict -s no -r name - {params.gtf_file} > {output_on_scratch} &&\
      date")
    cp_from_scratch(output, scratch)


rule htseq_counting:
  input: "{subsample}/dedup.{subsample}.sortedByName.bam"
  output: "{subsample}/htseq.{subsample}.sortedByName.out"
  params:
    name="htseq_counting",
    partition="quake,owners",
    qos="normal",
    time="02:00:00",
    mem="10000",
    gtf_file=parameters.ix['GTF_reference_annotation','entry']
  threads: 2
  version: "1.0"
  run:
    # manage input
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # run process
    #shell("source activate py2env")
    shell("source activate py2env &&\
      date &&\
      htseq-count -m intersection-strict -s yes -r name -f bam -i gene_name -o {wildcards.subsample}/htseq.sam.out --secondary-alignments ignore --supplementary-alignments ignore {input_on_scratch[0]} {params.gtf_file} > {output_on_scratch} &&\
      date")
    cp_from_scratch(output, scratch)


rule cufflinks:
  input: 
    "{subsample}/dedup.{subsample}.sortedByCoord.bam",
    "{subsample}/dedup.{subsample}.sortedByCoord.bam.bai"
  output: "{subsample}/cufflinks.{subsample}.genes_fpkm_tracking.out"
  params:
    name="cufflinks",
    partition="quake,owners",
    qos="normal",
    time="02:00:00",
    mem="10000",
    gtf_file=parameters.ix['GTF_reference_annotation','entry']
  threads: 2
  version: "1.0"
  run: 
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output,scratch)
    cp_to_scratch(input, scratch)
    # Need to transfer the bai file to scratch as well
    shell("bash {code_dir}/snakehelper_cufflinks.sh {input_on_scratch[0]} {output_on_scratch} {threads} \
      {params.gtf_file} {scratch}")
    cp_from_scratch(output,scratch)


rule run_rseqc:
  input: 
    "{subsample}/dedup.{subsample}.sortedByCoord.bam",
    "{subsample}/dedup.{subsample}.sortedByCoord.bam.bai"
  output: 
    "{subsample}/rseqc.{subsample}.read_distribution.txt",
    "{subsample}/rseqc.{subsample}.FPKM.xls",
    "{subsample}/rseqc.{subsample}.fragment_size.txt",
    "{subsample}/rseqc.{subsample}.geneBodyCoverage.txt", # must use geneBodyCoverage
    "{subsample}/rseqc.{subsample}.junction.xls" # must use junction
  params:
    name="run_rseqc",
    partition="quake,owners",
    qos="normal",
    time="02:00:00",
    mem="58000",
    bed_file=parameters.ix['BED_reference_annotation','entry']
  threads: 11
  version: "1.0"
  run:
    # Manage input
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    # Need to transfer bai file to scratch
    cp_to_scratch(input, scratch)
    # Create prefix
    FPKM_prefix = '.'.join(output_on_scratch[1].split('.')[0:-2])
    geneBody_prefix = '.'.join(output_on_scratch[3].split('.')[0:-2])
    junction_prefix = '.'.join(output_on_scratch[4].split('.')[0:-2])
    shell("bash {code_dir}/snakehelper_rseqcProcessing.sh {input_on_scratch[0]} {params.bed_file} \
      {output_on_scratch[0]} {FPKM_prefix} {output_on_scratch[2]} {geneBody_prefix} {junction_prefix} {scratch}")
    shell("cp {scratch}/rseqc.{wildcards.subsample}.* {wildcards.subsample}/")


