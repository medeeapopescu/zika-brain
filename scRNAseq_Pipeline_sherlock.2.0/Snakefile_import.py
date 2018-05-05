##########################################################
# Rules for importing, trimming, and clustering fastq
#### Update for sherlock 2017.03.16
# shell("source activate py3env")
# shell("fastx_trimmer -Q33 -l {params.trim_to_read_length} -i {input[0]} -o {input_on_scratch[0]}")
# shell("fastx_trimmer -Q33 -l {params.trim_to_read_length} -i {input[1]} -o {input_on_scratch[1]}")
####
##########################################################

#####################################
# rules
#####################################

# use an input file to specify and cat the fastq files into local10G
rule concatenate:
  # for each subsample, find fastq files from all sequencing runs and combine them.
  # need to be able to handle certain subsamples appearing on in some sequencing runs but not others
  # input: No Input
  # temp removes the files once it has been used
  output: temp("{subsample}/read1.{subsample}.fastq"), temp("{subsample}/read2.{subsample}.fastq")
  params: 
    name="combine_fastq", 
    partition="quake,owners",
    qos="normal",
    time="02:00:00",
    mem="5000" # Don't change this
  threads: 1
  version: "1.0"
  run:
    scratch = os.environ["L_SCRATCH_JOB"]
    output_on_scratch = names_on_scratch(output, scratch)
    print(scratch)
    s = sample_table.ix[wildcards.subsample,:] 
    # if only one row, automatically become a column so need to handle that
    if len(s.shape) == 1:
      print('Sample '+wildcards.subsample+' exists in only one sequencing run.')
      fastq_dir = s.ix["biosample"] + "/" + s.ix["sequencingrun"] + "/" + s.ix["subsamplename"]
      shell("bash {code_dir}/snakehelper_combine_fastq.sh {fastq_dir} {scratch} {output_on_scratch[0]} {output_on_scratch[1]}")
    else:
      numseqruns = len(s.index)
      print('Sample '+wildcards.subsample+' exists in '+str(numseqruns)+' sequencing runs.')
      for i in range(numseqruns):
        fastq_dir = s.ix[i,"biosample"] + "/" + s.ix[i,"sequencingrun"] + "/" + s.ix[i,"subsamplename"]
        shell("bash {code_dir}/snakehelper_combine_fastq.sh {fastq_dir} {scratch} {output_on_scratch[0]} {output_on_scratch[1]}")
    assert(file_empty(output_on_scratch)),"Output fastq files are empty."
    assert(check_lines(output_on_scratch[0],output_on_scratch[1])),"Output fastq files have different number of lines."
    assert(check_fastq_ids(output_on_scratch[0],output_on_scratch[1])),"Output fastq ids are not in the same order."
    cp_from_scratch(output, scratch)


rule overrep_seq:
  # use "wild card" to specify for multiple parallel processes
  input: "{subsample}/{type}.{subsample}.fastq"
  output: "{subsample}/{type}.{subsample}.fastqc_results.txt"
  params: 
    name="fastqc", 
    partition="quake,owners",
    qos="normal",
    time="02:00:00",
    mem="3000" # Don't change this
  threads: 1
  version: "1.0"
  run: 
    # Managing files and obtain scratch location
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    output_on_scratch = names_on_scratch(output, scratch)
    cp_to_scratch(input, scratch)
    # Quality checks
    shell("bash {code_dir}/snakehelper_fastqc.sh {scratch} {input_on_scratch} {output_on_scratch} {tool_dir}")
    cp_from_scratch(output, scratch)


rule preprocess:
  input: 
    "{subsample}/read1.{subsample}.fastq",
    "{subsample}/read2.{subsample}.fastq"
  output:
    temp("{subsample}/P1.{subsample}.fastq"),
    temp("{subsample}/P2.{subsample}.fastq")
  params:
    name="reads_preprocess",
    partition="quake,owners",
    qos="normal",
    time="02:00:00",
    mem="10000",
    trim_to_read_length=str(parameters.ix['Desired_Read_Length','entry']),
    downsample_read_number=str(4*int(parameters.ix['Down_Sample_Read_Number','entry']))
  threads: 2
  version: "1.0"
  run:
    # organize input
    scratch = os.environ["L_SCRATCH_JOB"]
    input_on_scratch = names_on_scratch(input, scratch)
    # trim fastq to desired read length (ie. 75 bp)
    #shell("source activate py3env")
    #shell("fastx_trimmer -Q33 -l {params.trim_to_read_length} -i {input[0]} -o {input_on_scratch[0]}")
    #shell("fastx_trimmer -Q33 -l {params.trim_to_read_length} -i {input[1]} -o {input_on_scratch[1]}")
    #shell("bash snakehelper_trimmomatic.sh {input[0]} {input_on_scratch[0]} {input_on_scratch[1]}")
    # randomize fastq using a python script
    print (int(params.downsample_read_number))
    if int(params.downsample_read_number) <= 0:
      print('Do not downsample reads')
      shell("""mv {input[0]} {output[0]}
        mv {input[1]} {output[1]}""")
    else:
      shell("""source activate py2env
        python {code_dir}/randomizeFastq.py {input_on_scratch} {scratch}/temp""")
      # take the frist set number of lines
      print('Down Sample Read Number of Lines is: '+params.downsample_read_number)
      shell("head -n {params.downsample_read_number} {scratch}/temp1.fastq > {output[0]} &&\
        head -n {params.downsample_read_number} {scratch}/temp2.fastq > {output[1]} &&\
        wc -l {output[0]} &&\
        wc -l {output[1]}")


