###############################################
# Snakemake rules associated with combining data
# from multiple single cells in order to be used
# for downstream processing.
# this file must be included into another 
# Snakemake file
###############################################

# functions

def combine_gene_count(input, output):
  """
  This function takes in a list of input htseq counts files
  and joins them. I also checks that all the gene names
  are correctly ordered and matched.
  You need to make sure you are keeping consistent sample and gene orders.
  I'm assuming that in each input name, the 3rd last string
  after splitting by '.' is the sample name (ILxxxx-N7xx-N5xx)
  """
  print(len(input))
  output_file = output[0]
  print(output_file)
  assert(file_empty(input)),"At least one of the input files are empty."
  samplename = [x.split('.')[-3] for x in input] # order same as input
  # find number of genes 
  with open(input[1],'r') as f:
    totalgenes = 0
    genenames = []
    for l in f:
      totalgenes += 1
      genenames.append(l.split()[0])
  genedict = dict(zip(genenames, range(len(genenames))))
  # Initialize combined array
  # GeneMat = np.zeros([totalgenes, len(input)])
  GeneMat = pd.DataFrame(index=genenames,columns=samplename)
  # Fill matrix
  for i in range(len(input)):
    # sample name should be the first part before slash
    with open(input[i],'r') as f:
      print(samplename[-1])
      for l in f:
        key = l.split()[0]
        GeneMat.ix[key, samplename[i]] = int(l.split()[1])
  # Generate output
  with open(output_file,'w') as f:
    a = f.write('Sample_names,')
    a = f.write(','.join(GeneMat.columns))
    a = f.write('\n')
    # now write data
    for key in GeneMat.index:
      a = f.write(key+',')
      a = f.write(','.join([str(x) for x in GeneMat.ix[key,:]]))
      a = f.write('\n')
  assert(file_empty(output)),"Output file is empty."
  print("Combining htseq files completed")


def organize_read_proportion(input,output):
  """
  This file tabulates how many reads are lost 
  in each step of the pipeline.
  This function keeps the same ordering of 
  sample names
  """
  input_file = input[0]
  output_file = output[0]
  print(input_file)
  print(output_file)
  # bring in input file
  expression = pd.read_table(input_file,header=0,index_col=0,delimiter=',')
  samplename = expression.columns.values
  tmp = expression.drop(['__no_feature','__ambiguous','__too_low_aQual','__not_aligned','__alignment_not_unique'])
  aligned_to_genes = []
  no_feature = expression.ix['__no_feature']
  ambiguous = expression.ix['__ambiguous']
  low_qual = expression.ix['__too_low_aQual']
  not_aligned = expression.ix['__not_aligned']
  not_unique = expression.ix['__alignment_not_unique']
  for i in range(len(samplename)):
    aligned_to_genes.append(str(sum(tmp.ix[:,i])))
  # Write to output file
  samplename = expression.columns
  with open(output_file,'w') as f:
    # write header
    a = f.write('Sample_names,')
    a = f.write(','.join(samplename))
    a = f.write('\n')
    # write not aligned
    t = [str(not_aligned[i]+no_feature[i]+low_qual[i]) for i in range(len(not_aligned))]
    a = f.write('Not_Aligned,')
    a = f.write(','.join(t))
    a = f.write('\n')
    # write not uniquely aligned
    t = [str(ambiguous[i]+not_unique[i]) for i in range(len(ambiguous))]
    a = f.write('Not_Unique,')
    a = f.write(','.join(t))
    a = f.write('\n')
    # write aligned
    a = f.write('Properly_Aligned,')
    a = f.write(','.join(aligned_to_genes))
    a = f.write('\n')
  assert(file_empty(output)),"Output file is empty."
  print("Combining read proportions files completed")


# rules

rule combine_gene_count:
  input: expand("{subsample}/htseq.{subsample}.sortedByName.out", subsample=subsampleIDs)
  output: "Combined_Analysis/expression_counts_{tag}.csv"
  params:
    name="combine_gene_count",
    partition="quake,owners",
    qos="normal",
    mem="15000",
    time="48:00:00"
  threads: 3
  version: "1.0"
  run:
    combine_gene_count(input, output)


rule organize_regular_read_distribution:
  input: "Combined_Analysis/expression_counts_{tag}.csv"
  output: "Combined_Analysis/read_proportion_{tag}.csv"
  params:
    name="organize_regular_read_distribution",
    partition="quake,owners",
    qos="normal",
    mem="5000",
    time="05:00:00"
  threads: 1
  version: "1.0"
  run:
    organize_read_proportion(input,output)


rule combine_genebody:
  input: expand("{subsample}/rseqc.{subsample}.geneBodyCoverage.txt", subsample=subsampleIDs)
  output: "Combined_Analysis/combined_gene_body_{tag}.csv"
  params:
    name="combine_genebody_results",
    partition="quake,owners",
    qos="normal",
    mem="10000"
  threads: 2
  version: "1.0"
  run:
    print(len(input))
    output_file = output[0]
    print(output_file)
    assert(file_empty(input)),"At least one of the input files are empty."
    samplename = subsampleIDs
    # find number of a genebody entries 
    with open(input[1],'r') as f:
      l = f.readline()
      percentile = l.split()[1:]
    # Initialize combined array
    genebody = pd.DataFrame(index=percentile,columns=samplename)
    # Fill matrix
    for i in range(len(input)):
      # sample name should be the first part before slash
      with open(input[i],'r') as f:
        header = f.readline()
        l = f.readline() # all the info is on the second line
        print(l)
        l = l.split()
        key = input[i].split('.')[-3]
        genebody.ix[:,key] = l[1:]
    # Generate output
    with open(output_file,'w') as f:
      a = f.write('Percentile,')
      a = f.write(','.join(genebody.columns))
      a = f.write('\n')
      # now write data
      for key in genebody.index:
        a = f.write(key+',')
        a = f.write(','.join([str(x) for x in genebody.ix[key,:]]))
        a = f.write('\n')
    assert(file_empty(output)),"Output file is empty."
    print("Combining genebody files completed")
    
