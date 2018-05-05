###############################################
# Snakemake rules associated with combining data
# from multiple single cells in order to be used
# for downstream processing.
# this file must be included into another 
# Snakemake file
# In particular for rarefaction curves
###############################################

# functions

def organize_multiple_rarefaction(input,output,downsampleList,geneDetectionThresh):
  """
  This function organizes read counts
  from multiple rarefaction down sampled
  expression count files.
  Using geneDetectionThresh to mean detected.
  """
  # downsampleList is a list of strings not ints
  print(len(input))
  print(input) # check if input is in the proper order
  output_file = output[0]
  print(output_file)
  assert(file_empty(input)),"At least one of the input files are empty."
  assert(len(input) == len(downsampleList)),"Number of input files and number downsample points are different."
  # Changing gene detection threshold to an int
  geneDetectionThresh = int(geneDetectionThresh)
  # Create an empty data frame
  tmp = pd.read_table(input[0],header=0,index_col=0,delimiter=',')
  samplename = tmp.columns
  rowlabels = ['down_sample_to_'+str(x) for x in downsampleList]
  rarefaction_data = pd.DataFrame(index=rowlabels, columns=samplename)
  input_order = [int(name.split('.')[-2]) for name in input]
  print(input_order)
  # Tabulate number of genes found
  for i in range(len(downsampleList)):
    # First find the correct input file to open up Combi.../rarefaction_expression_counts.xxxxxx.csv
    input_file_name = input[input_order.index(int(downsampleList[i]))]
    # Open up the input file
    expression = pd.read_table(input_file_name,header=0,index_col=0,delimiter=',')
    numrow = len(expression.index)
    print(str(numrow))
    numcol = len(expression.columns)
    # Check that the number of samples are correct
    assert(numcol == len(samplename)),"Sample names have different lengths."
    expression = expression.drop(['__no_feature','__ambiguous','__too_low_aQual','__not_aligned','__alignment_not_unique'])
    for k in samplename:
      # This statement works for pandas data frames but not for typical arrays (ints here)
      rarefaction_data.ix[rowlabels[i],k] = sum(expression.ix[:,k] > geneDetectionThresh)
  print(rarefaction_data)
  # Write to output file
  with open(output_file,'w') as f:
    a = f.write('Sample_names,')
    a = f.write(','.join(rarefaction_data.columns))
    a = f.write('\n')
    # now write data
    for key in rarefaction_data.index:
      a = f.write(key+',')
      a = f.write(','.join([str(x) for x in rarefaction_data.ix[key,:]]))
      a = f.write('\n')
  assert(file_empty(output)),"Output file is empty."
  print("Combining rarefaction files completed")
  

# rule

rule combine_rarefaction_gene_count:
  input: expand("{subsample}/htseq_downsample.{subsample}.{{downsample}}.out", subsample=subsampleIDs)
  output: "Combined_Analysis/rarefaction_expression_counts.{downsample}.csv"
  params:
    name="combine_rarefaction_gene_count",
    partition="general",
    mem="15000"
  threads: 3
  version: "1.0"
  run: 
    combine_gene_count(input, output)


rule organize_rarefaction:
  input: expand("Combined_Analysis/rarefaction_expression_counts.{downsample}.csv", downsample=downsampleList)
  output: "Combined_Analysis/rarefaction_gene_numbers.csv"
  params:
    name="organize_rarefaction",
    partition="general",
    mem="15000",
    geneDetectionThresh=parameters.ix['Gene_Detection_Threshold','entry']
  threads: 3
  version: "1.0"
  run: 
    organize_multiple_rarefaction(input,output,downsampleList,int(params.geneDetectionThresh))


rule organize_rarefaction_read_distribution:
#  input: expand("Combined_Analysis/rarefaction_expression_counts.{tag}.csv", tag=[downsampleList[-1]])
  input: "Combined_Analysis/rarefaction_expression_counts.{tag}.csv"
  output: "Combined_Analysis/rarefaction_read_proportion.{tag}.csv"
  params:
    name="organize_rarefaction_read_distribution",
    partition="general",
    mem="5000"
  threads: 1
  version: "1.0"
  run: 
    organize_read_proportion(input,output)



