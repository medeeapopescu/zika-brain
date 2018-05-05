"""
Combine htseq output files column wise
Checks if genes are the same


Edit History:
2014.11.19: Created
2014.12.07: Edited to include seedfile
"""

import os, glob, numpy, argparse

usage = "usage: python combine_gene_count.py exp_folder seedfile outputfile"

# check in put arguments
p = argparse.ArgumentParser(usage=usage)
p.add_argument(dest='exp_dir', action='store', type=str)
p.add_argument(dest='seedfile', action='store', type=str)
p.add_argument(dest='outfile', action='store', type=str)
a = p.parse_args()

# Change to folder with all sample folders
print a.exp_dir
os.chdir(a.exp_dir)

# Import folder list
folderfid = open(a.seedfile, 'r')
folderlist = []
for l in folderfid:
    folderlist.append(l.split()[0])
folderfid.close()
numfile = len(folderlist)

# find number of genes in folder list
os.chdir(folderlist[0] + '/CIRM_pipeline/htseq_output')
# os.chdir(folderlist[0] + '/htseq_output')
filename = glob.glob('*htseq.out')[0]
genefid = open(filename,'r')
totalgenes = 0
genenames = []
for l in genefid:
    totalgenes += 1
    genenames.append(l.split()[0])
genefid.close()
genedict = dict(zip(genenames, range(len(genenames))))
os.chdir('../../..')
# os.chdir('../..')

# initializae combined array
GeneMat = numpy.zeros([totalgenes, numfile])

# Fill matrix
samplename = []
for i in range(numfile):
    # os.chdir(folderlist[i] + '/htseq_output')
    os.chdir(folderlist[i] + '/CIRM_pipeline/htseq_output')
    filename = glob.glob('*htseq.out')[0]
    genefid = open(filename, 'r')
    samplename.append(filename.split('.')[0])
    print samplename[-1]
    for l in genefid:
        key = l.split()[0]
        GeneMat[genedict[key], i] = int(l.split()[1])
    genefid.close()
    # os.chdir('../..')
    os.chdir('../../..')

# Open and prepare output file
outfid = open(a.outfile,'w')
# print header
outfid.write('Sample_Names')
for name in samplename:
    outfid.write(',%s' %name)
outfid.write('\n')
# write genes
for i in range(totalgenes):
    outfid.write('%s' %genenames[i])
    temparr = GeneMat[i,:]
    for genecnt in temparr:
        outfid.write(',%1.0f' %genecnt)
    outfid.write('\n')
outfid.close()

print "Data Extraction Completed"



