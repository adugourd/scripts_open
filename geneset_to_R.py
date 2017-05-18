import re
import sys

geneSetsFile = open(sys.argv[1], 'r')

geneSets = geneSetsFile.readlines()

targetGeneSet = []

for geneSet in geneSets :
	geneSet = geneSet.split("\t")
	print(geneSet.pop(0))
	geneSetName = geneSet.pop(0)
	for gene in geneSet :
		targetGeneSet.append(geneSetName+"\t"+gene)

targetGeneSet = "\n".join(targetGeneSet)
targetGeneSet = re.sub("[\n]{2}","\n",targetGeneSet)
targetGeneSet = re.sub(">[\s]","",targetGeneSet)

targetGeneSetFile = open(re.sub("[.]","_target.",sys.argv[1]), 'w')
targetGeneSetFile.write(targetGeneSet)
