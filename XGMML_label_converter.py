import re
import sys

XGMMLFile = open(sys.argv[1], 'r')

XGMML = XGMMLFile.readlines()

print(XGMML[0])

mapTableFile = open(sys.argv[2], 'r')

mapTable = mapTableFile.readlines()

mapTable = [ re.sub('\n','',line).split('\t') for line in mapTable]
mapTable.pop(0)

i = 0

for line in XGMML:
	if re.search('[<]node',line):
		label = re.sub('"','',re.search('["][^"]*["]',line).group(0))
		for protein in mapTable:
			if re.match(label,protein[0]):
				XGMML[i+1] = re.sub(label,protein[1],XGMML[i+1])
	i = i+1

newXGMMFile = open(sys.argv[3], 'w')

newXGMMFile.write(''.join(XGMML))
				


