import csv

Fname="/media/Data/Akademisches_data/vervetpopgen/DistanceMatrix_pylogeny/PairwiseDistance_Method40.2012.Oct.9T155710/folderReduce/aggregateDistanceMatrix.tsv"
distMat=open(Fname).readlines()
distMat2 = map(lambda s: s.strip(), distMat)

idls=map(lambda s: s.split('_')[2] if len(s.split('_'))>1 else s.split('_')[0],distMat2[-132].split())

print idls

distMat3=distMat2[-131:]

distMat4=[]
for line in distMat3:
	distMat4.append(line.split()[1:])
	
outFnameMat='/media/Data/Akademisches_data/vervetpopgen/DistanceMatrix_pylogeny/PairwiseDistance_Method40.2012.Oct.9T155710/folderReduce/aggregateDistanceMatrix_hs_format.tsv'
outFnameID='/media/Data/Akademisches_data/vervetpopgen/DistanceMatrix_pylogeny/PairwiseDistance_Method40.2012.Oct.9T155710/folderReduce/aggregateDistanceMatrix_ID_ls_hs_format.tsv'
writerMat = csv.writer(open(outFnameMat, 'w'), delimiter='\t')
writerID =csv.writer(open(outFnameID, 'w'), delimiter='\t')
#writerID = open(outFnameID, 'w')

writerMat.writerows(distMat4)
writerID.writerow(idls)
#for id in idls:
#	writerID.write(id)
#	writerID.write('\n')
