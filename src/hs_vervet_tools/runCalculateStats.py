import hsCalculateStatsForSubPop as hs
reload(hs)

test1=hs.CalculateStatsForSubPop([],**{'outputFname': '/home/GMI/hannes.svardal/Desktop/lab/vervetpopgen/analyses/testplink1.tsv', 'hostname': 'dl324b-1.cmb.usc.edu', 'dataDir': '/home/GMI/hannes.svardal/Desktop/lab/vervetpopgen/genotype_file', 'genotypeMethodID': '39', 'db_passwd': '1o9p2a4', 'db_user': 'hannes', 'statname': 'aaf', 'chromosome': 'Contig917'})

test1.loadMetadataMat("/home/GMI/hannes.svardal/Desktop/lab/vervetpopgen/AfricaExtractedData/AfriCaribMetaData.tsv")

idls=test1.selectUCLAids(test1.metadata)

(indx,ind_id_ls)=test1.getVCFInd(idls)

Fname="/home/GMI/hannes.svardal/Desktop/lab/media/Data/Akademisches_data/vervetpopgen/AfricaExtractedData/Method27ContigNames.tsv"
contiglist=open(Fname).readlines()
contiglist = map(lambda s: s.strip(), contiglist)

filenamels=test1.produceContigFileNameLS(contiglist)


datastr=test1.selectSubPopNoDB(indx,ind_id_ls,filenamels)



#datastruct=test1.selectSubPop(idls)

#datastruct2=test1.thinSNPs(datastruct,1000)

#print datastruct.data.shape

#print datastruct2.data.shape

#test1.producePLINK(datastruct2)

