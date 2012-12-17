'''
Created on Nov 7, 2012

@author: hannes.svardal
'''
import hsVCFExtractionTools

obj=hsVCFExtractionTools.hsVCFExtractionTools("TestAnalyses1")
print '1'

#desc_str='This time I append test project 1\n'
#obj.createProjectDescription(desc_str,'a')

uclaid_ls=obj.loadUCLA_id_ls()

print uclaid_ls



#(idx_ls,ucla_id_ls)=obj.getVCFInd(uclaidlist=['VWP00312','VSAB2023','VSAE3003','VWP00456'],genotypeMethodID=40,chromosome='Contig917')
#print idx_ls,ucla_id_ls

