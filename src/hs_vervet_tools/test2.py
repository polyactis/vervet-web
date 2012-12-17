#import os
#filename = os.path.join("/home/GMI", 'hannes.svardal/Documents')

#print filename
import hsVCFExtractionTools

obj=hsVCFExtractionTools.hsVCFExtractionTools("TestAnalyses1")
genotypeMethodID=40

(VCFfilename_ls,contig_ls)=obj.loadVCFfilename_ls(genotypeMethodID=genotypeMethodID)
print VCFfilename_ls
print contig_ls
