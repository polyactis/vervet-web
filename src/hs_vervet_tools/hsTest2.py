import csv
import numpy as np

fname='writetest.txt'
fname2='/home/feilchenfeldt/writetest.txt'

snppos=np.array(['a','b','c'])
data=np.array([[1,1,0,2],[0,2,1,1],[1,0,0,0]])
ref=np.array(['A','G','T'])
alt=np.array(['C','A','A'])
#foo = np.asanyarray(range(100))
#mask = (foo < 40) | (foo > 60)

#print mask

plist=np.tile('p',(len(ref),))

print ref.shape,plist.shape

print np.core.defchararray.add(ref,plist)



#a=np.tile('bla',(3,1))
#zygotemat=np.random.rand(data.shape[0],2*data.shape[1])

#print np.column_stack((snppos,data,snppos,a))