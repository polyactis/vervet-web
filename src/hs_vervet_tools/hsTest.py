import csv

fh = open( "/media/Data/Akademisches_data/vervetpopgen/DistanceMatrix_pylogeny/PairwiseDistance_Method27.2012.9.5T1552/aggregateDistForm1.csv" );

print range(1,len([1,2,3]))

x = []
for line in fh.readlines():
    y = [value for value in line.split()]
    x.append( y )

fh.close()


for i in range(len(x[0])):
        spl=x[0][i].split('_')
        if len(spl)>2:
            x[0][i]=spl[2]

# move header 1 step to right to empty left upper corner
x[0][:0]='-'

for i in range(1,len(x)):
        spl=x[i][0].split('_')
        if len(spl)>2:
            x[i][0]=spl[2]

fw = csv.writer(open( "/media/Data/Akademisches_data/vervetpopgen/DistanceMatrix_pylogeny/PairwiseDistance_Method27.2012.9.5T1552/aggregateDistForm2.csv",'w'),delimiter='\t');

for line in x:
    fw.writerow(line)


#print x
#print distline.split(line.split(Arr