import wxversion
#wxversion.select('2.8')

import matplotlib
matplotlib.use('Agg')
import pylab
import numpy as np
#matplotlib.use('Qt4Agg')


n = 5
X = np.random.normal(0,1,n)
Y = np.random.normal(0,1,n)

matplotlib.pyplot.scatter(X,Y)
matplotlib.pyplot.show()
#
#import numpy
#import pylab 
#t = numpy.arange(0.0, 1.0+0.01, 0.01)
#s = numpy.cos(2*2*numpy.pi*t)
#pylab.plot(t, s) 
#pylab.xlabel('time (s)')
#pylab.ylabel('voltage (mV)')
#pylab.title('About as simple as it gets, folks')
#pylab.grid(True)
#pylab.savefig('simple_plot')
# 
#pylab.show()
