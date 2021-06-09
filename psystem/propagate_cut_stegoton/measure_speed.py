import numpy as np

def get_error(data):
    # left and right indeces
    iL0 = int(3*len(data)/5.0)
    iR0 = int(5*len(data)/5.0)
    # get the time that minimizes the L1-error
    t = data[iL0:iR0,0]
    e = data[iL0:iR0,2]
    amin = e.argmin()
    print 'Time, error, speed: ', t[amin], e[amin], 20.0/t[amin]
#

data = np.genfromtxt('file.csv',delimiter=' ')
get_error(data)
