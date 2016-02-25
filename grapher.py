#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This routine reads in pickles created by grabber.py, rotates the data into
# the coordinates we want, and plots it. 

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

try:
  import cPickle as pickle
except ImportError:
  import pickle
import matplotlib.pyplot as plt
import numpy as np
import os

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# Load all of the pickles in a directory into a dictionary. 
def load(datadir):
  data = {}
  for pklname in os.listdir(datadir):
    with open(datadir + pklname, 'rb') as handle:
      data[ pklname[:-4] ] = pickle.load(handle)
  return data

# Based on the pickle directory name, return the probe name and the start and
# end time of the event, in seconds from midnight. Note that all events are ten
# minutes long by construction. 
def ptt(pkldir):
  probe, date, time = pkldir.split('_')
  hh, mm, ss = int( time[0:2] ), int( time[2:4] ), int( time[4:6] )
  t0 = ss + 60*mm + 3600*hh
  return probe, t0, t0 + 600

# Take a rolling average of the magnetic field to estimate the background
# field. The average gets truncated at the edges -- notably, none of Lei's
# events happen within ten minutes of midnight. 
def getbg(t, B, tavg=600):
  B0 = np.empty(B.shape)
  for i in range( B.shape[1] ):
    # Find indeces five minutes in the past and five minutes in the future. 
    ipas = np.argmax( t > t[i] - tavg/2 )
    ifut = np.argmax( t > t[i] + tavg/2 )
    # If we're looking for something past the edge of the array, argmax will
    # return 0. That a problem for the upper bound. 
    ifut = ifut if ifut>0 else B.shape[1]
    # Take the average. 
    B0[:, i] = np.average(B[:, ipas:ifut], axis=1)
  return B0

# Dot product of a pair of vectors or a pair of arrays of vectors. 
def dot(v0, v1, axis=0):
  return np.sum(v0*v1, axis=axis)

# Scale a vector, or an array of vectors, to unit vectors. 
def unit(v):
  return v/np.sqrt( dot(v, v) )

# The background magnetic field defines zhat. The azimuthal direction yhat is
# normal to the plane defined by zhat and rhat (as long as zhat and rhat are
# not parallel). The "radial" direction xhat is normal to zhat and yhat. 
def unitvecs(X, B0):
  rhat, zhat = unit(X), unit(B0)
  yhat = unit( np.cross(zhat, rhat, axis=0) )
  return np.cross(yhat, zhat, axis=0), yhat, zhat

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  # Location of the output directories. 
  outdir = '/media/My Passport/RBSP/pickles/'

  # Each event has its own directory full of pickles. 
  for pkldir in os.listdir(outdir):

    print pkldir + '/'

    # Load all of the pickles in this directory. Offset the time array to start
    # at midnight (rather than in 1970). 
    data = load(outdir + pkldir + '/')
    data['time'] = data['time'] - data['time'][0]

    for key, arr in data.items():
      print '\t' + key, arr.shape

    # Get the probe name and the timestamp from the directory name. From the
    # timestamps, find the data indeces corresponding to the start and end of
    # the event. 
    probe, t0, t1 = ptt(pkldir)
    i0, i1 = np.argmax(data['time'] > t0), np.argmax(data['time'] > t1)

    # Slice the event out of the full day of data. Also compute a ten-minute
    # rolling average of the magnetic field to estimate the background field. 
    t     = data['time'][i0:i1]
    BGSE  = data['bgse'][:, i0:i1]
    B0GSE = getbg( data['time'], data['bgse'] )[:, i0:i1]
    XGSE  = data['xgse'][:, i0:i1]

    # Compute the dipole coordinate directions. The zhat unit vector lines up
    # with the background magnetic field, yhat is azimuthally eastward, and
    # xhat completes the orthonormal coordinate system. 
    xhat, yhat, zhat = unitvecs(XGSE, B0GSE)

    # Compute the electric and magnetic field components in dipole coordinates.
    # The background magnetic field -- which lies in the zhat direction by
    # construction -- is kept separate. 
    B0 = dot(zhat, B0GSE)
    Bx = dot(xhat, BGSE - B0GSE)
    By = dot(yhat, BGSE - B0GSE)
    Bz = dot(zhat, BGSE - B0GSE)

    plt.plot(t, Bx, label='Bx')
    plt.plot(t, By, label='By')
    plt.plot(t, Bz, label='Bz')

    plt.legend()
    plt.show()




#    By = np.sum(yhat*B, axis=0)
#    Bz = np.sum(zhat*B, axis=0)
#    Bz0 = np.sum(zhat*B0, axis=0)



    print 'xhat'
    print xhat[:, 0]
    print np.sum( xhat[:, 0]**2 )

    print 'yhat'
    print yhat[:, 0]
    print np.sum( yhat[:, 0]**2 )

    print 'zhat'
    print zhat[:, 0]
    print np.sum( zhat[:, 0]**2 )

    print 'xhat dot yhat = ', np.sum( xhat[:, 0]*yhat[:, 0] )
    print 'yhat dot zhat = ', np.sum( yhat[:, 0]*zhat[:, 0] )
    print 'zhat dot xhat = ', np.sum( zhat[:, 0]*xhat[:, 0] )



    return



#    # Do a linear fit to get the background magnetic field. 
#    B0 = np.zeros( B.shape )
#    for i in range(3):
#      slope, inter = np.polyfit(t, B[i, :], 1)
#      B0[i, :] = inter + slope*t




    return

    lshell = data['lshell'][i0:i1]
    mlt = data['mlt'][i0:i1]
    mlat = data['mlat'][i0:i1]

  return




# Started on an event class... I think that's overkill. We're not going to be
# doing any serious computation here. 
'''
class eventObj:

  def __init__(self, datadir):

    # The data directory gives which probe this is, as well as the timestamp. 
    self.probe, self.date, self.time = datadir.split('/')[-1].split('_')

    # Read in the pickles. This is the whole day of data, in GSE coordinates. 
    self.load(datadir)

    # Use a ten minute running average to get the background magnetic field. 
    # That's the zhat direction. 
    self.setzhat()

    pickles = self.load(datadir)


  # Based on the pickle directory name, return the probe name and the start and
  # end time of the event, in seconds from midnight. Note that all events are
  # ten minutes long by construction. 
  def ptt(pkldir):
    probe, date, time = pkldir.split('_')
    hh, mm, ss = int( time[0:2] ), int( time[2:4] ), int( time[4:6] )
    t0 = ss + 60*mm + 3600*hh
    return probe, t0, t0 + 600

  def load(self, datadir):
    data = {}
    for pklname in os.listdir(datadir):
      with open(datadir + pklname, 'rb') as handle:
        data[ pklname[:-4] ] = pickle.load(handle)
'''

'''
# Make sure the given variable is a 3-by-N array.  
def is3d(x):
  return isinstance(x, np.ndarray) and len(x.shape)==2 and x[0]==3
'''

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


