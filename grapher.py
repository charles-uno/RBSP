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


from plotmod import *

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
def pdtt(pkldir):
  probe, date, time = pkldir.split('_')
  hh, mm, ss = int( time[0:2] ), int( time[2:4] ), int( time[4:6] )
  t0 = ss + 60*mm + 3600*hh
  return probe, date, t0, t0 + 600

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

# Put the slashes back into a date string. 
def calendar(date):
  return date[0:4] + '/' + date[4:6] + '/' + date[6:8]

# Compute clock time from a count of seconds from midnight. 
def clock(sfm):
  hh = sfm/3600
  mm = (sfm%3600)/60
  ss = sfm%60
  return znt(hh, 2) + ':' + znt(mm, 2) + ':' + znt(ss, 2)

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  # Location of the output directories. 
  outdir = '/media/My Passport/RBSP/pickles/'

  # Each event has its own directory full of pickles. 
  for pkldir in os.listdir(outdir):

    # Load all of the pickles in this directory. Offset the time array to start
    # at midnight (rather than in 1970). 
    data = load(outdir + pkldir + '/')
    data['time'] = data['time'] - data['time'][0]

    # Get the probe name and the timestamp from the directory name. From the
    # timestamps, find the data indeces corresponding to the start and end of
    # the event. 
    probe, date, t0, t1 = pdtt(pkldir)
    i0, i1 = np.argmax(data['time'] > t0), np.argmax(data['time'] > t1)

    # Slice the event out of the full day of data. Also compute a ten-minute
    # rolling average of the magnetic field to estimate the background field. 
    t     = data['time'][i0:i1]
    BGSE  = data['bgse'][:, i0:i1]
    B0GSE = getbg( data['time'], data['bgse'] )[:, i0:i1]
    EGSE  = data['egse'][:, i0:i1]
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

    Ex = dot(xhat, EGSE)
    Ey = dot(yhat, EGSE)
    Ez = dot(zhat, EGSE)

    PW = plotWindow(nrows=2, ncols=-2)

    title = notex( 'RBSP-' + probe.upper() + ' on ' + calendar(date) + ' from ' + clock(t0) + ' to ' + clock(t1) )

    tcoord = (t - t[0])/60

    print 'min tcoord = ', tcoord[0]
    print 'max tcoord = ', tcoord[-1]

    PW.setParams( x=tcoord, xlabel=notex('Time (min)'), xlims=(0, 10), xticks=range(11), xticklabels=('$0$', '', '$2$', '', '$4$', '', '$6$', '', '$8$', '', '$10$'), title=title )
#    PW.setParams( x=tcoord, xlabel=notex('Time (s)'), title=title )

    PW[0].setParams(ylabel=notex('Magnetic Field (nT)'))
    PW[0].setLine(Bx, 'b')
    PW[0].setLine(By, 'r')
    PW[0].setLine(Bz, 'g')

    PW[1].setParams(ylabel=notex('Electric Field (\\frac{mV}{m})'))
    PW[1].setLine(Ex, 'b')
    PW[1].setLine(Ey, 'r')
    PW[1].setLine(Ez, 'g')

    PW.render()

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


