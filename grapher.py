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

from random import randrange

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
  return date[0:4] + '--' + date[4:6] + '--' + date[6:8]

# Compute clock time from a count of seconds from midnight. 
def clock(sfm, seconds=False):
  hh = sfm/3600
  if seconds:
    mm = (sfm%3600)/60
    ss = sfm%60
    return znt(hh, 2) + ':' + znt(mm, 2) + ':' + znt(ss, 2)
  else:
    mm = int( format( (sfm%3600)/60., '.0f') )
    return znt(hh, 2) + ':' + znt(mm, 2)


# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  # Location of the output directories. 
  outdir = '/media/My Passport/RBSP/pickles/'

  # Each event has its own directory full of pickles. 
  for pkldir in os.listdir(outdir)[2:3]:

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
    Bx = dot(xhat, BGSE - B0GSE)
    By = dot(yhat, BGSE - B0GSE)
    Bz = dot(zhat, BGSE - B0GSE)
    Ex = dot(xhat, EGSE)
    Ey = dot(yhat, EGSE)
    Ez = dot(zhat, EGSE)

    # Also compute the probe's average position. 
    L = np.average( data['lshell'][i0:i1] )
    MLT = np.average( data['mlt'][i0:i1] )
    mlat = np.average( data['mlat'][i0:i1] )

    # Set up a plot window. Use double-wide columns. 
    ylabels = ( notex('Poloidal (\\frac{mV}{m}) (nT)'), 
                notex('Poloidal FFT (\\frac{mV}{m}) (nT)'), 
                notex('Toroidal (\\frac{mV}{m}) (nT)'), 
                notex('Field-Aligned (\\frac{mV}{m}) (nT)') )
    PW = plotWindow(nrows=len(ylabels), ncols=-2)

    # Set the title and labels. 
    title = notex('In Situ Electric (Blue) and Magnetic (Red) Fields')
    xlabel = notex( 'Time (hh:mm) on ' + calendar(date) )
    collabel = ( notex( 'RBSP--' + probe.upper() ) + ' \\qquad L \\!=\\! ' +
                 format(L, '.1f') + ' \\qquad ' + notex( clock(3600*MLT) +
                 ' MLT') + ' \\qquad ' + format(mlat, '.1f') + '^\\circ' +
                 notex(' Magnetic Latitude') )
    PW.setParams(title=title, collabels=[collabel], xlabel=xlabel)
    [ PW[i].setParams(ylabel=ylbl) for i, ylbl in enumerate(ylabels) ]

    # Manually set the ticks and tick labels on the x axis. 
    xticks, xticklabels = range(11), ['']*11
    for i in range(1, 11, 2):
      xticklabels[i] = '$' + notex( clock(t[0] + 60*i) ) + '$'
    PW.setParams(x=( t - t[0] )/60, xticks=xticks, xticklabels=xticklabels)

    # Add the field components to the plot. 
    PW[0].setLine(Ey, 'b')
    PW[0].setLine(Bx, 'r')

    PW[2].setLine(Ex, 'b')
    PW[2].setLine(By, 'r')

    PW[3].setLine(Ez, 'b')
    PW[3].setLine(Bz, 'r')

    # Let's do a quick Fourier transform of the poloidal components and plot
    # the strongest component (in the Pc4 range) right under the poloidal
    # waveforms. That will let us estimate the phase offset. 
    dt, trange = t[1] - t[0], t[-1] - t[0]
    def harm(n):
      return np.exp(2j*np.pi*n*t/trange)

    # We only want frequencies within the Pc4 band: 45 to 150 seconds. 
    nmin, nmax = 600/150, 600/30
    Bweights = np.zeros(nmax, dtype=np.complex)
    Eweights = np.zeros(nmax, dtype=np.complex)

    for n in range(nmin, nmax):
      Bweights[n] = np.sum(Bx*harm(-n)*dt)/np.sum(harm(n)*harm(-n)*dt)
      Eweights[n] = np.sum(Ey*harm(-n)*dt)/np.sum(harm(n)*harm(-n)*dt)

    nB, nE = np.argmax( np.abs(Bweights) ), np.argmax( np.abs(Eweights) )

    tfine = np.linspace(t[0], t[-1], 1000)
    def harmfine(n):
      return np.exp(2j*np.pi*n*tfine/trange)

    PW[1].setParams(x=(tfine-tfine[0])/60)
    PW[1].setLine( np.real( harmfine(nE)*Eweights[nE] ) , 'b' )
    PW[1].setLine( np.real( harmfine(nB)*Bweights[nB] ) , 'r' )

    # Save the plot as an image. 
    PW.render()

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


