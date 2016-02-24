#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This routine uses IDL to interface with NASA's data server, grab the electric
# and magnetic field for a list of events, clean up the data, and output it as
# SAV files. Those files are then loaded and re-saved as pickles for later
# analysis in Python. 

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
from scipy import io
from subprocess import Popen, PIPE

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  module('load idl')

  src = '/home/user1/mceachern/Desktop/rbsp/'
  out = '/media/My Passport/RBSP/pickles/'

  for sat in ('a', 'b')[:1]:

    events = read(src + 'events/events_' + sat + '.txt')[:1]

    for event in events:

      print 'event: ', event
      name = event.replace('/', '_').replace('-', '').replace(':', '') + '/'

      '''
      if os.path.exists('temp.pro'):
        os.remove('temp.pro')
#      print '\tcreating pro script to download the data as a sav file'
      append( pro(sat=sat, event=event) )
      print '\tDownloading the data using SPEDAS... '
      out, err = bash('idl -e @temp -IDL_PATH +~/Desktop/RBSP/packages/:<IDL_DEFAULT>')
      print '\tReading the IDL data into Python...'
      sav = io.readsav('/home/user1/mceachern/Desktop/RBSP/temp.sav')
      '''

      # Save position to out + name + 'x.pkl'
      # etc

      print out + name + 'x.pkl'

      '''
      # Tuples are safer to pickle than dictionaries. 
      print '\tCreating a pickle: ', pklpath
      with open(pklpath, 'wb') as handle:
        pickle.dump(sav.items(), handle, protocol=-1)
#      print '\tsanity check... '
#      with open(pklpath, 'rb') as handle:
#        x = dict( pickle.load(handle) )
#      print all( sav['time'] == x['time'] )
#      print all( sav['bgse'].flatten() == x['bgse'].flatten() )

  # Now let's look at an event. 
  with open(pklpath, 'rb') as handle:
    x = dict( pickle.load(handle) )

  t = x['time']
  # Ten minutes compared to a day. 
  N = (600*t.size)/(24*3600)
  t = t[:N] - t[0]

  BX = x['bgse'][0][:N]
  BY = x['bgse'][1][:N]
  BZ = x['bgse'][2][:N]

  BX = BX - np.mean(BX)
  BY = BY - np.mean(BY)
  BZ = BZ - np.mean(BZ)

  plt.plot(t, BX)
  plt.plot(t, BY)
  plt.plot(t, BZ)

  plt.show()
  '''

  return

# #############################################################################
# ########################################################### 
# #############################################################################

def pro(sat, event):
  return ( 'rbsp_efw_init\n' +
           'timespan,\'' + event + '\'\n' +
           'rbsp_load_emfisis,probe=' + sat + ',coord=\'gse\',cadence=\'hires\',level=\'l3\'\n' + 
           'get_data,1,time,bgse\n' + 
           'save,time,bgse,filename=\'~/Desktop/RBSP/temp.sav\'' )

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

def append(text, filename=None):
  if filename is not None:
    with open(filename, 'a') as fileobj:
      fileobj.write(text + '\n')
  return text

def read(filename):
  with open(filename, 'r') as fileobj:
    return [ x.strip() for x in fileobj.readlines() ]

def bash(command, save='stdoe.txt'):
  out, err = Popen(command.split(), stdout=PIPE, stderr=PIPE).communicate()
  return append(out, save), append(err, save)

def module(command, save='stdoe.txt'):
  out, err = bash('/usr/bin/modulecmd python ' + command, save=save)
  exec out
  return err

# ########################################################### For Importability
# #############################################################################
# #############################################################################

if __name__=='__main__':
  main()


