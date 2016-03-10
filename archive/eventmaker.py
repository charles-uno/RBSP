#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# Filter Lei's raw events into a single, nicely organized file. 


import numpy as np

from plotmod import znt


# Note that Lei's events are rounded to the nearest multiple of 10 minutes. At most, this is a change of 28 seconds; in most cases, it's less than 5 seconds. 

def main():

  for probe in ('a', 'b'):

    events = read('events/events_' + probe + '.txt')

    for ev in events:
      
      date, time = ev.split('/')

      h0, m0, s0 = [ int(x) for x in time.split(':') ]

      sfm = 3600*h0 + 60*m0 + s0

      mfm = int( format(sfm/60., '.0f') )

      hh, mm, ss = mfm/60, mfm%60, 0

      append(probe + '\t' + date + '\t' + znt(hh, 2) + ':' + znt(mm, 2) + ':00', 'events.txt')

  return

# Append text to a file. 
def append(text, filename=None):
  if filename is not None:
    with open(filename, 'a') as fileobj:
      fileobj.write(text + '\n')
  return text


# Read in a file as a list of lines. 
def read(filename):
  with open(filename, 'r') as fileobj:
    return [ x.strip() for x in fileobj.readlines() ]


# Compute clock time from a count of seconds from midnight. 
def clock(sfm, seconds=False):
  if not np.isfinite(sfm):
    return '??:??'
  hh = sfm/3600
  if seconds:
    mm = (sfm%3600)/60
    ss = sfm%60
    return znt(hh, 2) + ':' + znt(mm, 2) + ':' + znt(ss, 2)
  else:
    mm = int( format( (sfm%3600)/60., '.0f') )
    return znt(hh, 2) + ':' + znt(mm, 2)


# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


