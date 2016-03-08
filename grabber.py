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

from event import timeint, timestr
try:
  import cPickle as pickle
except ImportError:
  import pickle
import numpy as np
import os
from plotmod import now
from scipy import io
from subprocess import Popen, PIPE
from sys import stdout

# #############################################################################
# ######################################################################## Main
# #############################################################################

# Keep everything nicely organized. 
srcdir = '/home/user1/mceachern/Desktop/rbsp/'
rundir = '/home/user1/mceachern/Desktop/rbsp/run_' + now() + '/'
outdir = '/media/My Passport/rbsp/pkls/'

def main():

  # Any files that get dumped should get dumped into the run directory. 
  if not os.path.exists(rundir):
    os.makedirs(rundir)
  os.chdir(rundir)

  # Load one chunk of dates at a time. We may need to clean out the scratch
  # directory in between chunks. 
  for d in listdates(onlydai=False)[100:200]:
    grabdate(d)

  return

# #############################################################################
# ############################################################# Grab Event Data
# #############################################################################

# =============================================================================
# ========================================================== List Dates to Grab
# =============================================================================

# By default, list the days from 01 October 2012 to 31 July 2014 -- a range
# over which RBSP's orbit precesses all the way around the planet. Optionally,
# instead, limit the dates to those on which Lei recorded a poloidal Pc4 event. 
def listdates(onlydai=False):
  # Grab Lei's event list and remove duplicate dates. 
  if onlydai:
    events = read(srcdir + 'events.txt')
    dates = sorted( set( line.split()[1] for line in events ) )
  # Grab all of the dates from October 2012 to July 2014. 
  else:
    dates = []
    d = '2012-10-01'
    while d < '2014-08-01':
      dates.append(d)
      d = timestr(timeint(date=d) + 86400)[0]
  return dates

# =============================================================================
# ======================================================= Grab Date for One Day
# =============================================================================

# Given a single date, use SPEDAS to load the RBSP data, dump it, read it, and
# save it as pickles. 
def grabdate(d):
  global rundir, srcdir, outdir
  # Limit output to one line per date. 
  status(d)
  # Loop over both probes. 
  for p in ('a', 'b'):
    status(p)
    # Nuke the run directory. Leave stdout and stderr. 
    [ os.remove(x) for x in os.listdir(rundir) if x not in ('stdoe.txt',) ]
    # Create and execute an IDL script to grab position, electric field, and
    # magnetic field data for the day and and dump it into a sav file. 
    out, err = spedas( idlcode(probe=p, date=d) )
    # Make sure there's somewhere to store the pickles. 
    pkldir = outdir + d.replace('-', '') + '/' + p + '/'
    if not os.path.exists(pkldir):
      os.makedirs(pkldir)
    # Read in the IDL output. 
    if not os.path.exists('temp.sav'):
      status('X')
      continue
    else:
      temp = io.readsav('temp.sav')
    # Rewrite the data as pickles. (Pickles are Python data files. They are
    # reasonably efficient in terms of both storage size and load time.)
    for key, arr in temp.items():
      with open(pkldir + key + '.pkl', 'wb') as handle:
        pickle.dump(arr, handle, protocol=-1)
    # Acknowledge successful date access. 
    status('OK')
  # Move to the next line. 
  return status()

# =============================================================================
# ========================================================= IDL Script Assembly
# =============================================================================

# This routine reads in a bunch of IDL commands from crib.pro then modifies
# them slightly and returns them. 
def idlcode(probe, date):
  # Read in the crib sheet. 
  crib = read('../crib.pro')
  # Find the lines that define the date and the probe. 
  idate = np.argmax( [ line.startswith('date = ') for line in crib ] )
  iprobe = np.argmax( [ line.startswith('probe = ') for line in crib ] )
  # Change those lines to describe this event. 
  crib[idate] = 'date = \'' + date + '\''
  crib[iprobe] = 'probe = \'' + probe + '\''
  # Return the list of IDL commands as a newline-delimited string. 
  return '\n'.join(crib)






'''
# #############################################################################
# ############################################################# Grab Event Data
# #############################################################################

# Let's grab all of the data from October 2012 to July 2014. That includes a
# full precession around the planet (hence Lei's analysis). Running will take
# a while, but the size of the data shouldn't really be a problem. 
def grabdate(dates):

  # Keep everything nicely organized. 
  srcdir = '/home/user1/mceachern/Desktop/rbsp/'
  rundir = '/home/user1/mceachern/Desktop/rbsp/run/'
  outdir = '/media/My Passport/rbsp/pkls/'

  # Any files that get dumped should get dumped into the run directory. 
  os.chdir(rundir)

  # For each of these dates -- working a chunk of days at a time in case we
  # need to clean out the scratch directory -- grab RBSP data for both probes. 
  for d in dates[:100]:
    # Limit output to one line per date. 
    status(d)

    # For each day, make sure we have data for both probes (even if an event
    # was only identified for one of them). 
    for p in ('a', 'b'):
      # Mark a status for each probe. 
      status(p)

      # Make a directory for this day of data. 
      pkldir = outdir + d.replace('-', '') + '/' + p + '/'
      if not os.path.exists(pkldir):
        os.makedirs(pkldir)

      # Nuke the run directory. Leave stdout and stderr. 
      [ os.remove(x) for x in os.listdir(rundir) if x not in ('stdoe.txt',) ]

      # Create and execute an IDL script to grab position, electric field, and
      # magnetic field data for the day and and dump it into a sav file. 
      out, err = spedas( idlcode(probe=p, date=d) )

      # Read the IDL output. 
      if not os.path.exists('temp.sav'):
        status('X')
        continue
      else:
        temp = io.readsav('temp.sav')

      # Rewrite the data as pickles. (Pickles are Python data files. They are
      # reasonably efficient in terms of both storage size and load time.)
      for key, arr in temp.items():
        with open(pkldir + key + '.pkl', 'wb') as handle:
          pickle.dump(arr, handle, protocol=-1)

      # Acknowledge successful date access. 
      status('OK')

    # Move to the next line. 
    status()

  return

# #############################################################################
# ############################################################# IDL Script Text
# #############################################################################
'''

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# Append text to a file. 
def append(text, filename=None):
  if filename is not None:
    with open(filename, 'a') as fileobj:
      fileobj.write(text + '\n')
  return text

# Make a call as if from the command line. 
def bash(command, save='stdoe.txt'):
  out, err = Popen(command.split(), stdout=PIPE, stderr=PIPE).communicate()
  return append(out, save), append(err, save)

# Load, unload, or list modules. 
def module(command, save='stdoe.txt'):
  out, err = bash('/usr/bin/modulecmd python ' + command, save=save)
  exec out
  return err

# Read in a file as a list of lines. 
def read(filename):
  with open(filename, 'r') as fileobj:
    return [ x.strip() for x in fileobj.readlines() ]

# Dump a bunch of commands into a (temporary) IDL batch file, load IDL, and
# execute that batch file. 
def spedas(command):
  if os.path.exists('temp.pro'):
    os.remove('temp.pro')
  module('load idl')
  os.environ['ROOT_DATA_DIR'] = '/export/scratch/users/mceachern/RBSP/'
  append('PREF_SET, \'IDL_DLM_PATH\', \'<IDL_DEFAULT>\' + ' + 
         'PATH_SEP(/SEARCH_PATH) + \'~/Desktop/rbsp/incl\', /COMMIT',
         'temp.pro')
  append(command, 'temp.pro')
  return bash('idl -e @temp -IDL_PATH +~/Desktop/rbsp/packages/:<IDL_DEFAULT>')

# Print to the terminal without advancing the line. 
def status(text=None):
  if text is None:
    print ''
    return
  else:
    print text + '\t',
    return stdout.flush()

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


