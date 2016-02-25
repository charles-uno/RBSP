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
import numpy as np
import os
from scipy import io
from subprocess import Popen, PIPE

# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  # Keep everything nicely organized. 
  srcdir = '/home/user1/mceachern/Desktop/rbsp/'
  rundir = '/home/user1/mceachern/Desktop/rbsp/run/'
  outdir = '/media/My Passport/RBSP/pickles/'

  # Any files that get dumped should get dumped into the run directory. Make
  # sure to always start with a fresh log. 
  os.chdir(rundir)
  if os.path.exists('log.txt'):
    os.remove('log.txt')

  # Loop over the two satellites. 
  for probe in ('a', 'b'):

    # Loop over the events seen by each one. 
    for event in read(srcdir + 'events/events_' + probe + '.txt'):

      # Nuke the run directory. 
      [ os.remove(x) for x in os.listdir(rundir) if x not in ('log.txt',) ]

      # The data directory is indexed by the event timestamp. 
      name = probe + '_' + event.replace('/', '_').translate(None, '-:') + '/'

      append(name, 'log.txt')

      # Check if we've already done this one. 
      if os.path.isdir(outdir + name):
        append('\tDATA ALREADY EXISTS', 'log.txt')
        continue
      else:
        os.mkdir(outdir + name)

      # Create and execute an IDL script to grab position, electric field, and
      # magnetic field data for the event, and dump it into a sav file. 
      date, time = event.split('/')
      out, err = spedas( idlcode(probe=probe, date=date) )
#      print out
#      print err

      # If IDL crashes, don't save any data. But leave the directory, so we
      # know not to bother with this event next time. 
      if 'Variable is undefined: RBSPX' in err:
        append('\tBAD DATA', 'log.txt')
        continue

      # Read in the IDL output. 
      if not os.path.exists('temp.sav'):
        append('\tNO DATA', 'log.txt')
        continue
      else:
        temp = io.readsav('temp.sav')

      # Re-write the data in pickle format. 
      print name
      for key, arr in temp.items():
        with open(outdir + name + key + '.pkl', 'wb') as handle:
          pickle.dump(arr, handle, protocol=-1)
        append('\tcreated ' + outdir + name + key + '.pkl', 'log.txt')
        print '\tcreated '+ outdir + name + key + '.pkl'

  return

# #############################################################################
# ############################################################# IDL Script Text
# #############################################################################

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

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


