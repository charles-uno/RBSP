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
# ######################################################################## Main
# #############################################################################


def main():

  # Keep everything nicely organized. 
  srcdir = '/home/user1/mceachern/Desktop/rbsp/'
  rundir = '/home/user1/mceachern/Desktop/rbsp/run/'
  outdir = '/media/My Passport/RBSP/pickles/'

  # Any files that get dumped should get dumped into the run directory. 
  os.chdir(rundir)

  # Loop over the two satellites. 
  for probe in ('a', 'b'):

    # Loop over the events seen by each one. 
    for event in read(srcdir + 'events/events_' + probe + '.txt'):

      # Nuke the run directory, except for the captured IDL output. 
      [ os.remove(x) for x in os.listdir(rundir) if x!='stdoe.txt' ]

      # Create a directory to hold the data. 
      name = probe + '_' + event.replace('/', '_').translate(None, '-:') + '/'

      print name

      if os.path.isdir(outdir + name):
        print '\tDATA ALREADY EXISTS'
        continue
      else:
        os.mkdir(outdir + name)

      # Create and execute an IDL script to grab position, electric field, and
      # magnetic field data for the event, and dump it into a sav file. 
      date, time = event.split('/')
      out, err = spedas( idlcode(probe=probe, date=date) )
#      print out
#      print err

      # Read in the IDL output. 
      if not os.path.exists('temp.sav'):
        print '\tNO DATA'
        continue
      else:
        temp = io.readsav('temp.sav')

      # Re-write the data in pickle format. 
      for key, arr in temp.items():
        with open(outdir + name + key + '.pkl', 'wb') as handle:
          pickle.dump(arr, handle, protocol=-1)
        print '\tcreated ' + outdir + name + key + '.pkl'

  '''

  # Let's look at a bit of data as a sanity check. 
  for pkldir in os.listdir(outdir):

    print pkldir + '/'

    # The directory gives the event's timestamp. 
    probe, date, time = pkldir.split('_')
    hh, mm, ss = int( time[0:2] ), int( time[2:4] ), int( time[4:6] )
    # Get the number of seconds from midnight to the start and end of the
    # event. Events are ten minutes long by construction. 
    t0 = ss + 60*mm + 3600*hh
    t1 = t0 + 600

    # Load the pickle files into a dictionary of arrays. 
    data = {}
    for pklname in os.listdir(outdir + pkldir):
      with open(outdir + pkldir + '/' + pklname, 'rb') as handle:
        data[ pklname[:-4] ] = pickle.load(handle)
        print '\tloaded ' + pklname

    # Find the indeces that correspond to the start and end of the event. 
    today = data['time'] - data['time'][0]
    i0 = np.argmax(today > t0)
    i1 = np.argmax(today > t1)

    print 'number of time steps for this ten-minute window: ', i1-i0


    t = today[i0:i1]
    bx = data['bgse'][0, i0:i1]
    by = data['bgse'][1, i0:i1]
    bz = data['bgse'][2, i0:i1]
    ex = data['egse'][0, i0:i1]
    ey = data['egse'][1, i0:i1]
    ez = data['egse'][2, i0:i1]
    x = data['xgse'][0, i0:i1]
    y = data['xgse'][1, i0:i1]
    z = data['xgse'][2, i0:i1]
    lshell = data['lshell'][i0:i1]
    mlt = data['mlt'][i0:i1]
    mlat = data['mlat'][i0:i1]

    # Plot a sanity check. 
    for index, label in enumerate( ('BX', 'BY', 'BZ') ):

      field = data['bgse'][index, i0:i1]
      slope, inter = np.polyfit(t, field, 1)
      linfit = slope*t + inter
      plt.plot(t, field - linfit, label=label)

    plt.legend()
    plt.show()

  '''
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

# ########################################################### For Importability
# #############################################################################
# #############################################################################

if __name__=='__main__':
  main()


