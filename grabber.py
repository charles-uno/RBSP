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

  # Keep everything nicely organized. 
  src = '/home/user1/mceachern/Desktop/rbsp/'
  run = '/home/user1/mceachern/Desktop/rbsp/run/'
  out = '/media/My Passport/RBSP/pickles/'

  # Any files that get dumped should get dumped into the run directory. 
  os.chdir(run)

  # Loop over the two satellites. 
  for probe in ('a', 'b')[0:1]:

    # Loop over the events seen by each one. 
    for event in read(src + 'events/events_' + probe + '.txt')[19:20]:

      # Nuke the run directory, except for the captured IDL output. 
      [ os.remove(x) for x in os.listdir(run) if x!='stdoe.txt' ]

      # Create a directory to hold the data. 
      name = probe + '_' + event.replace('/', '_').translate(None, '-:') + '/'

      print name

      if os.path.isdir(out + name):
        print '\tDATA ALREADY EXISTS'
        continue
      else:
        os.mkdir(out + name)

      # Create and execute an IDL script to grab position, electric field, and
      # magnetic field data for the event, and dump it into a sav file. 
      date, time = event.split('/')
      out, err = spedas( idlcode(probe=probe, date=date) )
      print out, err

      # Read in the IDL output. 
      if not os.path.exists('temp.sav'):
        print '\tNO DATA'
        continue
      else:
        temp = io.readsav('temp.sav')

      # Re-write the data in pickle format. 
      for key, arr in temp.items():
        with open(out + name + key + '.pkl', 'wb') as handle:
          pickle.dump(arr, handle, protocol=-1)
        print '\tcreated ' + out + name + key + '.pkl'


  '''
  key = 'time'
  with open(run + key + '.pkl', 'rb') as handle:
    t = pickle.load(handle)

  t = t - t[0]

  print t

  i0 = np.argmax(t > t0)
  i1 = np.argmax(t > t1)

  print t[i0:i1]
  print t[i0:i1].size'''


  return


  ''' 
  # Create the IDL script to grab this event. 
  append(idlcode(sat=sat, date=date), 'temp.pro')
  # Call the IDL routine to download the event data and save it in temp.sav. 
      out, err = bash('idl -e @temp -IDL_PATH +~/Desktop/RBSP/packages/:<IDL_DEFAULT>')



      if os.path.exists('temp.pro'):
        os.remove('temp.pro')
#      print '\tcreating pro script to download the data as a sav file'
      append( pro(sat=sat, event=event) )
      print '\tDownloading the data using SPEDAS... '
      out, err = bash('idl -e @temp -IDL_PATH +~/Desktop/RBSP/packages/:<IDL_DEFAULT>')
      print '\tReading the IDL data into Python...'
      sav = io.readsav('/home/user1/mceachern/Desktop/RBSP/temp.sav')
      '''



  hh, mm, ss = [ int(x) for x in time.split(':') ]
  t0 = ss + 60*mm + 3600*hh
  # Events are 10 minutes long by construction.  
  t1 = t0 + 600

  print event
  print date, time
  print hh, mm, ss
  print t0, t1




  return




  date = '2012-10-10'

  time = '09:50:01'

  hh, mm, ss = [ int(x) for x in time.split(':') ]

  t0 = ss + 60*mm + 3600*hh
  # Events are 10 minutes long by construction.  
  t1 = t0 + 600

  day = 86400

  key = 'time'
  with open(run + key + '.pkl', 'rb') as handle:
    t = pickle.load(handle)

  t = t - t[0]

  print t

  i0 = np.argmax(t > t0)
  i1 = np.argmax(t > t1)

  print t[i0:i1]
  print t[i0:i1].size




  return


  '''
  module('load idl')
  src = '/home/user1/mceachern/Desktop/rbsp/'
  out = '/media/My Passport/RBSP/pickles/'
  for sat in ('a', 'b')[:1]:
    events = read(src + 'events/events_' + sat + '.txt')[:1]
    for event in events:
      print 'event: ', event
      name = event.replace('/', '_').replace('-', '').replace(':', '') + '/'
      if os.path.exists('temp.pro'):
        os.remove('temp.pro')
#      print '\tcreating pro script to download the data as a sav file'
      append( pro(sat=sat, event=event) )
      print '\tDownloading the data using SPEDAS... '
      out, err = bash('idl -e @temp -IDL_PATH +~/Desktop/RBSP/packages/:<IDL_DEFAULT>')
      print '\tReading the IDL data into Python...'
      sav = io.readsav('/home/user1/mceachern/Desktop/RBSP/temp.sav')
      # Save position to out + name + 'x.pkl'
      # etc
      print out + name + 'x.pkl'
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
  return ( 'rbsp_efw_init\n' +
           'timespan,\'' + date + '\'\n' +
           'rbsp_load_emfisis,probe=' + sat + ',coord=\'gse\',cadence=\'hires\',level=\'l3\'\n' + 
           'get_data,1,time,bgse\n' + 
           'save,time,bgse,filename=\'~/Desktop/rbsp/run/temp.sav\'' )
  return
  '''

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


