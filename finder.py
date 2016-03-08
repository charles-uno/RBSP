#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# Lei's got some filters in his data which are designed to screen out
# fundamental mode Pc4 pulsations. Can we come up with our own list of events? 
# Specifically, a list of fundamental mode Pc4 pulsations? 

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

from day import *

from plotmod import *


# #############################################################################
# ######################################################################## Main
# #############################################################################

def main():

  today = day(probe='a', date='2012-10-02')

  ev = today.getslice('04:50', '05:10')

  if not ev.isok():
    print 'SKIPPING ' + ev.name + ' DUE TO BAD DATA. '
    return

  # Create plot window to hold waveforms and spectra. 
  PW = plotWindow(nrows=3, ncols=2)

  # Put title and labels on the plot. 
  title = notex('TITLE')
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red) ; E (Blue)'), 
                notex('Coh. (Black) ; CSD (Green) ; Par. (Violet)') )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)

  # Index the poloidal, toroidal, and field-aligned modes. 
  modes = ('p', 't', 'z')

  # Plot waveforms. 
  PW[:, 0].setParams( **ev.coords('waveform', cramped=True) )
  [ PW[i, 0].setLine(ev.get('B' + m), 'r') for i, m in enumerate(modes) ]
  [ PW[i, 0].setLine(ev.get('E' + m), 'b') for i, m in enumerate(modes) ]

  # Plot coherence. 
  PW[:, 1].setParams( **ev.coords('coherence', cramped=True) )
  [ PW[i, 1].setLine( ev.csd(m), 'g' ) for i, m in enumerate(modes) ]
  [ PW[i, 1].setLine( ev.coh(m), 'k' ) for i, m in enumerate(modes) ]
  [ PW[i, 1].setLine( 0.1 + 0.8*ev.isodd(m), 'm' ) for i, m in enumerate(modes) ]

  # Show the plot. 
  return PW.render()







# Break a given date into chunks and look for fundamental mode Pc4 events in
# that chunk. Each chunk is mpc minutes long -- minutes per chunk. 
def checkdate(d, mpc=20):

  print d

  for p in ('a', 'b'):

    today = day(probe=p, date=d)

    # Scroll through each 20 minute chunk. 
    for t in range(0, 86400, 60*mpc)[:3]:

      print '\t' + p + '\t' + timestr(t)[1], 

      ev = today.getslice(t, duration=60*mpc)

      if ev.isok():
        print '\tOK'
      else:
        print '\tX'



# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


