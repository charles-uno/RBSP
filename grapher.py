#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This script is for analysis of RBSP data, accessed through the event class
# defined in event.py. 

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

from event import *
from plotmod import *

# #############################################################################
# ######################################################################## Main
# #############################################################################

# Global variable indicating where images should be saved. 
savedir = '/home/user1/mceachern/Desktop/plots/' + now() + '/'

def main():

#  # Look at how much of RBSP's orbit is in a range where we would expect to see
#  # giant pulsations. 
#  return where()

#  return plotloc()

#  for evline in ('b\t2012-10-23\t19:15:00', 'a\t2012-10-23\t21:40:00'):
#    plot0(evline)
#  return

  # Flip through the events in random order. This makes it easier to debug -- a
  # different event pops up first every time. 
  for evline in np.random.permutation( read('goodevents.txt') ):

    if plot1(evline):
      break

#    if plot2(evline):
#      break


  return

# #############################################################################
# ########################################################## Plotting Functions
# #############################################################################

# =============================================================================
# ============================================================== Plot One Event
# =============================================================================

# Plot the waveform. Estimate frequency and phase offset with a Fourier series.
def plot1(evline):
  global savedir

#  # Force it to use a nice event. 
#  evline = 'a\t2014-05-02\t10:10:00'

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Grab and Check Event Data
  # ---------------------------------------------------------------------------

  # Based on the event record, create an event object. Make sure it's ok. 
  probe, date, time = evline.split()
  ev = event(probe=probe, date=date, time=time)
  if not ev.isok():
    print 'SKIPPING ' + ev.name + ' DUE TO BAD DATA. '
    return False

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Set Up Plot Window
  # ---------------------------------------------------------------------------

  # Create a plot window to look at both probes at the same time -- both the
  # waveforms and the spectra. 
  PW = plotWindow(nrows=3, ncols=2)

  # Set the title and labels. 
  title = notex( 'RBSP-' + probe.upper() + '  on  ' + date + '  from  ' +
                 notex( timestr(ev.t0)[1][:5] )  + '  to  ' +
                 notex( timestr(ev.t1)[1][:5] ) )

  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red) ; E (Blue)'), 
                notex('Coh. (Black) ; CSD (Green) ; Par. (Violet)') )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Add Waveforms and Spectra
  # ---------------------------------------------------------------------------

  # Index the poloidal, toroidal, and field-aligned modes. 
  modes = ('p', 't', 'z')

  # Plot waveforms. 
  PW[:, 0].setParams( **ev.coords('t', 'b', cramped=True) )
  [ PW[i, 0].setLine(ev.get('B' + m), 'r') for i, m in enumerate(modes) ]
  [ PW[i, 0].setLine(ev.get('E' + m), 'b') for i, m in enumerate(modes) ]

  # Plot coherence. 
  PW[:, 1].setParams( **ev.coords('f', 'c', cramped=True) )
  [ PW[i, 1].setLine( ev.csd(m), 'g' ) for i, m in enumerate(modes) ]
  [ PW[i, 1].setLine( ev.coh(m), 'k' ) for i, m in enumerate(modes) ]
  [ PW[i, 1].setLine( 0.1 + 0.8*ev.isodd(m), 'm' ) for i, m in enumerate(modes) ]

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------------- Screen Event
  # ---------------------------------------------------------------------------

  # Let's keep track of all the strong waves we find. 
  waves = []

  # Only look at the poloidal mode. Ignore any components of the spectra that
  # are not coherent or not spectrally dense. 
  for m in modes[:1]:
    mname = {'p':'Poloidal', 't':'Toroidal', 'z':'Parallel'}[m]

    # Find any harmonics worth talking about. 
#    harm = ev.harm(m)*ev.iscoh(m)*ev.iscsd(m)*ev.ispc4()*ev.has('Ey')
    harm = ev.harm(m)*ev.iscoh(m)*ev.ispc4()*ev.has('Bx')*ev.has('Ey')

    for i in np.nonzero(harm==1)[0]:
#    for i in np.nonzero(harm)[0]:

      waves.append( {'strength':ev.coh(m)[i]*ev.csd(m)[i], 
                     'frequency':ev.frq()[i],
                     'harmonic': harm[i], 
                     'lag':ev.lag(m)[i],
                     'mode':mname, 
                     'mlt':ev.avg('mlt'),
                     'lpp':ev.lpp,
                     'va':ev.va(m)[i],
                     'mlat':ev.avg('mlat'),
                     'lshell':ev.avg('lshell'),
                     'compression':ev.comp()[i] } )

  if len(waves) == 0:
    print 'SKIPPING ' + ev.name + ' DUE TO NO SUITABLE MODES. '
    return False

#  print 'CSD magnitudes: '
#  for c in ev.csd('p'):
#    print '\t', c

#  print 'FFT frequencies... note that there are twice as many compared to the CSD'
#  for x in ev.fft('Ep'):
#    print '\t', np.abs(x)

  print 'root unnormalized autospectral densities. '
  print 'freq\tBx ASD\tBx FFT\tEy ASD\tEy FFT'
  for i in range(10):
    print format(ev.frq()[i], '5.1f') + '\t' + format(np.sqrt( ev.asd('Bx') )[i], '5.1f') + '\t' + format(np.abs( ev.fft('Bx') )[i], '5.1f') + '\t' + format(np.sqrt( ev.asd('Ey') )[i], '5.1f') + '\t' + format(np.abs( ev.fft('Ey') )[i], '5.1f')


  print 'max poloidal E, B'

  print '\t', np.max( ev.get('Ep') ), np.max( ev.get('Bp') )

  sidelabels = []
  for w in waves:
    sidelabels.append( format(w['frequency'], '.0f') + notex('mHz') + ' \\qquad ' +
                       format(w['lag'], '+.0f') + '^\\circ' + ' \\qquad ' + 
                       ' | \\widetilde{B_z} / \\widetilde{B_x} | \\! = \\! ' +
                       format(100*w['compression'], '.0f')  + '\\% \\qquad ' + 
                       ' | \\widetilde{E_y} / \\widetilde{B_x} | \\! = \\! ' +
                       format(1e3*w['va'], '.0f')  + notex('km/s') )
  PW.setParams( sidelabel=' $\n$ '.join( [ ev.label() ] + sidelabels) )


  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Create the Plot
  # ---------------------------------------------------------------------------

  if '-i' not in argv:
    print 'Plotting ' + ev.name
    return PW.render()
  else:
    if not os.path.exists(savedir):
      os.mkdir(savedir)
    return not PW.render(savedir + ev.name + '.png')

# =============================================================================
# ============================================= Where Does RBSP Spend Its Time?
# =============================================================================

def where():

  # What range of lshell values do we want to look at? 
  lmin, lmax = 5.6, 6.5

  # Pick a random date and probe. We don't benefit from looping over them all,
  # since they are not randomly distributed in time. 
  date = choice( ulist( line.split()[1] for line in read('events.txt') ) )
  probe = choice( ('a', 'b') )

  # Grab the data. If the lshell values are bad, reroll. 
  ev = event(probe=probe, date=date)
  if not np.all( np.isfinite(ev.lshell) ):
    print probe, date, ' X'
    return where()

  # The orbit is about nine hours. If we pick the largest value between 09:00
  # and 15:00, that should be near apogee for an orbit that does not cross
  # midnight in either direction. 
  ibuff = int(ev.t.size*9./24)
  iapo = ibuff + np.argmax( ev.lshell[ibuff:-ibuff] )

  # Find the minimum in the 9 hours before the perigee. 
  iper0 = np.argmin( ev.lshell[iapo - ibuff:iapo] ) + (iapo - ibuff)
  iper1 = np.argmin( ev.lshell[iapo:iapo + ibuff] ) + iapo

  # Get the actual apogee. 
  iapo = np.argmax( ev.lshell[iper0:iper1] ) + iper0

  # What fraction of this orbit is spent between lmin and lmax? Assume time
  # steps of uniform size. 
  ttotal = ev.t[iper1] - ev.t[iper0]
  inrange = ( ev.lshell > lmin )*( ev.lshell < lmax )
  pinrange = np.sum( inrange[iper0:iper1] )*1./(iper1 - iper0)
  tinrange = ttotal*pinrange

  # Just use Matplotlib for now. We can clean it up later if necessary. 
  '''
  PW = plotWindow()
  PW.setLine(ev.t, ev.lshell, 'k')
  [ PW.setLine( np.ones(ev.t.shape)*lval, 'r:' ) for lval in (lmin, lmax) ]
  PW.render()
  '''

  # Plot the probe's path. 
  plt.plot(ev.t, ev.lshell, 'k')

  # Indicate prime lshell values for seeing giant pulsations. 
  plt.plot(ev.t, np.ones(ev.t.shape)*lmin, 'r:')
  plt.plot(ev.t, np.ones(ev.t.shape)*lmax, 'r:')

  # Indicate the orbit we're tallying. 
  plt.plot(np.ones(100)*ev.t[iper0], np.linspace(0, 7, 100), 'g:')
  plt.plot(np.ones(100)*ev.t[iper1], np.linspace(0, 7, 100), 'g:')

  # Give numbers for the orbit duration and how long is spent in Pg range. 
  plt.gca().text(s=timestr(tinrange)[1][:5] + ' (' +
                   format(100*pinrange, '.0f') + '%)' + ' with ' + str(lmin) +
                   ' < L < ' + str(lmax) , x=ev.t[iapo], y=0.5*(lmin + lmax),
                   horizontalalignment='center')
  plt.gca().text(s=timestr( ev.t[iper1] - ev.t[iper0] )[1][:5] + ' Orbit Time',
                 x=0.5*( ev.t[iper0] + ev.t[iper1] ), y=1, 
                 horizontalalignment='center')

  # Clean up the plot a little bit. 
  plt.gca().set_xlim( [0, 86400] )
  plt.gca().set_ylim( [0, 7] )

  return plt.show()

# =============================================================================
# ================================================= Plot of All Event Locations
# =============================================================================

def plotloc():

#  # Debug: force it to regenerate the list each time. 
#  if os.path.exists('locations.txt'):
#    os.remove('locations.txt')

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------ List Event Locations
  # ---------------------------------------------------------------------------

  if not os.path.exists('locations.txt'):

    # Flip through each event in the listing. 
    for evline in np.random.permutation( read('events.txt') ):

      # Create an event object. Make sure it's usable. 
      probe, date, time = evline.split()
      ev = event(probe=probe, date=date, time=time)
      if not ev.isok():
        print 'SKIPPING ' + ev.name + ' DUE TO BAD DATA. '
        continue

      # Find the strong Pc4 poloidal waves in this event. 
      waves = []
      for m in ('p', 't', 'z')[:1]:
        mname = {'p':'Poloidal', 't':'Toroidal', 'z':'Parallel'}[m]
        # Find any harmonics worth talking about. 
        harm = ev.harm(m)*ev.iscoh(m)*ev.iscsd(m)*ev.ispc4()
        for i in np.nonzero(harm)[0]:
          waves.append( {'strength':ev.coh(m)[i]*ev.csd(m)[i], 
                         'frequency':ev.frq()[i],
                         'harmonic': harm[i], 
                         'mode':mname, 
                         'mlt':ev.avg('mlt'),
                         'lpp':lpp(evline),
                         'mlat':ev.avg('mlat'),
                         'lshell':ev.avg('lshell') } )

      # If no modes were found, bail. 
      if len(waves) == 0:
        print 'SKIPPING ' + ev.name + ' DUE TO NO SUITABLE MODES. '
        continue

      # If there are multiple strong modes, keep only the strongest. It's
      # possible, in principle, to have a strong even harmonic and a strong odd
      # harmonic superimposed... but it's unlikely. 
      if len(waves) > 1:
        maxstrength = max( w['strength'] for w in waves )
        maxwave = [ w for w in waves if w['strength'] == maxstrength ][0]
      else:
        maxwave = waves[0]

      # Save the mode, its location, and the location of the plasmapause. 
      append(format(maxwave['harmonic'], '.0f') + '\t' +
             format(maxwave['lshell'], '.1f') + '\t' + 
             format(maxwave['mlt'], '.1f') + '\t' + 
             format(maxwave['mlat'], '.1f') + '\t' + 
             format(maxwave['lpp'], '.1f') + '\t' + 
             format(maxwave['frequency'], '.0f'), 'locations.txt')
      print 'Saving ' + ev.name

  # ---------------------------------------------------------------------------
  # ------------------------------------------------------ Load Event Locations
  # ---------------------------------------------------------------------------

  locs = read('locations.txt')
  locarr = np.array( [ [ float(x) for x in l.split('\t') ] for l in locs  ] )

  # Split each row into its own 1D array. 
  harm = locarr[:, 0]
  lshell = locarr[:, 1]
  mlt = locarr[:, 2]
  mlat = locarr[:, 3]
  lpp = locarr[:, 4]
  freq = locarr[:, 5]

  # ---------------------------------------------------------------------------
  # -------------------------------------------------- Plot Event Location Data
  # ---------------------------------------------------------------------------

  # Indicate which places are odd or even harmonics. 
  isodd = harm%2 == 1
  iseven = harm == 2

  # Split between large and small plasmasphere. 
  lppavg = np.mean(lpp)
  isbigpp = lpp > lppavg
  issmallpp = lpp < lppavg

  # Sanity check... scatter plot of frequency against lshell. 
  plt.scatter(lshell, freq)

  # Round lshell to the nearest half integer. 
  lround = 0.5*np.round(2*lshell)
  lbins = sorted( set(lround) )

  bigf, smallf = [], []

  # Let's also plot the average at each lshell. 
  for l in lbins:

    # Find the events that are at this lshell. 
    ishere = lround == l

    indbig = np.nonzero(isbigpp*ishere*isodd)[0]
    if indbig.size>0:
      bigf.append( np.mean( freq[indbig] ) )
    else:
      bigf.append(None)

    indsmall = np.nonzero(issmallpp*ishere*isodd)[0]
    if indsmall.size>0:
      smallf.append( np.mean( freq[indsmall] ) )
    else:
      smallf.append(None)

  plt.plot(lbins, bigf, 'r')
  plt.plot(lbins, smallf, 'b')

  plt.show()

  return

# =============================================================================
# =========================================================== Plot a Long Event
# =============================================================================

# Plot the waveform. Estimate frequency and phase offset with a Fourier series.
def plot0(evline):
  global savedir

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Grab and Check Event Data
  # ---------------------------------------------------------------------------

  # Based on the event record, create an event object. Make sure it's ok. 
  probe, date, time = evline.split()
  ev = event(probe=probe, date=date, time=time, mins=60)
  if not ev.isok():
    print 'SKIPPING ' + ev.name + ' DUE TO BAD DATA. '
    return False

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Set Up Plot Window
  # ---------------------------------------------------------------------------

  PW = plotWindow(nrows=3, ncols=-2)
  title = ev.label()
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  PW.setParams(title=title, rowlabels=rowlabels)

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Add Waveforms and Spectra
  # ---------------------------------------------------------------------------

  # Index the poloidal, toroidal, and field-aligned modes. 
  modes = ('p', 't', 'z')

  # Plot waveforms. 
  PW.setParams( **ev.coords('t', 'b') )
  [ PW[i].setLine(ev.get('B' + m), 'r') for i, m in enumerate(modes) ]
  [ PW[i].setLine(ev.get('E' + m), 'b') for i, m in enumerate(modes) ]

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Create the Plot
  # ---------------------------------------------------------------------------

  if '-i' not in argv:
    print 'Plotting ' + ev.name
    return PW.render()
  else:
    if not os.path.exists(savedir):
      os.mkdir(savedir)
    return not PW.render(savedir + ev.name + '.png')

# =============================================================================
# =========================================== Plot Both Probes During One Event
# =============================================================================

def plot2(evline):
  global savedir

#  # Force it to use a nice event. 
#  evline = 'b\t2012-10-23\t19:40:00'

  # ---------------------------------------------------------------------------
  # ------------------------------------------------- Grab and Check Event Data
  # ---------------------------------------------------------------------------

  # Based on the event record, create a pair of event objects. 
  probe, date, time = evline.split()
  ev = { 'a':event(probe='a', date=date, time=time),
         'b':event(probe='b', date=date, time=time) }

  # If the events have flawed data, start over. 
  if not ev['a'].isok() or not ev['b'].isok():
    print 'SKIPPING ' + ev[probe].name + ' DUE TO BAD DATA. '
    return False

  # ---------------------------------------------------------------------------
  # -------------------------------------------------------- Set Up Plot Window
  # ---------------------------------------------------------------------------

  # Create a plot window to look at both probes at the same time -- both the
  # waveforms and the spectra. 
  PW = plotWindow(nrows=4, ncols=3)

  # Set the title and labels. Label the rows with probe position data. 
  title = notex(probe.upper() + '  ' + date + '  ' + time)
  collabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  rowlabels = ( ev['a'].lab(), ev['b'].lab(), ev['a'].lab(), ev['b'].lab() )
  PW.setParams(title=title, collabels=collabels, rowlabels=rowlabels)

  # ---------------------------------------------------------------------------
  # ----------------------------------------- Add Waveforms and Spectra to Plot
  # ---------------------------------------------------------------------------

  # Index all probe/mode combinations. 
  pm = [ (p, m) for p in ('a', 'b') for m in ('p', 't', 'z') ]

  # For each probe and mode, compute coherence and (complex phase) lag. 
  frq = [ ev[p].coh(m)[0] for p, m in pm ]
  coh = [ ev[p].coh(m)[1] for p, m in pm ]
  lag = [ ev[p].lag(m)[1] for p, m in pm ]

  # Iterate over the probe/mode combinations. 
  for i, (p, m) in enumerate(pm):

    # Plot the fields. 
    PW[i/3, i%3].setParams( **ev[p].coords('t', 'b', cramped=True) )
    PW[i/3, i%3].setLine(ev[p].get('E' + m), 'b')
    PW[i/3, i%3].setLine(ev[p].get('B' + m), 'r')

    # Plot the spectra. Scale to the unit interval. 
    PW[i/3 + 2, i%3].setParams( **ev[p].coords('f', 'c', cramped=True) )
    PW[i/3 + 2, i%3].setLine( frq[i], coh[i], 'k')
#    PW[i/3 + 2, i%3].setLine( frq[i], lag[i]/360 + 0.5, 'g')

#    # Summarize the spectra. 
#    fr, ph, pc = ev[p].fpp(m)
#    PW[i/3, i%3].setParams(text=fr + ' \\qquad ' + ph)

  # ---------------------------------------------------------------------------
  # -------------------------------------- Look for Fundamental Mode Pc4 Events
  # ---------------------------------------------------------------------------

  fundlabels = []
  for i, f in enumerate( frq[0] ):

    # Fundamental FLRs should be around 10mHz. 
    if not 6 < f < 20:
      continue

    # Find places where the coherence is significantly above average. 
    iscoherent = np.array( [ c[i] > np.mean(c) + 2*np.std(c) for c in coh ] )

    # In fundamental modes, mlat and phase lag have opposite signs. 
    mlatsign = np.array( [ np.sign( ev[p].avg('mlat') ) for p, m in pm ] )
    lagsign = np.array( [ np.sign( l[i] ) for l in lag ] )

    # Skip anything within one degree of the magnetic equator. 
    isoffeq = np.array( [ np.abs( ev[p].avg('mlat') ) > 1 for p, m in pm ] )

    # Find the frequency/mode combinations -- if any -- which are coherent and
    # fundamental. 
    iscoherentfund = (mlatsign!=lagsign)*iscoherent*isoffeq

#    # Only keep events which trigger in the poloidal mode. 
#    if not any( iscoherentfund[::3] ):
#      continue

    # Can we find any simultaneous observations by both probes? 
    if not any( iscoherentfund[:3] ) or not any( iscoherentfund[3:] ):
      continue

    # Note the nontrivial fundamental modes, if any. 
    for x in np.nonzero(iscoherentfund)[0]:
      modenames = {'p':'Poloidal', 't':'Toroidal', 'z':'Parallel'}
      fundlabels.append( pm[x][0].upper() + '  ' + format(f, '.0f') + 'mHz  ' +
                         modenames[ pm[x][1] ] )

  # If no fundamental modes were found, bail. 
  if not fundlabels:
    print 'SKIPPING ' + ev[probe].name + ' DUE TO NO SUITABLE MODES. '
    return False

  # Otherwise, note the events on the plot. 
  PW.setParams( sidelabel=' $\n$ '.join( notex(x) for x in fundlabels ) )

  # ---------------------------------------------------------------------------
  # ----------------------------------------------------------- Create the Plot
  # ---------------------------------------------------------------------------

  if '-i' not in argv:
    print 'Plotting ' + ev[probe].name
    return PW.render()
  else:
    if not os.path.exists(savedir):
      os.mkdir(savedir)
    return not PW.render(savedir + ev[probe].name + '.png')


# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# Remove all duplicates from a list or generator. Does not preserve order. 
def ulist(x):
  return list( set( list(x) ) )

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


