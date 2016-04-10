#!/usr/bin/env python

# Charles McEachern

# Spring 2016

# #############################################################################
# #################################################################### Synopsis
# #############################################################################

# This script makes nice-looking plots of RBSP data. 

# #############################################################################
# ######################################################### Load Python Modules
# #############################################################################

from day import *
from plotmod import *

from socket import gethostname


# #############################################################################
# ######################################################################## Main
# #############################################################################

# A timestamped directory, in case we want to save any plots. 
if gethostname()=='charles-lenovo':
  plotdir = '/home/charles/Desktop/plots/' + now() + '/'
else:
  plotdir = '/home/user1/mceachern/Desktop/plots/' + now() + '/'

# What should the minimum amplitude be? Anything below there is noise. 
thresh = 0.001 if 'nothresh' in argv else 0.05 if 'THRESH' in argv else 0.01

# Bin all position plots the same way. 
pargs = {'dl':1, 'dm':1} if 'sharp' in argv else {'dl':2, 'dm':2, 'lmin':3, 'lmax':7}

# Make some titles based on arguments from the terminal. This is super kludgey -- sorry! 
ldiv = notex('Split at ') + ( 'L\\!=\\!L_{PP}' if 'lpp' in argv else 'L\\!=\\!5')
tname = tex('ImS') + ' \\geq ' + str(thresh) + tex('mW/m^2')

titlehelper = tname + notex(', ') + ldiv

label = ('_sharp' if 'sharp' in argv else '') + ('_nothresh' if 'nothresh' in argv else '_THRESH' if 'THRESH' in argv else '') + ('_lpp' if 'lpp' in argv else '') + ('_phase' if 'phase' in argv else '')

def main():

#  # Find an event from the list and plot it. 
#  return showevent(amp_ge=1)

#  return dungey()

  # Just tell me how many there are of each mode. 
  count()


  # Here are the plots we actually use. 

#  posplot(storm=None, save='-i' in argv)
#  allplot(storm=None, save='-i' in argv)
#  modeplot(storm=None, save='-i' in argv)
#  paramplot(name='amp', save='-i' in argv)
#  paramplot(name='f', save='-i' in argv)
#  paramplot(name='phase', save='-i' in argv)
#  modesbyparam(name='amp', save='-i' in argv)
#  modesbyparam(name='f', save='-i' in argv)
  modesbyparam(name='phase', save='-i' in argv)
#  modesbyparam(name='fwhm', save='-i' in argv)

#  # Location of the usable data. 
#  [ posplot(storm=s, save='-i' in argv) for s in (True, False, None) ]

#  # Location of all events, regardless of mode or harmonic. 
#  [ allplot(storm=s, save='-i' in argv) for s in (True, False, None) ]

#  # Location of events as a function of storm index. 
#  [ dstplot(mode, save='-i' in argv) for mode in ('p', 't') ]

#  # Location of all events, by parity and polarization. 
#  [ modeplot(storm=s, save='-i' in argv) for s in (True, False, None) ]

#  # Location of simultaneous poloidal-toroidal events. 
#  doubleplot(save='-i' in argv)

#  # Location of poloidal events by compressional coupling. 
#  [ azmplot(storm=s, save='-i' in argv) for s in (True, False, None) ]

#  # Location of poloidal events by spectral width. 
#  [ [ fwhmplot(mode, split=1.3, storm=s, save='-i' in argv) for mode in ('p', 't') ] for s in (True, False, None) ]

#  # Location of events inside or outside the plasmapause. 
#  [ llppplot(mode, save='-i' in argv) for mode in ('p', 't') ]
#  ppplot(save='-i' in argv)


  return

# How many events of each mode? 
def count():
  global pargs
  pos = getpos(**pargs)
  # Count modes individually. 
  for mode in ('P1', 'P2', 'T1', 'T2'):
    eh = eventhist(pos['hargs'], mode=mode)
    print mode, np.sum(eh)
  # Count double events. Actually, these get counted in the loader. 
  for pmode in ( ('P1', 'P2') ):
    for tmode in ( ('T1', 'T2') ):
      dh = doublehist(pos['hargs'], pmode=pmode, tmode=tmode)

  return

# #############################################################################
# ######################################################### Show a Single Event
# #############################################################################

# =============================================================================
# ================================================= Loop Over (Filtered) Events
# =============================================================================

# This is like the routine in finder.py, but not identical. It's meant to be
# prettier. 
def showevent(date=None, time=None, probe=None, **kargs):
  # Grab the list of events, filtered by event properties.  
  events = loadevents(**kargs)
  # Filter by date, time, probe. 
  if date is not None:
    dates = g2a( e[3:13] for e in events )
    inew = np.nonzero( g2a( d == date for d in dates ) )
    events = events[inew]
  if time is not None:
    times = g2a( e[17:25] for e in events )
    inew = np.nonzero( g2a( t == time for t in times ) )
    events = events[inew]
  if probe is not None:
    probes = g2a( e[0] for e in events )
    inew = np.nonzero( g2a( p.upper() == probe.upper() for p in probes ) )
    events = events[inew]
  # If we're looking, grab one event at random from those that have passed the
  # filters. If we're making images, increment over them all. 
  events = np.random.permutation(events)
  events = events if '-i' in argv else events[:1]
  for e in events:
    prb, dat, tim = e[0], e[3:13], e[17:25]
    plotevent( day(probe=prb, date=dat).getslice(tim, duration=1800) )
  return 

# =============================================================================
# ============================================================== Plot One Event
# =============================================================================


def evtitle(d):
  if 45 < np.abs( d['phase'] ) < 135:
    harm = 'Even' if d['harm']%2 else 'Odd'
  else:
    harm = 'Traveling'
  mode = 'Poloidal' if d['mode'].upper()=='P' else 'Toroidal'
  return notex(harm + ' ' + mode + ' Wave')

def evtop(d):
  if d is None:
    return ''
  freq = fmt(d['f'], digs=1) + tex('mHz')
  ampl = fmt(d['s'], digs=3) + tex('mW/m^2')
  phas = fmt(d['phase'], digs=0) + '^{\\circ}'
  return '\\;\\;\\;\\;'.join( (ampl, freq, phas) )

# Wrapper around np.log10 to avoid complaints about zero. 
def log10(x):
  return np.log10( np.maximum(x, 1e-20) )

# Basically like the function of the same name in finder.py, but prettier. 
def plotevent(ev):
  global plotdir
  # Create plot window to hold waveforms and spectra. 
  PW = plotWindow(nrows=3, ncols=2, footlabel=True)
  # Index the poloidal, toroidal, and field-aligned modes. 
  modes = ('p', 't', 'z')
  # Plot waveforms as a function of time. 
  PW[:, 0].setParams( **ev.coords('waveform', cramped=True) )
  [ PW[i, 0].setLine(ev.get('B' + m), 'r') for i, m in enumerate(modes) ]
  [ PW[i, 0].setLine(ev.get('E' + m), 'b') for i, m in enumerate(modes) ]
  # Grab the Fourier-domain Poynting flux. It's scaled by L^3 to give values at
  # the ionosphere. 
  scomplex = [ ev.sfft(m) for m in modes ]
  sreal = [ np.abs( np.real(s) ) for s in scomplex ]
  simag = [ np.abs( np.imag(s) ) for s in scomplex ]
  stotal = [ np.abs(s) for s in scomplex ]
  # Find the waves in this event. 
  waves = [ ev.wave(m, pc4=True) for m in modes ]
  [ PW[i, 1].setParams( toptext=evtop(w) ) for i, w in enumerate(waves) ]
  # Plot the spectra. 
  PW[:, 1].setParams( **ev.coords('spectra', cramped=True) )
  [ PW[i, 1].setLine(log10(s), 'm') for i, s in enumerate(simag) ]
  [ PW[i, 1].setLine(log10(s), 'g') for i, s in enumerate(sreal) ]
  # Plot the Gaussian fit of the total Poynting flux. 
  f = np.linspace(0, 50, 1000)
  for i, w in enumerate(waves):
    args = (0, 0, 1) if w is None else ( w['s'], w['f'], w['df'] )
    PW[i, 1].setLine(f, log10( ev.gauss(f, *args) ), 'k--')
  # Plot title and labels. 
  rowlabels = ( notex('Poloidal'), notex('Toroidal'), notex('Parallel') )
  collabels = ( notex('B (Red)   E (Blue)'), 
                tex('imag') + tex('S') + notex(' (Magenta)   ') + 
                tex('real') + tex('S') + notex(' (Green)') )
  PW.setParams(collabels=collabels, footer=ev.label(), rowlabels=rowlabels)
  # Information about the wave(s) goes in the title. 
  tlist = [ evtitle(w) for w in waves if w is not None ]
  PW.setParams( title=notex('Waveforms and Spectra: ') + notex(' and ').join(tlist) )
  # Show the plot, or save it as an image. 
  if '-i' in argv:
    return PW.render(plotdir + ev.name + '.pdf')
  else:
    return PW.render()










# #############################################################################
# ########################################################## Position Histogram
# #############################################################################

def posplot(storm=None, save=False):
  global pargs, plotdir

  if storm is None:
    pos = getpos(**pargs)
  elif storm is True:
    pos = getpos(dst_lt=-30, **pargs)
  elif storm is False:
    pos = getpos(dst_ge=-30, **pargs)
  else:
    print 'Plotting storm and calm together looks real bad! '
    return

  x, y, z = [ pos[key] for key in ('x', 'y', 'z') ]
  date0, date1 = pos['dates']
  dt = np.sum(z)/48.
  # Create the plot window using the bullseye params helper function. 
  PW = plotWindow( **bep(rate=False) )

  status = {True:'Storm ', False:'Quiet ', None:''}[storm]
  title = notex('Distribution of Usable ' + status + 'Data: ' + date0 + ' to ' + date1)

  PW.setParams( title=title, unitlabel=notex('days'), zmax=26, lcorner=notex('Total: ') + format(dt, '.0f') + notex(' days') )
  # Add the data to the plot. 
  PW.setMesh( x, y, zmask(z/48.) )

  print 'Days in each L bin:'
  for zrow in z/48.:
    print np.sum(zrow)

  # Show the plot, or save it as an image. 
  if save is True:
    return PW.render(plotdir + 'pos_' +  {True:'storm', False:'calm', None:'all'}[storm] + label + '.pdf')
  else:
    return PW.render()

# #############################################################################
# ################################################ Parameter Distribution Plots
# #############################################################################

# =============================================================================
# ======================================================== Get Parameter Values
# =============================================================================

def getstats(arr):
  return { 'mean':np.mean(arr), 'median':np.median(arr), 'min':np.min(arr),
           'max':np.max(arr), 'std':np.std(arr) }

def getparam(name, **kargs):
  global thresh
  # Grab the events that line up with our filters. 
  events = loadevents(**kargs)
  # Match name to column. 
  i = {'probe':0, 'date':1, 'time':2, 'lshell':3, 'mlt':4, 'mlat':5, 'lpp':6, 
       'mode':7, 'f':8, 'fwhm':9, 'phase':10, 'amp':11, 'comp':12,
       'dst':13}[name]
  # Grab the list of values. 
  arr = g2a( float( l.split()[i] ) for l in events )
  # Figure out the appropriate range for the histogram. 
  rang = {'probe':None, 'date':None, 'time':None, 'lshell':(1, 7),
          'mlt':(0, 24), 'mlat':(-20, 20), 'lpp':(3, 7), 'mode':None, 
          'f':(7, 25), 'fwhm':(0, 5), 'phase':(0, 180), 'amp':(-2, 1), 
          'comp':(0, 1), 'dst':(-150, 150)}[name]
  # Figure out an appropriate number of bins. 
  bins = {'probe':None, 'date':None, 'time':None, 'lshell':6, 'mlt':12,
          'mlat':10, 'lpp':4, 'mode':None, 'f':18, 'fwhm':20, 'phase':18, 
          'amp':12, 'comp':10, 'dst':30}[name]
  # We want the absolute value of the phase, and we want to plot amplitude on
  # a log scale. Make sure we do this in the right order. 
  if name=='phase':
    arr = np.abs(arr)
  stats = getstats(arr)
  if name=='amp':
    arr = np.log10(arr)
  # Compute the histogram. 
  vals, edges = np.histogram(arr, range=rang, bins=bins)
  # Put points at bin centers. 
  return 0.5*( edges[1:] + edges[:-1] ), vals, stats

# =============================================================================
# ============================================= Look at Parameter Distributions
# =============================================================================

# Define Gaussian function. 
def gauss(x, *args):
  if len(args)==3:
    amp, avg, std = args
  else:
    return 0*x
  return amp*np.exp( -(x - avg)**2/(2.*std**2) )

# Fit a Gaussian to data. The guess is amplitude, mean, spread. 
def gaussfit(x, y, guess=None):
  # Try to fit a peak around the highest value. 
  if guess is None:
    guess = ( np.max(y), x[ np.argmax(y) ], 0.1*( x[-1] - x[0] ) )
  # Suppress warnings, but don't return bad data. 
  try:
    with catch_warnings():
      simplefilter('ignore')
      return curve_fit(gauss, x, y, p0=guess)[0]
  except:
    return None

# Coordinate keyword dictionary to be passed to the plot window to make a
# histogram by parameter value look nice. 
def pcoords(name):
  if name not in ('phase', 'f', 'amp', 'dst', 'fwhm'):
    return {}
  xlabel = { 'f':notex('Frequency (mHz)'), 'fwhm':notex('FWHM (mHz)'), 'phase':notex('|Phase|'), 
             'amp':tex('S') + notex(' (\\frac{mW}{m^2})'), 
             'dst':notex('DST (nT)') }[name]
  xlims = { 'f':(7, 25), 'fwhm':(0, 5), 'phase':(0, 180), 'amp':(-2, 1), 
            'dst':(-150, 150) }[name]
  xticks = { 'f':np.mgrid[7:25:7j], 'fwhm':np.mgrid[0:5:11j], 'phase':np.mgrid[0:180:5j], 
             'amp':np.mgrid[-2:1:7j], 'dst':np.mgrid[-150:150:7j] }[name]
  xticklabels = g2a( '$' + str( int(t) ) + '$' for t in xticks )
  if name=='amp':
    xticklabels = g2a( '$' + tdp(10**t) + '$' for t in xticks )
  if name=='phase':
    xticklabels = g2a( '$' + str( int(t) ) + '^\\circ$' for t in xticks )

  xticklabels[1::2] = ''
  ylabel = notex('Events')
  ylims = (0, 160)
  yticks = np.mgrid[0:160:5j]
  yticklabels = g2a( '$' + str( int(t) ) + '$' for t in yticks )
  yticklabels[1::2] = ''
  return {'xlabel':xlabel, 'xlims':xlims, 'xticks':xticks, 
          'xticklabels':xticklabels, 'ylabel':ylabel, 'ylims':ylims, 
          'yticks':yticks, 'yticklabels':yticklabels, 'ylabelpad':-2}

# Let's take a look at a plot of how FWHM depends on mode. 
def paramplot(name, save=False):
  global plotdir
  # Set up the window. 
  PW = plotWindow(ncols=2, nrows=2, colorbar=None)
  # Figure out what we're looking at, title-wise. 
  ttl = {'probe':notex('Probe'), 'date':notex('Date'), 'time':notex('Time'),
         'lshell':'L' + notex('-Shell'), 'mlt':notex('MLT'),
         'mlat':notex('Magnetic Latitude'), 'lpp':'L_{PP}', 
         'mode':notex('Mode'), 'f':notex('Frequency'), 'fwhm':notex('FWHM'),
         'phase':notex('Phase'), 'amp':notex('Amplitude'), 
         'comp':notex('Compressional Coupling'), 'dst':notex('DST')}[name]
  # For labeling, what units are we looking at?
  unit = {'probe':'', 'date':'', 'time':'', 'lshell':'', 'mlt':notex('hours'),
         'mlat':'^\\circ', 'lpp':'', 'mode':'', 'f':notex('mHz'), 
         'fwhm':notex('mHz'), 'phase':'^\\circ', 
         'amp':notex('\\frac{mW}{m^2}'), 'comp':'', 'dst':notex('nT')}[name]
  # Title and labels. 
  title = ttl + notex(' Distribution of Pc4 Events by Mode')
  rlabs = ( notex('Odd'), notex('Even') )
  clabs = ( notex('Poloidal'), notex('Toroidal') )
  PW.setParams(collabels=clabs, rowlabels=rlabs, title=title, **pcoords(name) )
  # For each mode, grab the parameter histogram and plot it. 
  modes = ('P1', 'T1', 'P2', 'T2')
  xy = [ getparam(name, mode=mode, phase=60) for mode in modes ]
  # The bin width is needed because bins are listed at their center but plotted
  # from the left edge. 
  dx = xy[0][0][1] - xy[0][0][0]
  [ PW[i].setBars(x - dx/2, y, width=dx) for i, (x, y, s) in enumerate(xy) ]
  # For the phase, do a Gaussian fit. 
  if name in ('phase', 'fwhm'):
    gfit = [ gaussfit(x, y) for x, y, stats in xy ]
    xg = np.linspace(np.min(x), np.max(x), 1000)
    yg = [ gauss(xg, *gf) if gf is not None else None for gf in gfit ]
    [ PW[i].setLine(xg, y, 'r') for i, y in enumerate(yg) if y is not None ]
    # Add a label indicating the mean and spread. 
    ymax = max( stats['median'] for x, y, stats in xy )

    print 'YMAX = ', ymax

    dp = 2 if ymax < 1 else 1 if ymax < 10 else 0
    def label(gf):

      print 'LABEL WITH DECIMAL PLACES = ', dp

      return ( format(gf[1], '.' + str(dp) + 'f') + unit + '\\pm' +
               format(gf[2], '.' + str(dp) + 'f') + unit )
    [ PW[i].setParams( toptext=label(gf) ) for i, gf in enumerate(gfit) if gf is not None ]
  # Otherwise, just indicate mean and median. 
  else:
    def label(stats, meandigs=2, meddigs=3):
      return ( notex('Mean: ') + fmt( stats['mean'], digs=meandigs ) + unit +
               notex('  Med: ') + fmt( stats['median'], digs=meddigs ) + unit )
    meandigs, meddigs = (2, 3) if name=='amp' else (0, 0)
    [ PW[i].setParams( toptext=label(s, meandigs, meddigs) ) for i, (x, y, s) in enumerate(xy) ]
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + name + '.pdf')
  else:
    return PW.render()

# =============================================================================
# ============================================ Events Sliced by Parameter Value
# =============================================================================

def modesbyparam(name, save=False):
  global pargs, plotdir
  # Grab the position histogram for normalization. Ignore DST. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Create a plot window to show different subsets of the events. 
  PW = plotWindow( ncols=3, nrows=4, **bep() )
  # Set title and labels. 
  rlabs = ( notex('Odd\nPoloidal'), 
            notex('Even\nPoloidal'), 
            notex('Odd\nToroidal'), 
            notex('Even\nToroidal') )
  funit, aunit, to = notex('mHz'), notex('\\frac{mW}{m^2}'), notex('---')
  # How are we splitting the columns? 
  amps, fs, fwhms, phases = (0.01, 0.03, 0.1), (7, 11, 17, 25), (1.7, 1.4, 1.1), (60, 75, 85, None)
  # Column labels. 

  deg = '^\\circ'

  clabs = {'amp':[ '\\geq ' + str(a) + aunit for a in amps ], 
           'f':[ str(f0) + funit + to + str(f1) + funit for f0, f1 in zip( fs[:-1], fs[1:] ) ],
           'fwhm':[ '<' + str(f) + funit for f in fwhms ],
#           'phase':[ str(p) + '^\\circ' + to + str(180 - p) + '^\\circ' for p in phases ]}[name]
           'phase':[ str(phases[0]) + deg + to + str(phases[1]) + deg + notex(', ') + str(180 - phases[1]) + deg + to + str(180 - phases[0]) + deg,
    str(phases[1]) + deg + to + str(phases[2]) + deg + notex(', ') + str(180 - phases[2]) + deg + to + str(180 - phases[1]) + deg,
    str(phases[2]) + deg + to + str(180 - phases[2]) + deg ]}[name]

  title = notex( 'Distribution of Pc4 Events by Mode and ' + {'amp':'Amplitude', 'f':'Frequency', 'fwhm':'Spectral Width', 'phase':'Phase'}[name] )

  PW.setParams(collabels=clabs, rowlabels=rlabs, title=title)
  # Setting up filters to grab the events for each column. 
  filters = { 'phase':[ {'phase': p0, 'antiphase':p1} for p0, p1 in zip( phases[:3], phases[1:] ) ], 
              'f':[ {'f_ge': f0, 'f_lt':f1} for f0, f1 in zip( fs[:-1], fs[1:] ) ],
              'fwhm':[ {'fwhm_lt': f} for f in fwhms ], 
              'amp':[ {'amp_ge': a} for a in amps ] }[name]
  # Modes, for iterating over the rows. 
  modes = ('P1', 'P2', 'T1', 'T2')
  for row, mode in enumerate(modes):
    for col, filt in enumerate(filters):
      # Grab the event histogram, filtered appropriately. 
      eh = eventhist(hargs, mode=mode, **filt)
      # Indicate event count and overall rate in the corners. 
      eventcount = np.sum(eh)
      pct = tdp( 100.*np.mean(eh/z) ) + '\\%'
      count = znt(eventcount)
      PW[row, col].setParams(lcorner=count, rcorner=pct)
      # Add the mesh to the plot. 
      PW[row, col].setMesh(x, y, 100*zmask(eh)/z)
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'mode_' +  name + '.pdf')
  else:
    return PW.render()








# #############################################################################
# ############################################################ Event Histograms
# #############################################################################

# =============================================================================
# ================================================================== All Events
# =============================================================================

def allplot(storm=None, save=False):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  if storm is None:
    pos = getpos(**pargs)
  elif storm is True:
    pos = getpos(dst_lt=-30, **pargs)
  elif storm is False:
    pos = getpos(dst_ge=-30, **pargs)
  else:
    print 'Plotting storm and calm together looks real bad! '
    return
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Create a plot window to show different subsets of the events. 
  PW = plotWindow( **bep() )
  # Title and labels. 
  status = {True:'Storm ', False:'Quiet ', None:''}[storm]
#  title = notex( status + 'Pc4 Observation Rate: All Modes, All Phases, 0.01\\frac{mW}{m^2} and Larger')
  title = notex( status + 'Pc4 Observation Rate')
  PW.setParams(title=title)
  # Grab the events histogram. 
  sargs = { True:{'dst_lt':-30}, False:{'dst_ge':-30}, None:{} }[storm]
  eh = eventhist(hargs, **sargs)
  # Normalize it by the sampling rate and plot it. 
  PW.setMesh(x, y, 100*zmask(eh)/z)

  print 'Events in each L bin:'
  for erow in eh:
    print np.sum(erow)

  eventcount = np.sum(eh)
  pct = notex('Rate: ') + tdp( 100.*np.mean(eh/z) ) + '\\%'
  count = notex('Count: ') + znt(eventcount)
  PW.setParams(lcorner=count, rcorner=pct)

  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'rate_' + {True:'storm', False:'calm', None:'all'}[storm] + label + '.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================================================== All Events by Mode
# =============================================================================

def modeplot(storm=None, save=False):
  global pargs, plotdir

  # Set up the grid and 2D histogram based on probe position. 
  if storm is None:
    pos = getpos(**pargs)
  elif storm is True:
    pos = getpos(dst_lt=-30, **pargs)
  elif storm is False:
    pos = getpos(dst_ge=-30, **pargs)
  else:
    print 'Plotting storm and calm together looks real bad! '
    return

  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]

  # Create a plot window to show different subsets of the events. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  # Title and labels. 
  status = {True:'Storm ', False:'Quiet ', None:''}[storm]
#  title = notex(status + 'Pc4 Observation Rate by Mode: All Phases, 0.01\\frac{mW}{m^2} and Larger')
  title = notex(status + 'Pc4 Observation Rate by Mode')
  collabels = ( notex('Poloidal'), notex('Toroidal') )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title)
  # Iterate over the filters. 
  mfilt, hfilt = ('P', 'T'), ('1', '2')
  for row, hf in enumerate(hfilt):
    for col, mf in enumerate(mfilt):
      # Grab a histogram of the events, filtered by mode and harmonic. 
      sargs = { True:{'dst_lt':-30}, False:{'dst_ge':-30}, None:{} }[storm]
      eh = eventhist(hargs, mode=mf+hf, **sargs)

      eventcount = np.sum(eh)
      pct = tdp( 100.*np.mean(eh/z) ) + '\\%'
      count = znt(eventcount)
      PW[row, col].setParams(lcorner=count, rcorner=pct)

      # Add the mesh to the plot. Normalize by the sampling rate. 
      PW[row, col].setMesh(x, y, 100*zmask(eh)/z)
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'mode_' +  {True:'storm', False:'calm', None:'all'}[storm] + label + '.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================= Poloidal Compressional and Non-Compressional Events
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def azmplot(storm=None, save=False, split=0.2):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  if storm is None:
    pos = getpos(**pargs)
  elif storm is True:
    pos = getpos(dst_lt=-30, **pargs)
  elif storm is False:
    pos = getpos(dst_ge=-30, **pargs)
  else:
    print 'Plotting storm and calm together looks real bad! '
    return
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  status = {True:'Storm ', False:'Quiet ', None:''}[storm]
  title = notex(status + 'Poloidal Pc4 by Compressional Coupling: ' + titlehelper )
  rlabs = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )

  clabs = ( tex('BBp') + ' \\geq ' + format(split, '.1f'), tex('BBp') + ' < ' + format(split, '.1f') )

  PW.setParams(collabels=clabs, rowlabels=rlabs, title=title)
  # Iterate over the filters. 
  for row, mode in enumerate( ('P1', 'P2') ):
    # Grab a histogram of the events, filtered by spectral width. 

    sargs = { True:{'dst_lt':-30}, False:{'dst_ge':-30}, None:{} }[storm]

    smallm  = eventhist(hargs, mode=mode, comp_ge=split, **sargs)
    bigm = eventhist(hargs, mode=mode, comp_lt=split, **sargs)
    # Normalize by how long each region was sampled. 
    rates  = 100*zmask(smallm)/z, 100*zmask(bigm)/z

    smallmcount = np.sum(smallm)
    pct = tdf( 100.*np.mean(smallm/z) ) + '\\%'
    count = znt(smallmcount)
    PW[row, 0].setParams(lcorner=count, rcorner=pct)

    bigmcount = np.sum(bigm)
    pct = tdf( 100.*np.mean(bigm/z) ) + '\\%'
    count = znt(bigmcount)
    PW[row, 1].setParams(lcorner=count, rcorner=pct)

    # Add the mesh to the plot. 
    [ PW[row, i].setMesh(x, y, r) for i, r in enumerate(rates) ]

  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'azm_rate_' + {True:'storm', False:'calm', None:'all'}[storm] + label + '.pdf')
  else:
    return PW.render()
















# =============================================================================
# =============================================================== Double Events
# =============================================================================

def doubleplot(save=False, split=-30.):
  global pargs, plotdir, thresh
  # Create a plot window to show different subsets of the events. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  # Title and labels. 
  title = notex( 'Rate of Double Events by Parity and Storm Index' )

  rlabs = ( notex('Odd'), notex('Even') )

  clabs = ( notex('DST') + ' \\geq ' + znt(split) + notex('nT'), 
            notex('DST') + ' < ' + znt(split) + notex('nT') )

  PW.setParams(collabels=clabs, rowlabels=rlabs, title=title)

  for col, key in enumerate( ('dst_ge', 'dst_lt') ):
    # Set up the grid and 2D histogram based on probe position. 
    pos = getpos( **dict( pargs.items() + {key:split}.items() ) )
    x, y, z, hargs = [ pos[k] for k in ('x', 'y', 'z', 'hargs') ]

    # Odd and Even. 
    for row, harm in enumerate( ('1', '2') ):

      # Grab a histogram of appropriately-filtered double events. 
      dh = doublehist(hargs, pmode='P' + harm, tmode='T' + harm, **{key:split})

      eventcount = np.sum(dh)
      pct = tdp( 100.*np.mean(dh/z) ) + '\\%'
      count = znt(eventcount)
      PW[row, col].setParams(lcorner=count, rcorner=pct)

      # Build it into a histogram. Normalize based on sampling. 
      rate = 100*zmask(dh)/z
      PW[row, col].setMesh(x, y, rate)
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'double_rate' + label + '.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================================= Poloidal or Toroidal Events by FWHM
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def fwhmplot(mode, split=1., storm=None, save=False):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  if storm is None:
    pos = getpos(**pargs)
  elif storm is True:
    pos = getpos(dst_lt=-30, **pargs)
  elif storm is False:
    pos = getpos(dst_ge=-30, **pargs)
  else:
    print 'Plotting storm and calm together looks real bad! '
    return
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  status = {True:'Storm ', False:'Quiet ', None:''}[storm]
  modename = 'Poloidal' if mode=='p' else 'Toroidal'
  title = notex(status + modename + ' Pc4 by Spectral Width: ' + titlehelper )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( notex('FWHM') + ' \\geq ' + str(split) + notex('mHz'), notex('FWHM') + ' < ' + str(split) + notex('mHz') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title)
  # Iterate over the filters. 
  for row, mname in enumerate( ('P1', 'P2') if mode=='p' else ('T1', 'T2') ):
    # Grab a histogram of the events, filtered by spectral width. 
    sargs = { True:{'dst_lt':-30}, False:{'dst_ge':-30}, None:{} }[storm]

    broad = eventhist(hargs, mode=mname, fwhm_ge=split, **sargs)
    narrow  = eventhist(hargs, mode=mname, fwhm_lt=split, **sargs)

    broadcount = np.sum(broad)
    pct = tdf( 100.*np.mean(broad/z) ) + '\\%'
    count = znt(broadcount)
    PW[row, 0].setParams(lcorner=count, rcorner=pct)

    narrowcount = np.sum(narrow)
    pct = tdf( 100.*np.mean(narrow/z) ) + '\\%'
    count = znt(narrowcount)
    PW[row, 1].setParams(lcorner=count, rcorner=pct)

    # Normalize by how long each region was sampled. 
    rates  = 100*zmask(broad)/z, 100*zmask(narrow)/z
    # Add the mesh to the plot. 
    [ PW[row, i].setMesh(x, y, r) for i, r in enumerate(rates) ]
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'fwhm_rate_' + mode + '_' + {True:'storm', False:'calm', None:'all'}[storm] + '.pdf')
  else:
    return PW.render()

# =============================================================================
# ========================================== Poloidal or Toroidal Events by Dst
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def dstplot(mode, split=-30., save=False):
  global pargs, plotdir
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  modename = 'Poloidal' if mode=='p' else 'Toroidal'
  title = notex(modename + ' Pc4 by Dst: ' + titlehelper )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( notex('Dst') + ' \\geq ' + znt(split) + notex('nT'), 
                notex('Dst') + ' < ' + znt(split) + notex('nT') )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title)
  # Iterate over the filters. 
  for row, mname in enumerate( ('P1', 'P2') if mode=='p' else ('T1', 'T2') ):
    # Calm and storm. 
    for col, key in enumerate( ('dst_ge', 'dst_lt') ):
      # Set up the grid and 2D histogram based on probe position. 
      pos = getpos( **dict( pargs.items() + {key:split}.items() ) )
      x, y, z, hargs = [ pos[k] for k in ('x', 'y', 'z', 'hargs') ]
      # Get a histogram of events which are the desired mode and storm phase. 
      eh = eventhist( hargs, mode=mname, **{key:split} )

      eventcount = np.sum(eh)
      pct = tdf( 100.*np.mean(eh/z) ) + '\\%'
      count = znt(eventcount)
      PW[row, col].setParams(lcorner=count, rcorner=pct)

      # Normalize the event count based on sample time. 
      rate  = 100*zmask(eh)/z
      # Add the mesh to the plot. 
      PW[row, col].setMesh(x, y, rate)
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'dst_rate_' + mode + label + '.pdf')
  else:
    return PW.render()

# =============================================================================
# ====================================== Poloidal or Toroidal Events by L - LPP
# =============================================================================

# Let's take a look at a plot of how FWHM depends on mode. 
def llppplot(mode, split=0., save=False):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=2, **bep() )
  modename = 'Poloidal' if mode=='p' else 'Toroidal'
  title = notex(modename + ' Pc4 by Plasmapause Location: ' + titlehelper )
  rowlabels = ( notex('Odd\nHarmonic'), notex('Even\nHarmonic') )
  collabels = ( 'L - L_{PP} < ' + str(split), 'L - L_{PP} > ' + str(split) )
  PW.setParams(collabels=collabels, rowlabels=rowlabels, title=title)
  # Iterate over the filters. 
  for row, mname in enumerate( ('P1', 'P2') if mode=='p' else ('T1', 'T2') ):
    # Grab a histogram of the events, filtered by spectral width. 
    inside  = eventhist(hargs, mode=mname, llpp_lt=split)
    outside = eventhist(hargs, mode=mname, llpp_ge=split)

    print mname + ' overall inside  rate: ' + format(100*np.sum(inside)/np.sum(z), '.1f') + '%'
    print mname + ' overall outside rate: ' + format(100*np.sum(outside)/np.sum(z), '.1f') + '%'

    # Normalize by how long each region was sampled. 
    rates  = 100*zmask(inside)/z, 100*zmask(outside)/z
    # Add the mesh to the plot. 
    [ PW[row, i].setMesh(x, y, r) for i, r in enumerate(rates) ]
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'llpp_rate_' + mode + '.pdf')
  else:
    return PW.render()

# =============================================================================
# ======================================================= All Events by L - LPP
# =============================================================================

def ppplot(split=0., save=False):
  global pargs, plotdir
  # Set up the grid and 2D histogram based on probe position. 
  pos = getpos(**pargs)
  x, y, z, hargs = [ pos[key] for key in ('x', 'y', 'z', 'hargs') ]
  # Set up the window. 
  PW = plotWindow( ncols=2, nrows=1, **bep() )
  title = notex('Pc4 Rate by Plasmapause Location: ' + titlehelper )
  collabels = ( 'L - L_{PP} < ' + str(split), 'L - L_{PP} > ' + str(split) )
  PW.setParams(collabels=collabels, title=title)
  # Grab a histogram of the events, filtered by L-LPP. 
  inside  = eventhist(hargs, llpp_lt=split)
  outside = eventhist(hargs, llpp_ge=split)
  # Normalize by how long each region was sampled. 
  rates  = 100*zmask(inside)/z, 100*zmask(outside)/z
  # Add the mesh to the plot. 
  [ PW[i].setMesh(x, y, r) for i, r in enumerate(rates) ]
  # Show or save the plot. 
  if save is True:
    return PW.render(plotdir + 'llpp' + label + '.pdf')
  else:
    return PW.render()

# #############################################################################
# ################################################################## Data Input
# #############################################################################

# =============================================================================
# ==================================================================== Position
# =============================================================================

def getpos(dl=0.5, dm=1, lmin=None, lmax=None, dst_ge=None, dst_lt=None):
  # Scale from five-minute chunks to half-hour events. 
  secondsper = 1800.
  # The orbit of both RBSP paths has been broken into five-minute chunks. Grab
  # the position of each chunk which gives good data for E dot B = 0. 
  poslines = g2a( line for line in read('pos.txt') if 'ok' in line )
  # If we only want the position when Dst is above a certain threshold. 
  if dst_ge is not None:
    dst = g2a( int( line.split()[7] ) for line in poslines )
    inew = np.nonzero(dst >= dst_ge)[0]
    poslines = poslines[inew]
  # If we only want the position when Dst is below a certain threshold.
  if dst_lt is not None:
    dst = g2a( int( line.split()[7] ) for line in poslines )
    inew = np.nonzero(dst < dst_lt)[0]
    poslines = poslines[inew]
  # Get the date range. 
  dates = ( poslines[0].split()[1], poslines[-1].split()[1] )
  # Arrange the positions as an array of floats. 
  pos = np.array( [ [ float(x) for x in p.split()[3:6] ] for p in poslines ] )
  # Figure out the histogram bounds. 
  if lmin is None:
    lmin = np.floor( np.min( pos[:, 0] ) )
  if lmax is None:
    lmax = np.ceil( np.max( pos[:, 0] ) )
  # Center MLT bins on the hour, at least at midnight. 
  mmin, mmax = -dm/2., 24 - dm/2.
  # We want a bin to be centered at zero. That means anything between (24-dm/2)
  # and 24 should be mapped to the range (-dm/2) to 0. 
  posm = np.where( pos[:, 1] > mmax, pos[:, 1] - 24, pos[:, 1] )
  # Number of bins in each direction. 
  lbins, mbins = int( (lmax - lmin)/dl ) + 1, int( (mmax - mmin)/dm ) + 1
  # Keyword arguments for the histogram2d call. 
  hargs = { 'range':( (lmin, lmax), (mmin, mmax) ), 'bins':(lbins-1, mbins-1) }
  # Bin bounds in terms of L and MLT. 
  l, m = np.mgrid[lmin:lmax:lbins*1j, mmin:mmax:mbins*1j]
  # Map to GSE coordinates. Put midnight at the bottom. 
  x, y = -l*np.sin(2*pi*m/24.), -l*np.cos(2*pi*m/24.)
  # Sometimes we don't want our radial binning to depend on L, but rather on
  # L - LPP. Map inside and outside the plasmapause to lmin and lmax.  
  posl = pos[:, 0]
  if 'lpp' in argv: 
    poslpp = np.array( [ float( p.split()[8] ) for p in poslines ] )
    posl = np.where(posl < poslpp, lmin, lmax)
  # Bin the position data into a 2D histogram, then scale it to days. 
  h = np.histogram2d(posl, posm, **hargs)[0]
  z = 300*h/secondsper
  # Return the position data. Total amount of usable time too. 
  return {'dates':dates, 'l':l, 'm':m, 'x':x, 'y':y, 'z':z, 'hargs':hargs}

# =============================================================================
# ================================================== Loading and Binning Events
# =============================================================================

# Returns a filtered list of events. 
def loadevents(mode=None, fwhm_ge=None, fwhm_lt=3., amp_ge=None, amp_lt=None, f_ge=None, f_lt=None, comp_ge=None, comp_lt=None, double=False, llpp_lt=None, llpp_ge=None, dst_lt=None, dst_ge=None, phase=None, antiphase=None):
  global thresh
  # Grab the events file as an array of strings. Skip the header. 
  events = g2a( line for line in read('events.txt') if 'probe' not in line )
  # Filter for simultaneous events. 
  if double is True:
    events = g2a( line for line in events if 'BIG' in line or 'SMALL' in line )
  # Filter based on mode. 
  if mode is not None:
    events = g2a( line for line in events if mode in line )
  # Filter on amplitude. 
  if thresh > 0:
    amp = g2a( float( line.split()[11] ) for line in events )
    inew = np.nonzero(amp >= thresh)[0]
    events = events[inew]

  # For phase, we look at how close the absolute value is to 90 degrees. A cutoff of 45 means anything with an absolute value between 45 degrees and 135 degrees is fair game. 
  if phase is not None:
    ph = g2a( float( line.split()[10] ) for line in events )
    inew = np.nonzero( np.logical_and(np.abs(ph) > phase, np.abs(ph) < 180 - phase) )
    events = events[inew]
  if antiphase is not None:
    ph = g2a( float( line.split()[10] ) for line in events )
    inew = np.nonzero( np.logical_or(np.abs(ph) < antiphase, np.abs(ph) > 180 - antiphase) )
    events = events[inew]

  # Filter on frequency. 
  if f_ge is not None:
    f = g2a( float( line.split()[8] ) for line in events )
    inew = np.nonzero(f >= f_ge)[0]
    events = events[inew]
  if f_lt is not None:
    f = g2a( float( line.split()[8] ) for line in events )
    inew = np.nonzero(f < f_lt)[0]
    events = events[inew]


  # Filter on amplitude. 
  if amp_ge is not None:
    amp = g2a( float( line.split()[11] ) for line in events )
    inew = np.nonzero(amp >= amp_ge)[0]
    events = events[inew]
  if amp_lt is not None:
    amp = g2a( float( line.split()[11] ) for line in events )
    inew = np.nonzero(amp < amp_lt)[0]
    events = events[inew]


  # Filter on compressional coupling (lower bound). 
  if comp_ge is not None:
    comp = g2a( float( line.split()[12] ) for line in events )
    inew = np.nonzero(comp >= comp_ge)[0]
    events = events[inew]
  # Filter on compressional coupling (upper bound). 
  if comp_lt is not None:
    comp = g2a( float( line.split()[12] ) for line in events )
    inew = np.nonzero(comp < comp_lt)[0]
    events = events[inew]
  # Filter on spectral width (lower bound). 
  if fwhm_ge is not None:
    fwhm = g2a( float( line.split()[9] ) for line in events )
    inew = np.nonzero(fwhm >= fwhm_ge)[0]
    events = events[inew]
  # Filter on compressional coupling (upper bound). 
  if fwhm_lt is not None:
    fwhm = g2a( float( line.split()[9] ) for line in events )
    inew = np.nonzero(fwhm < fwhm_lt)[0]
    events = events[inew]
  # Filter on position relative to the plasmapause (lower bound). 
  if llpp_ge is not None:
    lpp = g2a( float( line.split()[6] ) for line in events )
    lshell = g2a( float( line.split()[3] ) for line in events )
    inew = np.nonzero(lshell - lpp >= llpp_ge)[0]
    events = events[inew]
  # Filter on position relative to the plasmapause (upper bound). 
  if llpp_lt is not None:
    lpp = g2a( float( line.split()[6] ) for line in events )
    lshell = g2a( float( line.split()[3] ) for line in events )
    inew = np.nonzero(lshell - lpp < llpp_lt)[0]
    events = events[inew]
  # Filter on storm index (upper bound). 
  if dst_ge is not None:
    dst = g2a( float( line.split()[13] ) for line in events )
    inew = np.nonzero(dst >= dst_ge)[0]
    events = events[inew]
  # Filter on storm index (upper bound). 
  if dst_lt is not None:
    dst = g2a( float( line.split()[13] ) for line in events )
    inew = np.nonzero(dst < dst_lt)[0]
    events = events[inew]
  # If we're not looking for double events, and we're not filtering based on
  # mode, make sure not to double-count double events. For each big event, make
  # sure there's no small event at that same timestamp. 
  if mode is None and double is False:
    bigstamps = g2a( e[:26] for e in events if 'BIG' in e )
    inew = np.nonzero( g2a( 'SMALL' not in e or e[:26] not in bigstamps for e in events ) )[0]
    events = events[inew]
  # Return the remaining events. 
  return events

# Returns a histogram of filtered events. Each count is half an hour. 
def eventhist(hargs, **kargs):
  # Grab the events. 
  events = loadevents(**kargs)
  # If there are no events, return an empty histogram. 
  if events.size == 0:
    return np.histogram2d([], [], **hargs)[0]
  # Get the position from each. 
  pos = g2a( [ float(x) for x in line.split()[3:6] ] for line in events )
  lshell, mlt, mlat = pos[:, 0], pos[:, 1], pos[:, 2]

  # Sometimes we want to use L-LPP for the radial coordinate. 
  if 'lpp' in argv:
    poslpp = g2a( float( line.split()[6] ) for line in events )
    lshell = np.where( lshell < poslpp, hargs['range'][0][0], hargs['range'][0][1] )

  # Shift mlt a bit... we want our first bin centered at midnight. 
  mltmax = hargs['range'][-1][-1]
  mlt = np.where(mlt>mltmax, mlt - 24, mlt)
  # Assemble the positions of these events into a histogram, using the bins
  # defined by hargs. 
  return np.histogram2d(lshell, mlt, **hargs)[0]

# Returns a histogram of double events: those that trigger the poloidal and 
# toroidal channels simultaneously. 
def doublehist(hargs, pmode=None, tmode=None, **kargs):
  # Grab the poloidal and toroidal events. 
  pev = loadevents(mode=pmode, double=True, **kargs)
  tev = loadevents(mode=tmode, double=True, **kargs)
  # The number of events should be pretty small, so just find matches by brute
  # force. The first 26 characters of the line give probe, date, time. 
  doubles = g2a( p + '\n' + t for p in pev for t in tev if p[:26]==t[:26] )
  # If there are no events, return an empty histogram. 
  if doubles.size == 0:
    return np.histogram2d([], [], **hargs)[0]

  dates = g2a( d[3:13] for d in doubles )

#  for d in sorted( set(dates) ):
#    print '\t' + d + '\t', len( list(p for p in doubles if d in p ) )


  # Tally up the days. 
  print pmode, tmode, kargs
  print '\tnumber of events: ', dates.size
  print '\tnumber of dates:  ', g2a( set(dates) ).size

  # Get the position for each event, then return a histogram of those events,
  # scaled to units of days. 
  pos = g2a( [ float(x) for x in line.split()[3:6] ] for line in doubles )

  lshell, mlt, mlat = pos[:, 0], pos[:, 1], pos[:, 2]
  mltmax = hargs['range'][-1][-1]
  mlt = np.where(mlt>mltmax, mlt - 24, mlt)

  hist = np.histogram2d(lshell, mlt, **hargs)[0]

  print '\tnumber of points in histogram: ', np.sum(hist)



  return np.histogram2d(lshell, mlt, **hargs)[0]

# #############################################################################
# #################################### Outer Magnetosphere and the Dungey Cycle
# #############################################################################

'''
def gauss(x, amp=0, avg=np.pi/2, std=1):
  return amp*np.exp( -(x - avg)**2/(2.*std**2) )

def lshell(L, stretch=1., qstretch=0., *args, **kargs):
  q0 = np.arcsin( np.sqrt( 1./np.abs(L) ) )
  q = np.linspace(q0, np.pi - q0, 100)
  r = L*np.sin(q)**2 + gauss(q, *args, **kargs)
  x, z = r*np.sin(q), r*np.cos(q)
  dx = x - x[0]
  xnew = x[0] + stretch*dx + qstretch*dx**2
  return xnew, z

def straighten(x, y):
  n = x.size
  xp, yp = x, y  
  x0, y0 = x[n/2-1], y[n/2-1]
  xp[n/2:] = np.linspace(x0, 30*np.sign(x0), n/2)
  yp[n/2:] = y0
  return xp, yp

def q(L):
  q0 = np.arcsin( np.sqrt( 1./np.abs(L) ) )
  return np.linspace(q0, np.pi - q0, 100)
'''

# L-shell, with a quadratic stretch. 
def lshell(L, qs=0.):
  x = np.linspace(-20, 20, 1000)
  x0 = np.sqrt( np.abs(1./L) ) if L!=0 else 0
  if L>0:
    xp = np.where( x<x0, x, x + qs*(x - x0)**2  )
  else:
    xp = np.where( x>x0, x, x - qs*(x - x0)**2  )
  zz = np.where( np.sign(L)==np.sign(xp), np.abs(L*xp**2)**(2./3) - xp**2, 0)
  z = np.sqrt( np.maximum(zz, 0 ) )
  return x, np.where( x**2 + z**2 > 1, z, 0)

def openline(x0):
  x = np.linspace(-20, 20, 1000)
  c = 0.1*(20 - x0)**2 / np.sqrt(20 - x0)
  if x0 < 0: 
    z = 100*c*np.sqrt( np.maximum(0, x0 - x) )
  else: 
    z = c*np.sqrt( np.maximum(x - x0, 0) )
  return x, z


def dungey():
  PW = plotWindow(ncols=-2, colorbar=None)
  PW.setParams( xlims=(-20, 20), ylims=(-8, 8) )
  PW.setParams(earth='left')
  dstretch = 1.2e-3
  nstretch = 1.2e-3
  [ PW.setLine( *lshell(L, qs=-nstretch*L), color='b' ) for L in range(2, 12, 4) ]
  [ PW.setLine( *lshell(L, qs=-dstretch*L), color='b' ) for L in range(-10, 0, 4) ]
  [ PW.setLine( *openline(x0), color='r' ) for x0 in (-18, -14, 16, 18) ]
  for L, x0 in ( (-12, -10), (-20, -6), (-32, -1.5), (12, 14), (18, 10), (32, 2) ):
    x, zc = lshell(L, qs=-dstretch*L)
    x, zo = openline(x0)
    p = 3.
    PW.setLine(x, (zo**p + zc**p)**(1/p), 'm')
  return PW.render()

  '''
  for i, L in enumerate( (2, 4, 6, 8, 10) ):
    PW.setLine( *lshell(-L), color='b' )
    PW.setLine( *lshell(-L, qstretch=i/200.), color='r' )
    PW.setLine( *lshell(L), color='b' )
    x, z = straighten( *lshell(L, qstretch=i/100., std=0.1, amp=i/10.) )
    PW.setLine( x, z, color='r' )
#  for i, x0 in enumerate( range(0, 20, 2) ):
#    z = (10 - i)*np.sqrt( np.maximum(x - x0, 0)/(20. - x0) )
#    PW.setLine(x, z, 'g')
  z = np.linspace(-8, 8, 100)
  for x0 in range(-20, 20, 2):
    x = x0 + 10*(20 - x0)**-2 * z**2
    PW.setLine(x, z, 'g')
  return PW.render()
#  PW.setLine( *lshell(4, amp=0.5, std=1), color='b' )
#  PW.setLine( *lshell(6, amp=1, std=0.5), color='b' )
#  PW.setLine( *lshell(8,  amp=1, std=0.1), color='b' )
#  PW.setLine( *lshell(10,  amp=2, std=0.1), color='b' )
  x, y0 = straighten( *lshell(12, amp=2, std=0.1) )
  PW.setLine( x, y0, color='m' )
#  y1 = np.sqrt( np.maximum(x - 10, 0) )
#  PW.setLine( x, y1, color='g' )
#  PW.setLine( x, np.sqrt(y0**2 + y1**2), color='orange' )
  y1 = 2*np.sqrt( np.maximum(x - 10, 0)/10. )
  PW.setLine(x, y1, 'r')
  PW.setLine(x, np.sqrt(y0**2 + y1**2), 'g')
  y2 = 1*np.sqrt( np.maximum(x - 14, 0)/6. )
  PW.setLine(x, y2, 'r')
#  for i, x0 in enumerate( (16, 14, 12, 10, 8) ):
#    y = (i+1)*np.sqrt( np.maximum(x - x0, 0)/(20. - x0) )
#    PW.setLine(x, y, 'r')
  PW.render()
  '''
  return

# #############################################################################
# ############################################################ Helper Functions
# #############################################################################

# Turns a generator expression into a list, if necessary, then turns the list
# into an array. 
def g2a(expr):
  return np.array( list(expr) )

# Bullseye params for the plotter. 
def bep(rate=True):
  tls = ('$-8$', '', '$-4$', '', '$0$', '', '$+4$', '', '$+8$')
  colorbar = 'pct' if rate is True else 'pos' 

  ncolors = 11 if rate is True else 12

  return {'earth':'top', 'flipx':True, 'grid':True, 'square':True, 
          'xlabel': 'Y' + notex(' (R_E)'), 'xlims':(-8, 8),
          'xticks':np.mgrid[-8:8:9j], 'xticklabels':tls, 
          'ylabel': 'X' + notex(' (R_E)'), 'ylims':(-8, 8), 
          'yticks':np.mgrid[-8:8:9j], 'yticklabels':tls, 'ylabelpad':-2, 
          'colorbar':colorbar, 'ncolors':ncolors }

# Space out a nice right-justified column. 
def col(x, width=12, unit='', digs=2):
  d = str(digs)
  if isinstance(x, float):
    return format(x, str(width - len(unit) - 1) + '.' + d + 'f') + unit + ' '
  else:
    return str(x).rjust(width - len(unit) - 1) + unit + ' '

# #############################################################################
# ########################################################### For Importability
# #############################################################################

if __name__=='__main__':
  main()


