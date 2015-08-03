# plotting.py
# Jacob Hummel
"""
This module contains classes for various types of plots.
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

def _hist2d(ax, x, y, **kwargs):
    N = kwargs.pop('gridsize',500)
    binning = kwargs.pop('binning','log')
    xscale = kwargs.pop('xscale','log')
    yscale = kwargs.pop('yscale','log')
    xlims = kwargs.pop('xlims', (-3.5,12))
    cmap = kwargs.pop('cmap', plt.cm.Blues_r)
    if binning == 'log':
        norm = LogNorm()
    else:
        norm=None
    if xscale == 'log':
        xbins = np.logspace(np.log10(x.min()), np.log10(x.max()), N)
    elif xscale == 'linear':
        xbins = np.linspace(x.min(), x.max(), N)
    else:
        raise KeyError
    if yscale == 'log':
        ybins = np.logspace(np.log10(y.min()), np.log10(y.max()), N)
    elif yscale == 'linear':
        ybins = np.linspace(y.min(), x.max(), N)
    else:
        raise KeyError
    heatmap, xedges, yedges = np.histogram2d(y, x, bins=(ybins,xbins))
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    ax.imshow(heatmap, origin='lower', norm=norm, extent=extent, cmap=cmap)
    ax.set(xscale=xscale, yscale=yscale, aspect='auto')
    return ax

def temp(snapshot, ax, **kwargs):
    dens = snapshot.gas.get_number_density()
    temp = snapshot.gas.get_temperature()
    select = kwargs.pop('select', None)
    cmbline = kwargs.pop('cmbline', True)
    if select is not None:
        dens = dens[select]
        temp = temp[select]

    ax = _hist2d(ax, dens, temp, **kwargs)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, 1e12)
    ax.set_ylim(10, 2e4)

    if cmbline:
        ax.axhline(2.725 * (snapshot.header.Redshift + 1),
                   linestyle='--', color='k')
    ax.set_xlabel('n [cm$^{-3}$]')
    ax.set_ylabel('Temperature [K]')
    plt.draw()
    return ax

def radial_temp(snapshot, ax, **kwargs):
    snapshot.gas.calculate_spherical_coords(unit='pc', centering='avg')
    r = snapshot.gas.spherical_coords[:,0]
    temp = snapshot.gas.get_temperature()
    select = kwargs.pop('select', None)
    if select is not None:
        r = r[select]
        temp = temp[select]

    virialized = numpy.where(r <= 70.)[0]
    r = r[virialized]
    temp = temp[virialized]
    snapshot.gas.cleanup()

    gridsize = kwargs.get('gridsize', None)
    if gridsize is None:
        kwargs['gridsize'] = 250
    ax = _hist2d(ax, r, temp, **kwargs)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(5e-5,70)
    ax.set_ylim(10, 2e4)
    ax.axhline(2.725 * (snapshot.header.Redshift + 1),
               linestyle='--', color='k')
    ax.set_xlabel('Radius [pc]')
    ax.set_ylabel('Temperature [K]')
    return ax

def electron_frac(snapshot, ax, **kwargs):
    dens = snapshot.gas.get_number_density()
    efrac = snapshot.gas.get_electron_fraction()
    select = kwargs.pop('select', None)
    if select is not None:
        dens = dens[select]
        efrac = efrac[select]

    ax = _hist2d(ax, dens, efrac, **kwargs)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, 1e12)
    ax.set_ylim(2e-11, 5e-2)
    ax.set_xlabel('n [cm$^{-3}$]')
    ax.set_ylabel('f$_{e^-}$')
    plt.draw()
    return ax

def h2frac(snapshot, ax, **kwargs):
    dens = snapshot.gas.get_number_density()
    h2frac = snapshot.gas.get_H2_fraction()
    select = kwargs.pop('select', None)
    if select is not None:
        dens = dens[select]
        h2frac = h2frac[select]

    ax = _hist2d(ax, dens, h2frac, **kwargs)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, 1e12)
    ax.set_ylim(5e-7,2)
    ax.set_xlabel('n [cm$^{-3}$]')
    ax.set_ylabel('f$_{H_2}$')
    plt.draw()
    return ax

def HDfrac(snapshot, ax, **kwargs):
    dens = snapshot.gas.get_number_density()
    HDfrac = snapshot.gas.get_HD_fraction()
    select = kwargs.pop('select', None)
    if select is not None:
        dens = dens[select]
        HDfrac = HDfrac[select]

    ax = _hist2d(ax, dens,HDfrac, **kwargs)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, 1e12)
    ax.set_ylim(5e-11,1e-4)
    ax.set_xlabel('n [cm$^{-3}$]')
    ax.set_ylabel('f$_{HD}$')
    plt.draw()
    return ax