# -*- coding: utf-8 -*-
# Copyright 2011-2018 Kwant authors.
#
# This file is modified starting from the plotter.py provided as a  part of Kwant.  

# It is subject to the license terms in the file
# LICENSE.rst found in the top-level directory of this distribution and at
# http://kwant-project.org/license.  A list of Kwant authors can be found in
# the file AUTHORS.rst at the top-level directory of this distribution and at
# http://kwant-project.org/authors.

""" Zak compute Module 

This module contains some helped functions to get the eigenvectors of a 
1D periodic system 
"""

import sys
import itertools
import functools
import warnings
import cmath
import numpy as np
import tinyarray as ta
from scipy import spatial, interpolate
from math import cos, sin, pi, sqrt
from kwant import system,_common 


def zak_bands(sys, args=(), momenta=65, file=None, show=True, dpi=None,
                fig_size=None, ax=None, *, params=None,mom_start=0):
    """Plot band structure of a translationally invariant 1D system.

    Parameters
    ----------
    sys : kwant.system.InfiniteSystem
        A system bands of which are to be plotted.
    args : tuple, defaults to empty
        Positional arguments to pass to the ``hamiltonian`` method.
        Mutally exclusive with 'params'.
    momenta : int or 1D array-like
        Either a number of sampling points on the interval [-pi, pi], or an
        array of points at which the band structure has to be evaluated.
    file : string or file object or `None`
        The output file.  If `None`, output will be shown instead.
    show : bool
        Whether ``matplotlib.pyplot.show()`` is to be called, and the output is
        to be shown immediately.  Defaults to `True`.
    dpi : float
        Number of pixels per inch.  If not set the ``matplotlib`` default is
        used.
    fig_size : tuple
        Figure size `(width, height)` in inches.  If not set, the default
        ``matplotlib`` value is used.
    ax : ``matplotlib.axes.Axes`` instance or `None`
        If `ax` is not `None`, no new figure is created, but the plot is done
        within the existing Axes `ax`. in this case, `file`, `show`, `dpi`
        and `fig_size` are ignored.
    params : dict, optional
        Dictionary of parameter names and their values. Mutually exclusive
        with 'args'.

    Returns
    -------
    fig : matplotlib figure
        A figure with the output if `ax` is not set, else None.

    Notes
    -----
    See `~kwant.physics.Bands` for the calculation of dispersion without plotting.
    """


    syst = sys  # for naming consistency inside function bodies
    momenta = np.array(momenta)
    if momenta.ndim != 1:
        #momenta = np.linspace(-np.pi, np.pi, momenta)
        momenta = np.linspace(mom_start, np.pi, momenta)

    # expand out the contents of 'physics.Bands' to get the H(k),
    # because 'spectrum' already does the diagonalisation.
    ham = syst.cell_hamiltonian(args, params=params)
    if not np.allclose(ham, ham.conjugate().transpose()):
        raise ValueError('The cell Hamiltonian is not Hermitian.')
    _hop = syst.inter_cell_hopping(args, params=params)
    hop = np.empty(ham.shape, dtype=complex)
    hop[:, :_hop.shape[1]] = _hop
    hop[:, _hop.shape[1]:] = 0

    def h_k(k):
        # H_k = H_0 + V e^-ik + V^\dagger e^ik
        mat = hop * cmath.exp(-1j * np.array(k))
        mat +=  mat.conjugate().transpose() + ham
        return mat

    return spectrum(h_k, ('k', momenta), file=file, show=show, dpi=dpi,
                    fig_size=fig_size, ax=ax)
    #return h_k


def spectrum(syst, x, y=None, params=None, mask=None, file=None,
             show=True, dpi=None, fig_size=None, ax=None):
    """Plot the spectrum of a Hamiltonian as a function of 1 or 2 parameters

    Parameters
    ----------
    syst : `kwant.system.FiniteSystem` or callable
        If a function, then it must take named parameters and return the
        Hamiltonian as a dense matrix.
    x : pair ``(name, values)``
        Parameter to ``ham`` that will be varied. Consists of the
        parameter name, and a sequence of parameter values.
    y : pair ``(name, values)``, optional
        Used for 3D plots (same as ``x``). If provided, then the cartesian
        product of the ``x`` values and these values will be used as a grid
        over which to evaluate the spectrum.
    params : dict, optional
        The rest of the parameters to ``ham``, which will be kept constant.
    mask : callable, optional
        Takes the parameters specified by ``x`` and ``y`` and returns True
        if the spectrum should not be calculated for the given parameter
        values.
    file : string or file object or `None`
        The output file.  If `None`, output will be shown instead.
    show : bool
        Whether ``matplotlib.pyplot.show()`` is to be called, and the output is
        to be shown immediately.  Defaults to `True`.
    dpi : float
        Number of pixels per inch.  If not set the ``matplotlib`` default is
        used.
    fig_size : tuple
        Figure size `(width, height)` in inches.  If not set, the default
        ``matplotlib`` value is used.
    ax : ``matplotlib.axes.Axes`` instance or `None`
        If `ax` is not `None`, no new figure is created, but the plot is done
        within the existing Axes `ax`. in this case, `file`, `show`, `dpi`
        and `fig_size` are ignored.

    Returns
    -------
    fig : matplotlib figure
        A figure with the output if `ax` is not set, else None.
    """


    if isinstance(syst, system.FiniteSystem):
        def ham(**kwargs):
            return syst.hamiltonian_submatrix(params=kwargs, sparse=False)
    elif callable(syst):
        ham = syst
    else:
        raise TypeError("Expected 'syst' to be a finite Kwant system "
                        "or a function.")

    params = params or dict()
    keys = (x[0],) if y is None else (x[0], y[0])
    array_values = (x[1],) if y is None else (x[1], y[1])

    # calculate spectrum on the grid of points
    spectrum = []
    wf = []
    bound_ham = functools.partial(ham, **params)
    for point in itertools.product(*array_values):
        p = dict(zip(keys, point))
        if mask and mask(**p):
            spectrum.append(None)
        else:
            h_p = np.atleast_2d(bound_ham(**p))
            energies,eig_vec = np.linalg.eigh(h_p)
            spectrum.append(energies)
            wf.append(eig_vec)
    # massage masked grid points into a list of NaNs of the appropriate length
    n_eigvals = len(next(filter(lambda s: s is not None, spectrum)))
    nan_list = [np.nan] * n_eigvals
    spectrum = [nan_list if s is None else s for s in spectrum]
    # make into a numpy array and reshape
    new_shape = [len(v) for v in array_values] + [-1]
    spectrum = np.array(spectrum).reshape(new_shape)
    return spectrum,wf
