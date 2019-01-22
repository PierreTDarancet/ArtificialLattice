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

_p = _common.lazy_import('_plotter')

def _make_figure(dpi, fig_size, use_pyplot=False):
    if 'matplotlib.backends' not in sys.modules:
        warnings.warn(
            "Kwant's plotting functions have\nthe side effect of "
            "selecting the matplotlib backend. To avoid this "
            "warning,\nimport matplotlib.pyplot, "
            "matplotlib.backends or call matplotlib.use().",
            RuntimeWarning, stacklevel=3
        )
    if use_pyplot:
        # We import backends and pyplot only at the last possible moment (=now)
        # because this has the side effect of selecting the matplotlib backend
        # for good.  Warn if backend has not been set yet.  This check is the
        # same as the one performed inside matplotlib.use.
        from matplotlib import pyplot
        fig = pyplot.figure()
    else:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        fig = _p.Figure()
        fig.canvas = FigureCanvasAgg(fig)
    if dpi is not None:
        fig.set_dpi(dpi)
    if fig_size is not None:
        fig.set_figwidth(fig_size[0])
        fig.set_figheight(fig_size[1])
    return fig

def _maybe_output_fig(fig, file=None, show=True):
    """Output a matplotlib figure using a given output mode.

    Parameters
    ----------
    fig : matplotlib.figure.Figure instance
        The figure to be output.
    file : string or a file object
        The name of the target file or the target file itself
        (opened for writing).
    show : bool
        Whether to call ``matplotlib.pyplot.show()``.  Only has an effect if
        not saving to a file.

    Notes
    -----
    The behavior of this function producing a file is different from that of
    matplotlib in that the `dpi` attribute of the figure is used by defaul
    instead of the matplotlib config setting.
    """
    if fig is None:
        return

    if file is not None:
        fig.canvas.print_figure(file, dpi=fig.dpi)
    elif show:
        # If there was no file provided, pyplot should already be available and
        # we can import it safely without additional warnings.
        from matplotlib import pyplot
        pyplot.show()

def zak_bands(sys, args=(), momenta=65, file=None, show=True, dpi=None,
                fig_size=None, ax=None, *, params=None):
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

    if not _p.mpl_available:
        raise RuntimeError("matplotlib was not found, but is required "
                           "for bands()")

    syst = sys  # for naming consistency inside function bodies
    momenta = np.array(momenta)
    if momenta.ndim != 1:
        momenta = np.linspace(-np.pi, np.pi, momenta)

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
        mat = hop * cmath.exp(-1j * k)
        mat +=  mat.conjugate().transpose() + ham
        return mat

    return spectrum(h_k, ('k', momenta), file=file, show=show, dpi=dpi,
                    fig_size=fig_size, ax=ax)


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

    if not _p.mpl_available:
        raise RuntimeError("matplotlib was not found, but is required "
                           "for plot_spectrum()")
    if y is not None and not _p.has3d:
        raise RuntimeError("Installed matplotlib does not support 3d plotting")

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

    # set up axes
    if ax is None:
        fig = _make_figure(dpi, fig_size, use_pyplot=(file is None))
        if y is None:
            ax = fig.add_subplot(1, 1, 1)
        else:
            warnings.filterwarnings('ignore',
                                    message=r'.*mouse rotation disabled.*')
            ax = fig.add_subplot(1, 1, 1, projection='3d')
            warnings.resetwarnings()
        ax.set_xlabel(keys[0])
        if y is None:
            ax.set_ylabel('Energy')
        else:
            ax.set_ylabel(keys[1])
            ax.set_zlabel('Energy')
        ax.set_title(', '.join('{} = {}'.format(*kv) for kv in params.items()))
    else:
        fig = None

    # actually do the plot
    if y is None:
        ax.plot(array_values[0], spectrum)
    else:
        if not hasattr(ax, 'plot_surface'):
            msg = ("When providing an axis for plotting over a 2D domain the "
                   "axis should be created with 'projection=\"3d\"")
            raise TypeError(msg)
        # plot_surface cannot directly handle rank-3 values, so we
        # explicitly loop over the last axis
        grid = np.meshgrid(*array_values)
        for i in range(spectrum.shape[-1]):
            spec = spectrum[:, :, i].transpose()  # row-major to x-y ordering
            ax.plot_surface(*(grid + [spec]), cstride=1, rstride=1)

    _maybe_output_fig(fig, file=file, show=show)
    return [fig,wf]
