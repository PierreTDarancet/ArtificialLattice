# -*- coding: utf-8 -*-
# Copyright 2011-2018 Kwant authors.
#
# This file is modified starting from the plotter.py provided as a  part of Kwant.  

# The bands method is rewritten to return the hamiltonian as a function of k


# It is subject to the license terms in the file
# LICENSE.rst found in the top-level directory of this distribution and at
# http://kwant-project.org/license.  A list of Kwant authors can be found in
# the file AUTHORS.rst at the top-level directory of this distribution and at
# http://kwant-project.org/authors.

""" Zak compute Module 

This module contains 
"""
from utilities import get_width
import cmath
import numpy as np
import z2pack


def get_Hk(sys, args=(), momenta=65, file=None, *, params=None,dim=3):
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
    params : dict, optional
        Dictionary of parameter names and their values. Mutually exclusive
        with 'args'.
    dim: int , optional
        Dimension of the system
        Default = 3 

    Returns
    -------
        h_k: A callable function defition which returns the hamiltonian 
            as a function of 'k'
    """


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
        mat = hop * cmath.exp(-1j * np.array(2*np.pi*k[0]))
        mat +=  mat.conjugate().transpose() + ham
        return mat

    #return spectrum(h_k, ('k', momenta), file=file, show=show, dpi=dpi,
    #                fig_size=fig_size, ax=ax)
    return h_k


def calc_pol(syst,red_pos):   
    
    Hk = get_Hk(syst,dim=2)
    z2_system = z2pack.hm.System(Hk,dim=2,pos=red_pos,
                                 convention=1)
    result = z2pack.line.run(system=z2_system, 
                              line=lambda t1: [t1,0],
                              pos_tol=1e-3,iterator=range(3,501,2));
    return result.pol
    
    
    
