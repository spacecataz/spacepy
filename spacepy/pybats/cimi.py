'''
CIMI - A module for interfacing with output from the Comprehensive 
Inner Magnetosphere Ionosphere Model.

.. currentmodule:: spacepy.pybats.cimi

.. rubric:: Classes

.. autosummary::
    :template: clean_class.rst
    :toctree:

    CimiIono
    CimiEq

'''

import numpy as np

from spacepy.plot import set_target
from spacepy.pybats import IdlFile

def parse_cimihead(headstring):
    '''
    Given a CIMI-style header line that reads `Var1[Units1] Var2[Units2] ...` return
    a list of variables and a list of units.
    '''

    pass

class CimiIono(IdlFile):
    pass

class CimiEq(IdlFile):
    '''
    A class to handle 2D equatorial CIMI output files.
    '''
    def __repr__(self):
        return 'CimiEq single output file %s' % (self.attrs['file'])

    # Init by calling IdlFile init and then building qotree, etc.
    def __init__(self, filename, *args, **kwargs):
        
        # Read file.
        IdlFile.__init__(self, filename, header=None, gtype="Irregular",
                         *args, **kwargs)

        # Stash number of lons/lats in named attributes
        self.attrs['nLat'], self.attrs['nLon'] = self['grid']

    def add_slice(self, value, target=None, loc=111, xlim=[-15,15], ylim=[-15,15], 
                  zlim=None, **kwargs):
        '''
        Add a 2D plot of *value* in the equatorial plane.
        '''
        
        from spacepy.pybats import add_planet

        # Start by prepping values for plotting:
        x, y = np.array(self['x']), np.array(self['y'])
        z = np.array(self[value])

        #TEMP: 
        zlim=[z.min(), z.max()]

        #add points at "infinity" with values set zmin so that color
        #background fills in smoothly
        x = np.append(x,[xlim[0], xlim[0], xlim[1], xlim[1]])
        y = np.append(y,[ylim[0], ylim[1], ylim[0], ylim[1]])
        z = np.append(z, 4*[zlim[0]])

        # Generate target figure and axes if not given:
        fig, ax = set_target(target, figsize=(10, 10), loc=loc)
        
        cont = ax.tricontourf(x, y, z, **kwargs)
        
        cbar=None

        # Configure axes:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        add_planet(ax)

        return fig, ax, cont, cbar

    def extract(self):
        pass
