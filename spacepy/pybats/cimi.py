'''
CIMI - A module for interfacing with output from the Comprehensive 
Inner Magnetosphere Ionosphere Model.

.. currentmodule:: spacepy.pybats.cimi

.. rubric:: Classes

.. autosummary::
    :template: clean_class.rst
    :toctree:

    Cimi2d

'''

from spacepy.pybats import IdlFile

class Cimi2d(IdlFile):
    '''
    A class to handle 2D CIMI output files.
    '''

    # Init by calling IdlFile init and then building qotree, etc.
    def __init__(self, filename, *args, **kwargs):
        
        
        # Read file.
        IdlFile.__init__(self, filename, header=None, gridtype="Irregular",
                         *args, **kwargs)
