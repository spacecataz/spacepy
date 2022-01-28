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
import matplotlib.pyplot as plt

from spacepy.plot import set_target
from spacepy.pybats import IdlFile


def parse_cimihead(headstring):
    '''Given a CIMI-style header line that reads `Var1[Units1] Var2[Units2] ...`
       return a list of variables and a list of units.'''

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

    def add_contour(self, value, target=None, loc=111, xlim=[-15, 15],
                    ylim=[-15, 15], zlim=None, title=None, nlev=51,
                    add_cbar=False, clabel=None, dolog=False,
                    show_points=False, add_body=True,
                    *args, **kwargs):
        '''
        Create a contour plot of variable **value** in the equatorial plane.

        Simple example:

        >>> self.add_contour('Pot')

        If kwarg `target` is `None` (default), a new figure is
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object,
        the plot is placed into that axis.

        Four values are returned: the matplotlib Figure and Axes objects,
        the matplotlib contour object, and the matplotlib colorbar object
        (defaults to *False* if not used.)

        Extra args and kwargs are passed to matplotlib's contour function.

        Parameters
        ----------
        value : str
            The name of the variable to plot.

        Other Parameters
        ----------------
        target : matplotlib figure or axes, optional
            Figure or axes object on which to place resulting plot. Default
            behavior is to create a new figure and axes.
        loc : str or int, default=111
            Sets the subplot position of the axes if a new axes is created.
        xlim, ylim : array-like, default=[-15,15]
            Sets the range of the x and y-axes.
        zlim : array-like, optional
            Sets the contour range, e.g., [0,100]. Default behavior is to
            use the plotted variable's range.
        nlev : int, default=51
            Sets the number of contour levels.
        xlabel, ylabel : str, optional
            Labels for the x and y axes, respectively.  Defaults to 'X $R_E$'.
        title : str, optional
            Plot title; defaults to `value`.
        add_cbar : bool, default=False
           Toggle the color bar on or off.
        clabel : str, optional
           Sets colorbar label.  Defaults to `value` and units.
        add_body : bool, default=False
            Places planetary body in plot.
        dolog : bool, default=False
            Sets use of logarithmic scale for contour.
        show_points: bool, default=False
            Toggles plotting of CIMI cell center locations.

        Returns
        -------
        figure object
            The figure onto which the plot is placed.
        axes object
            The figure onto which the plot is placed.
        contour object
            The matplotlib contour object.
        colorbar object or None
            If a colorbar is created, it is returned.
        '''

        from matplotlib.colors import LogNorm
        from matplotlib.ticker import LogLocator, LogFormatterMathtext
        from spacepy.pybats import add_planet

        # Start by prepping values for plotting:
        x, y = np.array(self['x']), np.array(self['y'])
        z = np.array(self[value])

        # Set the contour range intelligently:
        if zlim is None:
            zlim = [z.min(), z.max()]
            if dolog and zlim[0] <= 0:
                zlim[0] = np.min([0.0001, zlim[1]/1000.0])

        # Create levels and set norm based on dolog.
        if dolog:
            # Create log-spaced contours; limit range of plotted data:
            levs = np.power(10, np.linspace(np.log10(zlim[0]),
                                            np.log10(zlim[1]), nlev))
            z = np.where(z > zlim[0], z, 1.01*zlim[0])
            norm, ticks, fmt = LogNorm(), LogLocator(), LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            norm, ticks, fmt = None, None, None

        # Add points at "infinity" with values set zmin so that color
        # background fills in smoothly
        x = np.append(x, [xlim[0], xlim[0], xlim[1], xlim[1]])
        y = np.append(y, [ylim[0], ylim[1], ylim[0], ylim[1]])
        z = np.append(z, 4*[zlim[0]])

        # Generate target figure and axes if not given:
        fig, ax = set_target(target, figsize=(10, 10), loc=loc)

        # Add contour to axes:
        cont = ax.tricontourf(x, y, z, levs, *args, norm=norm, **kwargs)

        # Add cbar if requested:
        cbar = None
        if add_cbar:
            cbar = plt.colorbar(cont, ax=ax, ticks=ticks, format=fmt, pad=0.01)
            if clabel is None:
                clabel = "%s (%s)" % (value, self[value].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar = None  # Need to return something, even if none.

        # Configure axes:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        if add_body:
            add_planet(ax)

        # Set labels:
        if title is not None:
            ax.set_title(title)
        else:
            ax.set_title(value)

        return fig, ax, cont, cbar

    def extract(self):
        pass
