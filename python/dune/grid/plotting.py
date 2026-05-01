"""DUNE Grid plotting utility.

This module provides a collection of functions for plotting.

Example:

    Simple plotting::

        from dune.grid.plotting import plot
        plot(uh)


Example:

    Multiple and customized visualizations in one figure::

        import matplotlib.pyplot as plt
        from dune.grid.plotting import plotData, plotGrid, plotContours, plotQuivers, plotColorbar

        fig = plt.figure()
        plotData(uh, fig, cmap='jet')
        plotGrid(grid, fig, color='white')
        plotContours(uh, fig, levels=10)
        plotQuivers(uh, fig, vectors=[0,1])
        plotColorbar(uh, fig)
        plt.show()


Todo:
    * Clarify function signatures
      - Is it fine to require a figure? We could alternatively create a new one and return it, but we wouldn't call plt.show().
      - Most of the parameters can be passed to pyplot by **kwargs. We have to figure what to pass where exactly.
      - level: Do we want/need level? Everywhere?
      - xlim, ylim: Can be controlled from outside via ax.margins(...). We should provide an example.
      - logscale: Do we want to provide this? If yes, everywhere?
      - Can't we simply use plot[Point]Data for plotCellData by converting the data?
      - What is 'mayaviPointData'?

    * 'plot' shortcut:
      - I think we could replace plotComponents by adding a loop over all components here.
      - Should this also work for grids? plotGrid requires a figure so far.

    * Should we also forward the plot methods on the objects? E.g. uh.plot() and grid.plot(fig).
    * Add more examples for various scenarios. Add one using the triangulation and (point-)data directly.

    * Deprecate and remove implementations in dune.common.plotting and dune.fem.plotting
"""

def plot(solution):
    fig = plt.figure()
    plotData(solution, fig)
    plotColorbar(solution, fig)
    plt.show()


def plotGrid(grid, figure, **kwargs):
    pass


def plotData(solution, figure, **kwargs):
    pass


def plotContours(solution, figure, **kwargs):
    pass


def plotQuiver(solution, figure, nofVectors=None, **kwargs):
    assert len(solution) == 2
    pass


def plotColorbar(solution, figure, ticks=11, orientation="vertical", **kwargs):
    pass
