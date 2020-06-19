import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import astropy.io.fits as fits

def main(plate=4458,mjd=55536,fiber=596):
    if plate >=10000:
        filename = 'spec-{0:05d}-{1:05d}-{2:04d}.fits'.format(plate,mjd,fiber)
    else:
        filename = 'spec-{0:04d}-{1:05d}-{2:04d}.fits'.format(plate,mjd,fiber)
    
    data = fits.open('BALQSOs_Spectra/'+filename)[1].data

    x, y = 10**data.loglam,data.flux
    fig, ax = plt.subplots()
    ax.plot(x, y,'.', color='black')
    highlighter = Highlighter(ax, x, y)
    plt.show()

    selected_regions = highlighter.mask
    # Print the points _not_ selected
    print x[~selected_regions], y[~selected_regions]


class Highlighter(object):
    def __init__(self, ax, x, y):
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.x, self.y = x, y
        self.mask = np.zeros(x.shape, dtype=bool)

        self._highlight = ax.scatter([], [], s=1, color='yellow', zorder=10)

        self.selector = RectangleSelector(ax, self, useblit=True)

    def __call__(self, event1, event2):
        self.mask |= self.inside(event1, event2)
        xy = np.column_stack([self.x[self.mask], self.y[self.mask]])
        self._highlight.set_offsets(xy)
        self.canvas.draw()

    def inside(self, event1, event2):
        """Returns a boolean mask of the points inside the rectangle defined by
        event1 and event2."""
        # Note: Could use points_inside_poly, as well
        x0, x1 = sorted([event1.xdata, event2.xdata])
        y0, y1 = sorted([event1.ydata, event2.ydata])
        mask = ((self.x > x0) & (self.x < x1) &
                (self.y > y0) & (self.y < y1))
        return mask

main()
