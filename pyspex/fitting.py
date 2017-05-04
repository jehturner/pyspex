from copy import deepcopy
import numpy as np


def fit_with_inverse(fitter, model, x, y, x1):

    # Fit x', y = f(x, y) and also its inverse, x, y = f'(x', y). This relies
    # on attributes defined in 2D model implementations.

    ifitter = deepcopy(fitter)  # not sure whether this is necessary?

    # Determine the forward fit as usual:
    fwd = fitter(model, x, y, x1)

    # Construct a grid on which to fit the inverse model. Extend the original
    # fitting domain by 2x, to obtain an inverse that's valid for slightly
    # extrapolated points at the edges of the data. Increase the order by one
    # term, for accurate approximation:
    halfx = 0.5 * (fwd.x_domain[1] - fwd.x_domain[0])
    xlow, xhigh = fwd.x_domain[0] - halfx, fwd.x_domain[1] + halfx
    xdeg = fwd.x_degree + 1
    ylow, yhigh = fwd.y_domain[0], fwd.y_domain[1]
    ydeg = fwd.y_degree
    density = 2  # sample points per unknown coefficient, min 1

    # Construct a sparsely-sampled grid (but not a sparse one in the sense of
    # being an implict vector product, since one of the axes needs inverting):
    ix, iy = np.meshgrid(np.linspace(xlow, xhigh, density*xdeg),
                         np.linspace(ylow, yhigh, density*ydeg), sparse=False)

    # Evaluate the forward transform on the sub-sampled grid:
    ix1 = fwd(ix, iy)

    # Fit an inverse model to the points from the forward one:
    inv = ifitter(type(fwd)(xdeg, ydeg), ix1, iy, ix)

    # # Quick test
    # f0, f100, f2000 = fwd(0, 0), fwd(100, 100), fwd(2000, 300)
    # i0, i100, i2000 = inv(f0, 0), inv(f100, 100), inv(f2000, 300)
    # ff0, ff100, ff2000 = fwd(i0, 0), fwd(i100, 100), fwd(i2000, 300)
    # print 'F0', f0, i0, ff0
    # print 'F100', f100, i100, ff100
    # print 'F2000', f2000, i2000, ff2000

    return fwd, inv


