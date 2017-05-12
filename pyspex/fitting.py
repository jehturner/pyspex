from copy import deepcopy
import numpy as np


def fit_with_inverse(fitter, model, x, y, z=None):
    """
    Fit x' = f(x) and also its inverse, x = f'(x')
     or x' = f(x, y)     and also       x = f'(x', y).

    Ideally, this would be more general, taking an arbitrary 2D model and
    returning a 2D model with its inverse already attached, but there are
    currently some impediments to doing that:

      https://github.com/astropy/astropy/issues/6012
      https://github.com/astropy/astropy/issues/6026
      https://github.com/astropy/astropy/issues/6038

    Returns
    -------

    (fwd : `~astropy.modeling.FittableModel`,
     inv : `~astropy.modeling.FittableModel`)

    """

    # This relies on attributes defined in 1D/2D model implementations.

    ifitter = deepcopy(fitter)  # not sure whether this is necessary?

    # Determine the forward fit as usual:
    fwd = fitter(model, x, y, z)

    # Determine a sparsely-sampled grid on which to fit the inverse model.
    # Extend the original fitting domain by 2x, to obtain an inverse that's
    # valid for slightly extrapolated points at the edges of the data and
    # increase the order by one term, for accurate approximation.
    xdom = fwd.domain if z is None else fwd.x_domain
    xdeg = (fwd.degree if z is None else fwd.x_degree) + 1
    halfx = 0.5 * (xdom[1] - xdom[0])
    xlow, xhigh = xdom[0] - halfx, xdom[1] + halfx

    density = 2  # sample points per unknown coefficient, min 1

    # Construct the sub-sampled grid and evaluate the forward transform on it:
    if z is None:
        coords = [np.linspace(xlow, xhigh, density*xdeg)]
        imodel = type(model)(xdeg)
    else:
        ylow, yhigh, ydeg = fwd.y_domain[0], fwd.y_domain[1], fwd.y_degree
        coords = np.meshgrid(np.linspace(xlow, xhigh, density*xdeg),
                             np.linspace(ylow, yhigh, density*ydeg),
                             sparse=False)
        imodel = type(model)(xdeg, ydeg)

    # Fit an inverse model to the points from the forward one:
    inv = ifitter(imodel, fwd(*coords), *reversed(coords))

    # Return the fits separately for now, since one can't be the inverse of
    # the other until they have the same number of inputs & outputs:
    return fwd, inv

