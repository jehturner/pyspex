#!/bin/env python

import numpy as np
import scipy.ndimage.interpolation as ndi
from astropy.modeling.models import Mapping, Identity, Chebyshev1D, Chebyshev2D
from astropy.modeling.fitting import LinearLSQFitter
from astropy import units as u
from ndmapper.data import *
from gwcs import wcs, coordinate_frames as cf
from qddb import read_iddb
from wcs import wcs_append
from fitting import fit_with_inverse

# import matplotlib.pyplot as plt


extracted_frame = cf.Frame2D(name="extracted", axes_names=("x", "y"),
                             unit=(u.pix, u.pix))
common_frame = cf.Frame2D(name="common_wave", axes_names=("refx", "y"),
                          unit=(u.pix, u.pix))

dist_fit = LinearLSQFitter()
wave_fit = LinearLSQFitter()  # re-use the same instance??
#or_fit = fitting.FittingWithOutlierRemoval(fit, sigma_clip,
#                                           niter=3, sigma=3.0)

df = DataFile('eprgS20120828S0005.fits')

# Determine the distortions for each extension and add a WCS object for them:
refpks = []
for n, ndd in enumerate(df):

    # Determine the transformation between each row of the 2D spectrum and a
    # reference row. For slit spectra, this model can be represented as a 2D
    # transform, whereas for IFU spectra it might want to be a collection of
    # 1D mappings (check Nadia's example when trying to implement that)?
    dist_model = Chebyshev2D(4, 4)

    # Want to fit the (currently 0) Chebyshev model to samples of the mapping
    # between line positions in this vector and the reference row.

    # PLACEHOLDER: Get line centre measurements from the IRAF database, just
    # to help implement a proof-of-concept sketch and then add real peak
    # detection and matching later:
    record_name = df.filename.root + '_{:03d}'.format(n+1)
    features = read_iddb(record_name)
    # feature_list = [features[key] for key in sorted(features.iterkeys())]
    # Convert to lists of peak centre and reference wavelength tuples:
    feature_list = [sorted([(float(mx), float(rwl)) for rwl, (mx, fwl) \
                            in features[key].iteritems()]) \
                    for key in sorted(features.iterkeys())]

    nspec = ndd.shape[0]
    if len(feature_list) != nspec:
        raise ValueError('Image / database aperture mismatch')

    # Use the middle aperture as the reference for now:
    refspec = nspec / 2
    x, y, x1 = [], [], []

    refpks.append(feature_list[refspec])
    nrefpks = len(refpks[-1])

    # Determine x_ref as a function of x_row wherever there is a close match
    # in their nominal wavelengths.
    for rownum in range(nspec):
        rowpks = feature_list[rownum]
        for npk, (refpkx, refpkwl) in enumerate(refpks[-1]):
            found = False
            for peakx, peakwl in rowpks:
                # For some reason the nominal line wavelengths differ slightly?
                if abs(peakwl - refpkwl) < 0.1:
                    found = True
                    break
            if found:
                x.append(peakx)
                y.append(rownum)
                x1.append(refpkx)

    x = np.array(x)
    y = np.array(y)
    x1 = np.array(x1)

    # Fit pixel position in the reference row or shift WRT it as a function of
    # pixel index in each row.
    spec_dist, ispec_dist = fit_with_inverse(dist_fit, dist_model, x, y, x1)

    # tx = np.arange(0, 3000)
    # ty = np.zeros_like(tx) + 100
    # tx1 = spec_dist(tx, ty)
    # tx2 = ispec_dist(tx1, ty)
    # tx3 = spec_dist(tx2, ty)
    # plt.plot(tx, tx1)
    # plt.plot(tx, tx3)
    # plt.show()

    #plt.plot(x[100])

    # Combine into a 2D transformation that just preserves the y co-ordinate:
    spec_dist_mapping = Mapping((0, 1, 1)) | spec_dist & Identity(1)
    spec_dist_mapping.inverse = Mapping((0, 1, 1)) | ispec_dist & Identity(1)

    # Ref. frames & mappings between them for constructing a gwcs object:
    transforms = [spec_dist_mapping, None]
    frames = [extracted_frame, common_frame]

    # Test: evaluate & write out mapping.
    # y, x = np.mgrid[0:ndd.data.shape[0], 0:ndd.data.shape[1]]
    # ndd.data, y = spec_dist_mapping(x, y)

    # The current API only allows setting the WCS at the time of creation,
    # raising a problem when one wants to load a dataset and determine and
    # attach its WCS, requiring a new instance:
    # df[n] = NDLater(data=ndd, wcs=wcs.WCS(forward_transform=spec_dist_mapping,
    #                                       input_frame=extracted_frame,
    #                                       output_frame=common_frame))
    df[n] = NDLater(data=ndd, wcs=wcs.WCS([(f, t) for f, t \
                                          in zip(frames, transforms)]))


rss_frame = cf.CoordinateFrame(1, 'SPATIAL', (1,), unit=u.pix, axes_names='y',
                               name='row')
spec_frame = cf.SpectralFrame((0,), unit=u.angstrom, axes_names='lambda',
                               name='wavelength')
wavelength_rss_frame = cf.CompositeFrame([spec_frame, rss_frame])

# wave_frame = cf.Frame2D(name="wavelength", axes_names=("refx", "y"),
#                           unit=(u.pix, u.pix))

# Determine the line list match for the reference row. We could also add
# features detected in other rows, using the distortion mapping to convert the
# centres to equivalent co-ordinates.
for n, ndd in enumerate(df):

    wave_model_ref = Chebyshev1D(4)

    # Dummy match to the line list (note that this depends on the above loop
    # to get the IRAF database info. in the current hack):
    x = np.array([xcen for xcen, wcen in refpks[n]])
    l = np.array([wcen for xcen, wcen in refpks[n]])
    print x[3], l[3]

    wave_dist, iwave_dist = fit_with_inverse(wave_fit, wave_model_ref, x, l)

    wave_mapping = wave_dist & Identity(1)
    wave_mapping.inverse = iwave_dist & Identity(1)

    # Append new stage to a copy of the existing WCS:
    #print type(ndd.wcs.get_transform("extracted", "common_wave"))
    df[n] = NDLater(data=ndd, wcs=wcs_append(ndd.wcs, wave_mapping,
                                             wavelength_rss_frame))

    # y, x = np.mgrid[0:ndd.data.shape[0], 0:ndd.data.shape[1]]
    # ndd.data[:], y = df[n].wcs(x, y)
    # w1 = df[n].wcs(1131.41, 370)
    # w2 = df[n].wcs(1114.49, 100)
    # w3 = df[n].wcs(1103.75, 745)
    # w1i = df[n].wcs.invert(*w1)
    # w2i = df[n].wcs.invert(*w2)
    # w3i = df[n].wcs.invert(*w3)
    # print 'w', w1, w2, w3
    # print 'wi', w1i, w2i, w3i

df.filename="test.fits"
#df.save()

# 370 6416.25 = 6416.3 at 1131.4
# 100 at 1114.49 (gives 6425)

