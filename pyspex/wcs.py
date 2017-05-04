from gwcs import wcs

def wcs_append(wcs_obj, forward_transform, output_frame):
    pipeline = [(getattr(wcs_obj, f) or f, m) for f, m in wcs_obj.pipeline]
    pipeline[-1] = (pipeline[-1][0], forward_transform)
    pipeline.append((output_frame, None))
    return wcs.WCS(pipeline)


