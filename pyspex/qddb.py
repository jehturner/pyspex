# Quick and dirty parse of IRAF database, since the gemini_python version
# wasn't working and will probably take longer to fix.

import re

def read_iddb(name):
    """
    Returns { aperture : { feature : (xcen, wlcen) } }.
    """

    with open('database/id' + name) as dbfile:

        lines = dbfile.readlines()

    apnum=-1
    features_block=False
    features = dict()

    for line in lines:

        words = line.split()

        if words:

            if words[0] == 'aperture':
                apnum = int(words[1])

            elif words[0] == 'features':
                features_block = True
                context_indent = line.replace('\t', ' '*8).find('features')
                features[apnum] = dict()
                # print (apnum)

            elif features_block:
                line = line.replace('\t', ' '*8)
                if len(line) - len(line.lstrip(' ')) > context_indent:
                    x, wl, refwl = words[0:3]
                    features[apnum][refwl] = x, wl
                else:
                    features_block = False

    return features

