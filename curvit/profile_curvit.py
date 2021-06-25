#!/usr/bin/env python3

import cProfile
import pstats

from curvit import makecurves


cProfile.run('makecurves(events_list = "AS1G06_084T01_9000000710uvtFIIPC00F1_l2ce.fits.gz",\
                         background = "auto")', 
             'profile.stats')


p = pstats.Stats('profile.stats')
p.strip_dirs().sort_stats('time').print_stats('curvit', 15)
