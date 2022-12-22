#!/usr/bin/env python3


'''Create light curves from UVIT data. 


   Copyright 2019 Prajwel Joseph
  
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
  
       http://www.apache.org/licenses/LICENSE-2.0
  
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.'''


import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import os
import sys
import ntpath
import random
import numpy as np
import astroalign as aa
import matplotlib.pyplot as plt

from glob import glob
from collections import Counter
from scipy.ndimage import zoom
from scipy.spatial import KDTree
from scipy.interpolate import interp1d
from matplotlib.colors import LogNorm
from photutils import DAOStarFinder, CircularAperture
from photutils.background import Background2D, MedianBackground
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma, sigma_clip, SigmaClip
from astroquery.astrometry_net import AstrometryNet
 
#######################################################################
# Initial set of parameters.

'''Window size Vs Framecount rate dictionary (approximate). 
The most accurate way to get the rate would be to take the value 
of (1 / INT_TIME). INT_TIME value can be found from the image header
window_rate_dict = {'512 x 512': 28.7185,
                    '350 x 350': 61.0,
                    '300 x 300': 82.0,
                    '250 x 250': 115.0,
                    '200 x 200': 180.0,
                    '150 x 150': 300.0,
                    '100 x 100' : 640.0}
'''

events_list = '' #events file
radius = 6  # radius of aperture in pixels.
sky_radius = 12 # radius of background aperture in pixels.
bwidth = 50 # bin width in seconds, change as you please. 
framecount_per_sec = 28.7185  # 28.7185 frames / second for 512x512 mode.

# The following is only required for makecurves. 
detection_method = 'daofind' # valid inputs are 'daofind' / 'kdtree'.
threshold = 4 # threshold ('daofind').
how_many = 4 # limit ('kdtree').

# The coordinates are only required for curves.
xp = None
yp = None

'''The following parameter affects how the background estimation is done.
The default value is 'None' and no background estimation is carried out.
If you prefer to manually specify a background region, then give 'manual' 
as the value. Also, you can provide 'auto' and background region 
will be automatically selected.'''
background = None  # valid inputs are None / 'manual' / 'auto'.

# If 'manual', PLEASE FILL the following.
x_bg = None # background X-coordinate.
y_bg = None # background Y-coordinate.


'''The following parameters determines whether corrections are 
applied to the CPF. They are aperture-correction and
saturation-correction.'''
aperture_correction = None # valid inputs are None / 'fuv' / 'nuv'.
saturation_correction = False # True or False.


# Following parameters need not be changed (unless you want to).
whole_figure_resolution = 256 # resolution of full figure.
sub_fig_size = 40 # size of sub figure.
fontsize = 9 #fontsize for plots.

# Encircled energy data (https://doi.org/10.3847/1538-3881/ab72a3).
radius_pixels = np.array([ 1.5, 2, 2.5, 3, 4,
                           5, 7, 9, 12, 15, 
                           20, 30, 40, 50, 70, 
                           80, 95 ])

nuv_energy_percentage = np.array([ 29.9, 42.0, 52.0, 59.3, 68.8,  
                                   74.5, 81.3, 85.1, 89.3, 92.1,  
                                   95.2, 97.6, 98.4, 98.8, 99.4,  
                                   99.6, 100.0 ])

fuv_energy_percentage = np.array([ 28.1, 40.7, 51.1, 59.1, 68.9,
                                   74.6, 81.4, 85.0, 88.6, 91.3,  
                                   94.5, 96.9, 97.7, 98.3, 99.1,  
                                   99.5, 100.0 ])
                                   
#ratio = measured CPF / actual CPF
nuv_ratio = nuv_energy_percentage / 100.
fuv_ratio = fuv_energy_percentage / 100.
fuv_ratio_function = interp1d(radius_pixels, fuv_ratio, kind = 'cubic')
nuv_ratio_function = interp1d(radius_pixels, nuv_ratio, kind = 'cubic')

# The parameters below are for combining events lists.

'''single_star option is useful only when there is only one star in
the field and no rotation between frames.'''
shift_algorithm = 'multiple_star' # 'single_star' or 'multiple_star'.

min_exptime = 100 # will ignore orbits with exptimes below limit.

'''Astroalign settings. Change only if you know what you are 
doing.'''
aa.NUM_NEAREST_NEIGHBORS = 7
aa.MIN_MATCHES_FRACTION = 0.1
aa.PIXEL_TOL = 2

# For Astrometry.
AstrometryNet_API_key = 'ujmrvwqqyelxmzcj'  
#######################################################################


def read_columns(events_list):
    # Reading few columns.
    f = fits.open(events_list)
    time = f[1].data['MJD_L2']
    fx = f[1].data['Fx']
    fy = f[1].data['Fy']
    photons = f[1].data['EFFECTIVE_NUM_PHOTONS']
    bad_flag = f[1].data['BAD FLAG']
    mask = photons > 0 
    mask = np.logical_and(mask, bad_flag)
    time = time[mask]
    fx = fx[mask]
    fy = fy[mask]
    photons = photons[mask]
    return time, fx, fy, photons

def tobe_or_notobe(time, bwidth, 
                   detection_method,
                   threshold,
                   how_many, 
                   background, 
                   x_bg, y_bg,
                   aperture_correction, radius,
                   saturation_correction):
    
    sanity = (time.max() - time.min()) / bwidth
    if int(sanity) < 1:
        print('\nEvents list contain little data OR check bwidth parameter.\n')

    if detection_method not in ['daofind', 'kdtree']:
        print('\nInvalid input for "detection_method" parameter.\n')
        sanity = 0 
        
    if threshold == 0:
        print('\nThe "threshold" parameter is set at 0.\n')
        sanity = 0 
        
    if how_many == 0:
        print('\nThe "how_many" parameter is set at 0.\n')
        sanity = 0 
    
    if background not in [None, 'auto', 'manual']:
        print('\nInvalid input for "background" parameter.\n')
        sanity = 0 
        
    if background == 'manual':
        if None in [x_bg, y_bg]:
            print('\nPlease provide values for both "x_bg" and "y_bg".\n')
            sanity = 0

    if aperture_correction not in [None, 'fuv', 'nuv']:
        print('\nInvalid input for "aperture_correction" parameter.\n')
        sanity = 0 
        
    if saturation_correction not in [True, False]:
        print('\nInvalid input for "saturation_correction" parameter.\n')
        sanity = 0
    
    if aperture_correction != None:
        if 1.5 <= radius <= 95:
            pass
        else:
            print('\nThe "radius" parameter should be in the range {1.5, 95}')
            sanity = 0 
    return int(sanity) 

def modify_string(events_list):
    if events_list[-5:] == '.fits':
        events_list = events_list[:-5]
    if events_list[-8:] == '.fits.gz':
        events_list = events_list[:-8]
    return events_list

# To automatically choose background region.
def auto_bg(fx, fy, time, photons, radius, framecount_per_sec, sky_radius): 
    weights = photons / framecount_per_sec    
    bins = np.arange(0, 4801, 16)    
    lowres_counts, lowres_xedges, lowres_yedges = np.histogram2d(fx, fy, 
                                                                 bins = (bins, bins), 
                                                                 weights = weights)  
    lowres_xcentres = (lowres_xedges[:-1] + lowres_xedges[1:]) / 2.
    lowres_ycentres = (lowres_yedges[:-1] + lowres_yedges[1:]) / 2.
    flat_counts = lowres_counts.flatten()
    x_mesh, y_mesh = np.meshgrid(lowres_ycentres, lowres_xcentres) #Notice the swap.
    x_mesh = x_mesh.flatten() 
    y_mesh = y_mesh.flatten()
    # To avoid the edges. 
    mask_radius = 1700
    mask = (x_mesh - 2400) ** 2 + (y_mesh - 2400) ** 2 <= mask_radius ** 2
    array = np.array([flat_counts, x_mesh, y_mesh]).T 
    polished_array = array[mask]
    bg_mask = sigma_clip(polished_array[:, 0], sigma = 3, maxiters = 5)
    bg_mask = np.logical_not(bg_mask.mask)
    
    bg_CPS_sample = []
    bg_CPS_e_sample = []
    sample_size = 100
    for d in range(sample_size):
        r_count, x_bg, y_bg = random.choice(polished_array[bg_mask])
        bg_CPS, bg_CPS_e = bg_estimate(fx, fy, time, photons, framecount_per_sec,
                                       radius, x_bg, y_bg, sky_radius)
        bg_CPS_sample.append(bg_CPS)
        bg_CPS_e_sample.append(bg_CPS_e)
        
    bg_CPS_sample = np.array(bg_CPS_sample)
    bg_CPS_e_sample = np.array(bg_CPS_e_sample)
    bg_CPS_mask = sigma_clip(bg_CPS_sample, sigma = 3, maxiters = 5)
    bg_CPS_mask = np.logical_not(bg_CPS_mask.mask)
    bg_CPS = np.mean(bg_CPS_sample[bg_CPS_mask])
    bg_CPS_e = np.mean(bg_CPS_e_sample[bg_CPS_mask])
    return lowres_counts, bg_CPS, bg_CPS_e   
    
# To estimate background CPS.
def bg_estimate(fx, fy, time, photons, framecount_per_sec, radius, x_bg, y_bg, sky_radius):

    weights = photons / framecount_per_sec    
    mask = ((fx - x_bg) ** 2 + (fy - y_bg) ** 2) <= sky_radius ** 2    
    T = time[mask]
    W = weights[mask]

    if len(T) != 0:
        scaled_events = (np.sum(W) * radius ** 2) / float(sky_radius ** 2)
        scaled_events_e = (np.sqrt(len(T)) * radius ** 2) / float(sky_radius ** 2)
    else:
        scaled_events = 0
        scaled_events_e = 0

    unique_time = np.unique(time)
    Number_of_frames = float(len(unique_time))
    bg_CPS = (scaled_events * framecount_per_sec) / Number_of_frames
    bg_CPS_e = (scaled_events_e * framecount_per_sec) / Number_of_frames
    return bg_CPS, bg_CPS_e

# To create subset images. 
def create_sub_image(pos_x, pos_y, 
                     sub_size, 
                     cir_rad, 
                     sub_name,
                     fx, fy,
                     path_to_events_list,
                     events_list):

    mask = np.logical_and(np.abs(fx - pos_x) <= sub_size,
                          np.abs(fy - pos_y) <= sub_size)

    obj_fx = fx[mask]
    obj_fy = fy[mask]
    
    obj_circle = plt.Circle((pos_x, pos_y), cir_rad, 
                            color = 'k', fill = False)

    plt.hist2d(obj_fx, obj_fy, bins = sub_size * 2, norm = LogNorm())
    plt.gcf().gca().add_artist(obj_circle)
    source_png_name = os.path.join(path_to_events_list, sub_name + events_list + '.png')
    plt.savefig(source_png_name, format = 'png', bbox_inches = 'tight')
    plt.clf()
    return source_png_name  

# To find positions of interest (using daofind algorithm).
def detect_sources_daofind(fx, fy, photons, threshold):
    mask_radius = 1700
    kernel = Gaussian2DKernel(x_stddev = 4 * gaussian_fwhm_to_sigma)

    aperture = CircularAperture((2400, 2400), r = mask_radius)
    mask = aperture.to_mask(method = 'center')
    mask = mask.to_image(shape = ((4800, 4800)))

    weights = photons / framecount_per_sec
    bins = np.arange(0, 4801)
    ndarray, yedges, xedges = np.histogram2d(fy, fx, bins = (bins, bins), weights = weights)  
    data = ndarray * mask

    mean, _, std = sigma_clipped_stats(data, sigma = 5., maxiters = 1)
    # to avoid std becoming zero. 
    if std == 0:
        mean = np.mean(data)
        std = np.std(data)
        
    data = convolve(data, kernel)
    daofind = DAOStarFinder(fwhm = 3.0, threshold = threshold * std, exclude_border = True)
    sources = daofind(data - mean)    
    sources.sort('mag')

    uA = np.array([sources['xcentroid'].data, sources['ycentroid'].data]).T
    uA = np.round(uA, 2)
    return uA

# To find positions of interest (positions with maximum events).
def detect_sources_kdtree(fx, fy, how_many):
    fxi = [int(round(s)) for s in fx]
    fyi = [int(round(s)) for s in fy]
    # Counting stuff to know who all are popular. 
    counter = Counter(zip(fxi, fyi))
    most_common_limit = int(100 * how_many)
    A = np.array(list(zip(*counter.most_common(most_common_limit)))[0])
    # Sieving out the duplicates. 
    uA = []
    while len(A) != 0:
        uA.append(A[0]) 
        A = np.array([x for x in A 
                      if x not in A[KDTree(A).query_ball_point(uA[-1], 16)]])

    uA = np.array(uA)[:how_many]
    uA = uA

    if len(uA) != 0:
        # To avoid sources which have coordinates 0 or 4800.
        mask = np.isin(uA, [0, 4800], invert = True)
        mask = mask[:, 0] * mask[:, 1]
        uA = uA[mask]
    return uA

def get_counts(fx, fy, time, photons, framecount_per_sec, xp, yp, radius):

    weights = photons / framecount_per_sec    
    # selecting events within a circular region.
    mask = ((fx - xp) ** 2 + (fy - yp) ** 2) <= radius ** 2    
    T = time[mask]
    W = weights[mask]

    # To find Counts per Frame (CPF).
    unique_time = np.unique(time)
    Number_of_frames = float(len(unique_time))
    CPF = np.sum(W) / Number_of_frames
    CPF_err = np.sqrt(len(T)) / Number_of_frames
    return CPF, CPF_err

# To change mission elapsed time in seconds to modified julian date.
def met_to_mjd(met):
    jan2010 = 55197.0  # 2010.0(UTC) expressed with MJD format and scale UTC.
    mjd = (met / 86400.0) + jan2010  # 1 julian day = 86400 seconds.
    return mjd

def apply_aperture_correction(CPF, CPF_err, radius, aperture_correction): 
    if aperture_correction == 'fuv':
        CPF = CPF / fuv_ratio_function(radius)
        CPF_err = CPF_err / fuv_ratio_function(radius)
    elif aperture_correction == 'nuv':
        CPF = CPF / nuv_ratio_function(radius)
        CPF_err = CPF_err / nuv_ratio_function(radius)
    return CPF, CPF_err
    
def apply_saturation_correction(CPF5, CPF5_err, saturation_correction):
    
    if saturation_correction == True:
        
        if np.sum(CPF5 >= 0.6) != 0:
            print("\nCounts per frame exeeds 0.6; saturation correction cannot be applied")
            return
        
        ICPF5 = -1 * np.log(1 - CPF5)
        ICPF5_err = CPF5_err / CPF5
        
        ICORR = ICPF5 - CPF5
        ICORR_err = np.sqrt((ICPF5_err ** 2) + (CPF5_err ** 2))
        
        RCORR = ICORR * (0.89 - (0.30 * (ICORR ** 2)))
        RCORR_err = RCORR * np.sqrt((ICORR_err ** 2) + ((0.30 * 2 * ICORR * ICORR_err) ** 2))
        
        CPF5 = CPF5 + RCORR
        CPF5_err = np.sqrt((CPF5_err ** 2) + (RCORR_err ** 2))       
    return CPF5, CPF5_err  
    
def makecurves(events_list = events_list,
               radius = radius,
               detection_method = detection_method,
               threshold = threshold,
               how_many = how_many,
               bwidth = bwidth,
               framecount_per_sec = framecount_per_sec,
               background = background,
               sky_radius = sky_radius,
               x_bg = x_bg,
               y_bg = y_bg,
               aperture_correction = aperture_correction,
               saturation_correction = saturation_correction,
               whole_figure_resolution = whole_figure_resolution,
               sub_fig_size = sub_fig_size, 
               fontsize = fontsize):

    """Automatically detect sources amd create light curves.

    Parameters
    ----------
    events_list : file path
        The name of the events list FITS file.
                       
    radius : float, optional
        The source aperture radius in pixels. 
        This parameter has a default value of 6.
        
    detection_method : {'daofind', 'kdtree'}, optional
        The parameter to choose between available detection methods. 
        
        * ``'daofind'``: To use the DAOFIND algorithm. This is the default method.
        * ``'kdtree'``: Source detection method based on a k-d tree implementation. 
        
    threshold : float, optional
        The threshold parameter associated with the ``'daofind'`` method. 
        The default value is 4.
        
    how_many : int, optional
        The limit for the number of sources to be detected using 
        the ``'kdtree'`` method. 
        The default value is 4. 
        
    bwidth : float, optional
        Time bin width in seconds. 
        the default value is 50.   
        
    framecount_per_sec : float, optional
        The framerate of the observation, with a default value of 28.7185
        frames per second for 512 x 512 window mode. 
        The most accurate way to get the framerate would be to take the value 
        of (``1 / INT_TIME``). 
        ``INT_TIME`` value can be found from the corresponding image header. 
        Approximate values of framerate for different window modes of UVIT 
        are given in the table below.

        +---------------+---------------------+
        | window mode   | frames per second   |
        +===============+=====================+
        | 512 x 512     | 28.7                |
        +---------------+---------------------+
        | 350 x 350     | 61                  |
        +---------------+---------------------+
        | 300 x 300     | 82                  |
        +---------------+---------------------+
        | 250 x 250     | 115                 |
        +---------------+---------------------+
        | 200 x 200     | 180                 |
        +---------------+---------------------+
        | 150 x 150     | 300                 |
        +---------------+---------------------+
        | 100 x 100     | 640                 |
        +---------------+---------------------+ 

    background : {'auto', 'manual', None}, optional
        The parameter affects how the background count-rate estimation is done. 
        
        * ``'auto'``: Automatic estimation of the background count-rate.
        * ``'manual'``: To manually specify a background region using **x_bg** and **y_bg** parameters.
        * ``None``: No background estimation is carried out. This is the default method.
     
    sky_radius: float, optional
        The background aperture radius in pixels. 
        The default value is 12.
        
    x_bg : float, optional
        The X-coordinate of the background region. 
        
    y_bg : float, optional
        The Y-coordinate of the background region. 
        
    aperture_correction : {'fuv', 'nuv', None}, optional
        The parameter affects how the aperture correction is done. 
        
        * ``'fuv'``: Aperture correction for the FUV channel is applied. 
        * ``'nuv'``: Aperture correction for the NUV channel is applied.
        * ``None``: No aperture correction is applied. This is the default method.
        
    saturation_correction : bool, optional
        If `True`, saturation correction is applied. 
        The default value is `False`. 


    Note
    ---- 
    It is essential to set the correct value of the framerate. 
    Most UVIT observations are carried out in 512 x 512 window mode.
        
    Example
    --------
    >>> import curvit
    >>> curvit.makecurves(events_list = 'AS1G06_084T01_9000000710uvtNIIPC00F2_l2ce.fits.gz',
    ...                   threshold = 5)
    
    :: 
    
        Detected source coordinates saved in file:
        * sources_AS1G06_084T01_9000000710uvtNIIPC00F2_l2ce.coo
        Detected sources are plotted in the image:
        * sources_AS1G06_084T01_9000000710uvtNIIPC00F2_l2ce.png

        ---------------------- light curves ----------------------
        * makecurves_3136.64_3651.08_AS1G06_084T01_9000000710uvtNIIPC00F2_l2ce.png
        * makecurves_2530.02_1442.18_AS1G06_084T01_9000000710uvtNIIPC00F2_l2ce.png
        * makecurves_2912.31_3657.17_AS1G06_084T01_9000000710uvtNIIPC00F2_l2ce.png
        ...
        ...

        Done!

    """


    time, fx, fy, photons = read_columns(events_list)
    weights = photons / framecount_per_sec

    sanity = tobe_or_notobe(time, bwidth, 
                            detection_method,
                            threshold,
                            how_many, 
                            background, 
                            x_bg, y_bg,
                            aperture_correction, radius,
                            saturation_correction)
                           
    if sanity < 1:
        return

    original_input = events_list
    path_to_events_list, events_list = ntpath.split(events_list)
    events_list = modify_string(events_list)

    if detection_method == 'daofind':
        uA = detect_sources_daofind(fx, fy, photons, threshold)
    else:
        uA = detect_sources_kdtree(fx, fy, how_many)


    if len(uA) == 0:
        print('No sources, try changing the detection threshold parameter.')
        return
    
    coo_file = os.path.join(path_to_events_list, 'sources_' + events_list +'.coo')
    np.savetxt(coo_file, uA, fmt = '%4.2f\t%4.2f')
    print('\nDetected source coordinates saved in file:\n* {}'.format(coo_file))

    # To automatically choose background region.
    plt.figure(figsize = (10.5, 10))
    if background == 'auto':
        lowres_counts, bg_CPS, bg_CPS_e = auto_bg(fx, fy, time, photons, radius,
                                                  framecount_per_sec, sky_radius)
       
    # To create a quick look figure marking sources and background.
    bins = np.arange(0, 4801, 4096 / whole_figure_resolution)    
    plt.hist2d(fx, fy, 
               bins = (bins, bins),
               weights = weights,
               norm = LogNorm())
                   
    plt.tick_params(axis = 'both', labelsize = fontsize)

    for u in uA:
        plt.annotate('Source', u, 
                     size = 13, color = 'black', fontweight = 'bold')

        obj_circle = plt.Circle(u, 100, color = 'k', fill = False)
        plt.gcf().gca().add_artist(obj_circle)

    if background == 'manual':
        plt.annotate('Background', (x_bg, y_bg),
                     size = 13, color = 'black', fontweight = 'bold')

        bg_circle = plt.Circle((x_bg, y_bg), 100, color = 'k', fill = False)
        plt.gcf().gca().add_artist(bg_circle)
        
    png_name = os.path.join(path_to_events_list, 'sources_' + events_list + '.png')
    plt.savefig( png_name, format = 'png', bbox_inches = 'tight')
    plt.clf()

    print('Detected sources are plotted in the image:\n* {}'.format(png_name))   
 
    if background != None:
        if background == 'auto':
            print('\nThe estimated background CPS = {:.5f} +/-{:.5f}'.format(bg_CPS, bg_CPS_e))
        if background == 'manual':
            # To estimate Background CPS.
            bg_CPS, bg_CPS_e = bg_estimate(fx, fy, time, photons, framecount_per_sec,
                                           radius, x_bg, y_bg, sky_radius)
            bg_png = create_sub_image(x_bg, y_bg, 
                                      sub_fig_size, 
                                      sky_radius, 
                                      'background_',
                                      fx, fy,
                                      path_to_events_list,
                                      events_list)
            print('\nThe estimated background CPS = {:.5f} +/-{:.5f}'.format(bg_CPS, bg_CPS_e))
            print('Region selected for background estimate:\n* {}'.format(bg_png))
    else:
        bg_CPS, bg_CPS_e = 0, 0 

    # Calculating number of bins.
    time_width = time.max() - time.min() 
    nbin = time_width / bwidth
    nbin = int(nbin)

    unique_time = np.unique(time)
    till_here = time.min() + (bwidth * nbin)

    # Changing mission elapsed time in seconds to modified julian date.
    till_here = met_to_mjd(till_here)
    time_start = met_to_mjd(time.min())
    unique_time = met_to_mjd(unique_time)

    # Getting the number of unique frames within a bin.
    u_counts, u_bin_edges = np.histogram(unique_time, bins = nbin,
                                         range = (time_start, till_here))

    # selecting events within a circular region.
    print('\n---------------------- light curves ----------------------')
    plt.figure(figsize = (8, 5))
    for uaxy in uA:
        xp, yp = uaxy
        
        # selecting events within a circular region.
        mask = ((fx - xp) ** 2 + (fy - yp) ** 2) <= radius ** 2    
        T = time[mask]
        W = weights[mask]      
        T = met_to_mjd(T)         

        plt.title("X = %s, Y = %s, bin = %ss, radius = %spx" \
                  %(xp, yp, bwidth, radius), fontsize = fontsize)

        plt.xlabel("Time (Julian Date)", fontsize = fontsize)
        plt.ylabel("Counts per second", fontsize = fontsize)
        plt.tick_params(axis = 'both', labelsize = fontsize)

        weighted_counts, bin_edges = np.histogram(T, bins = nbin, 
                                                  range = (time_start, till_here),
                                                  weights = W)    

        counts, _ = np.histogram(T, bins = nbin, 
                                 range = (time_start, till_here))

        bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2.

        if np.array_equal(bin_edges, u_bin_edges) == True:
            count_mask = counts != 0
        else:
            print('\nThis happens when bwidth is too small\n')
            return

        if np.sum(count_mask) != 0:
            mcentres =  bin_centres[count_mask]
            weighted_mcounts = weighted_counts[count_mask]
            mcounts = counts[count_mask]
            frames_in_bin = u_counts[count_mask]
        else:
            print('No counts for the source at %s' %uaxy)
            continue

        CPF = weighted_mcounts / frames_in_bin
        CPF_err = np.sqrt(mcounts) / frames_in_bin
        
        # Background subtraction.
        CPF = CPF - (bg_CPS / framecount_per_sec)
        CPF_err = np.sqrt(CPF_err ** 2 + (bg_CPS_e / framecount_per_sec) ** 2)   

        CPF, CPF_err = apply_aperture_correction(CPF, CPF_err, radius, aperture_correction)
        CPF, CPF_err = apply_saturation_correction(CPF, CPF_err, saturation_correction)

        
        CPS = CPF * framecount_per_sec
        CPS_err = CPF_err * framecount_per_sec

        plt.scatter(mcentres, CPS)
        plt.errorbar(mcentres, CPS, yerr = CPS_err, linestyle = "None")
        empty_space = (till_here - time_start) / 25.0
        plt.xlim(time_start - empty_space, till_here + empty_space)
        
        #To write the array to output.
        data_to_output = list(zip(mcentres, CPS, CPS_err))
        output_prefix = 'makecurves_' + str(xp) + '_' + str(yp) + '_' + events_list
        datname = os.path.join(path_to_events_list, output_prefix + '.dat')
        np.savetxt(datname, data_to_output,
                   fmt = '%10.11f\t%.5e\t%.5e',
                   header = 'MJD\t\t\tCPS (bin=%ss)\tCPS_error' %bwidth)

        output_prefix = 'makecurves_' + str(xp) + '_' + str(yp) + '_' + events_list
        figname = os.path.join(path_to_events_list, output_prefix + '.png')
        plt.savefig(figname, format = 'png', bbox_inches = 'tight', dpi = 150)
        print('* {}'.format(figname))
        
        plt.clf()

    print('\nDone!\n')
    plt.close('all')
        
def curve(events_list = events_list,
          xp = xp,
          yp = yp,
          radius = radius,
          bwidth = bwidth,
          framecount_per_sec = framecount_per_sec,
          background = background,
          sky_radius = sky_radius,
          x_bg = x_bg,
          y_bg = y_bg,
          aperture_correction = aperture_correction,
          saturation_correction = saturation_correction,
          whole_figure_resolution = whole_figure_resolution,
          sub_fig_size = sub_fig_size,
          fontsize = fontsize):
          
          
    """Create light curve for a source.

    Parameters
    ----------
    events_list : file path
        The name of the events list FITS file.
        
    xp : float
        The X-coordinate of the source.
        
    yp : float
        The Y-coordinate of the source. 
               
    radius : float, optional
        The source aperture radius in pixels. 
        This parameter has a default value of 6.
        
    bwidth : float, optional
        Time bin width in seconds. 
        the default value is 50.   
        
    framecount_per_sec : float, optional
        The framerate of the observation, with a default value of 28.7185
        frames per second for 512 x 512 window mode. 
        The most accurate way to get the framerate would be to take the value 
        of (``1 / INT_TIME``). 
        ``INT_TIME`` value can be found from the corresponding image header. 
        Approximate values of framerate for different window modes of UVIT 
        are given in the table below.

        +---------------+---------------------+
        | window mode   | frames per second   |
        +===============+=====================+
        | 512 x 512     | 28.7                |
        +---------------+---------------------+
        | 350 x 350     | 61                  |
        +---------------+---------------------+
        | 300 x 300     | 82                  |
        +---------------+---------------------+
        | 250 x 250     | 115                 |
        +---------------+---------------------+
        | 200 x 200     | 180                 |
        +---------------+---------------------+
        | 150 x 150     | 300                 |
        +---------------+---------------------+
        | 100 x 100     | 640                 |
        +---------------+---------------------+ 

    background : {'auto', 'manual', None}, optional
        The parameter affects how the background count-rate estimation is done. 
        
        * ``'auto'``: Automatic estimation of the background count-rate.
        * ``'manual'``: To manually specify a background region using **x_bg** and **y_bg** parameters.
        * ``None``: No background estimation is carried out. This is the default method.
     
    sky_radius: float, optional
        The background aperture radius in pixels. 
        The default value is 12.
        
    x_bg : float, optional
        The X-coordinate of the background region. 
        
    y_bg : float, optional
        The Y-coordinate of the background region. 
        
    aperture_correction : {'fuv', 'nuv', None}, optional
        The parameter affects how the aperture correction is done. 
        
        * ``'fuv'``: Aperture correction for the FUV channel is applied. 
        * ``'nuv'``: Aperture correction for the NUV channel is applied.
        * ``None``: No aperture correction is applied. This is the default method.
        
    saturation_correction : bool, optional
        If `True`, saturation correction is applied. 
        The default value is `False`. 


    Note
    ---- 
    It is essential to set the correct value of the framerate. 
    Most UVIT observations are carried out in 512 x 512 window mode.
        
    Example
    --------
    >>> curvit.curve(events_list = 'AS1G06_084T01_9000000710uvtFIIPC00F1_l2ce.fits.gz', 
    ...              xp = 2636.71, yp = 907.91,
    ...              radius = 15,
    ...              bwidth = 50, 
    ...              background = 'auto')
    
    ::
    
        The estimated background CPS = 0.02155 +/-0.00425

        -------------------------- curve --------------------------
        source: source_AS1G06_084T01_9000000710uvtFIIPC00F1_l2ce.png
                source_zoomed_AS1G06_084T01_9000000710uvtFIIPC00F1_l2ce.png
        data: curve_2636.71_907.91_AS1G06_084T01_9000000710uvtFIIPC00F1_l2ce.dat
        plot: curve_2636.71_907.91_AS1G06_084T01_9000000710uvtFIIPC00F1_l2ce.png

        Done!
    
    """

    time, fx, fy, photons = read_columns(events_list)
    weights = photons / framecount_per_sec

    if None in [xp, yp]:
        print('\nPlease provide values for both "xp" and "yp".\n')
        return 
    
    sanity = tobe_or_notobe(time, bwidth, 
                            detection_method,
                            threshold,
                            how_many, 
                            background, 
                            x_bg, y_bg,
                            aperture_correction, radius,
                            saturation_correction)
    if sanity < 1:
        return

    path_to_events_list, events_list = ntpath.split(events_list)
    events_list = modify_string(events_list)

    # To automatically choose background region.
    plt.figure(figsize = (10.5, 10))
    if background == 'auto':
        lowres_counts, bg_CPS, bg_CPS_e = auto_bg(fx, fy, time, photons, radius,
                                                  framecount_per_sec, sky_radius)
        
    # To create a quick look figure marking sources and background.
    bins = np.arange(0, 4801, 4096 / whole_figure_resolution)    
    plt.hist2d(fx, fy, 
               bins = (bins, bins),
               weights = weights,
               norm = LogNorm())

    plt.tick_params(axis = 'both', labelsize = fontsize)

    plt.annotate("Source", (xp, yp), 
                 size = 13, color = 'black', fontweight = 'bold')

    obj_circle = plt.Circle((xp, yp), 100, color = 'k', fill = False)
    plt.gcf().gca().add_artist(obj_circle)

    if background == 'manual':
        plt.annotate('Background', (x_bg, y_bg),
                     size = 13, color = 'black', fontweight = 'bold')

        bg_circle = plt.Circle((x_bg, y_bg), 100, color = 'k', fill = False)
        plt.gcf().gca().add_artist(bg_circle)
        
    png_name = os.path.join(path_to_events_list, 'source_' + events_list + '.png')
    plt.savefig(png_name, format = 'png', bbox_inches = 'tight')
    plt.clf()     
        
    source_png = create_sub_image(xp, yp, 
                                  sub_fig_size, 
                                  radius, 
                                  'source_zoomed_',
                                  fx, fy,
                                  path_to_events_list,
                                  events_list)
    
    if background != None:
        if background == 'auto':
            print('\nThe estimated background CPS = {:.5f} +/-{:.5f}'.format(bg_CPS, bg_CPS_e))
        if background == 'manual':
            # To estimate Background CPS.
            bg_CPS, bg_CPS_e = bg_estimate(fx, fy, time, photons, framecount_per_sec,
                                           radius, x_bg, y_bg, sky_radius)
            bg_png = create_sub_image(x_bg, y_bg, 
                                      sub_fig_size, 
                                      sky_radius, 
                                      'background_',
                                      fx, fy,
                                      path_to_events_list,
                                      events_list)
            print('\nThe estimated background CPS = {:.5f} +/-{:.5f}'.format(bg_CPS, bg_CPS_e))
            print('Region selected for background estimate:\n* {}'.format(bg_png))
    else:
        bg_CPS, bg_CPS_e = 0, 0 
        
    # selecting events within a circular region.
    mask = ((fx - xp) ** 2 + (fy - yp) ** 2) <= radius ** 2    
    T = time[mask]
    W = weights[mask]  

    # Calculating number of bins.
    time_width = time.max() - time.min() 
    nbin = time_width / bwidth
    nbin = int(nbin)

    unique_time = np.unique(time)
    till_here = time.min() + (bwidth * nbin)

    # Changing mission elapsed time in seconds to modified julian date.
    T = met_to_mjd(T)
    till_here = met_to_mjd(till_here)
    time_start = met_to_mjd(time.min())
    unique_time = met_to_mjd(unique_time)

    # Binning stuff, plotting stuff.
    plt.figure(figsize = (8, 5))
    plt.title('bin = %ss, radius = %spx' %(bwidth, radius), fontsize = fontsize)
    plt.xlabel("Time (Julian Date)", fontsize = fontsize)
    plt.ylabel("Counts per second", fontsize = fontsize)
    plt.tick_params(axis = 'both', labelsize = fontsize)

    u_counts, u_bin_edges = np.histogram(unique_time, bins = nbin,
                                         range = (time_start, till_here))

    weighted_counts, bin_edges = np.histogram(T, bins = nbin, 
                                              range = (time_start, till_here),
                                              weights = W)    

    counts, _ = np.histogram(T, bins = nbin, 
                             range = (time_start, till_here))

    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2.

    if np.array_equal(bin_edges, u_bin_edges) == True:
        count_mask = counts != 0
    else:
        print('\nThis happens when bwidth is too small\n')
        return

    if np.sum(count_mask) != 0:
        mcentres =  bin_centres[count_mask]
        weighted_mcounts = weighted_counts[count_mask]
        mcounts = counts[count_mask]
        frames_in_bin = u_counts[count_mask]
    else:
        print('No counts inside the aperture!')

    CPF = weighted_mcounts / frames_in_bin
    CPF_err = np.sqrt(mcounts) / frames_in_bin
    
    # Background subtraction.
    CPF = CPF - (bg_CPS / framecount_per_sec)
    CPF_err = np.sqrt(CPF_err ** 2 + (bg_CPS_e / framecount_per_sec) ** 2)    

    CPF, CPF_err = apply_aperture_correction(CPF, CPF_err, radius, aperture_correction)
    CPF, CPF_err = apply_saturation_correction(CPF, CPF_err, saturation_correction)
    
    CPS = CPF * framecount_per_sec
    CPS_err = CPF_err * framecount_per_sec

    # Let us get on with the plot. 
    plt.scatter(mcentres, CPS)
    plt.errorbar(mcentres, CPS, yerr = CPS_err, linestyle = "None")
    # To make the plot look good.
    empty_space = (till_here - time_start) / 25.0
    plt.xlim(time_start - empty_space, till_here + empty_space)

    #To write the array to output.
    data_to_output = list(zip(mcentres, CPS, CPS_err))
    output_prefix = 'curve_' + str(xp) + '_' + str(yp) + '_' + events_list
    datname = os.path.join(path_to_events_list, output_prefix + '.dat')
    np.savetxt(datname, data_to_output,
               fmt = '%10.11f\t%.5e\t%.5e',
               header = 'MJD\t\t\tCPS (bin=%ss)\tCPS_error' %bwidth)

    figname = os.path.join(path_to_events_list, output_prefix + '.png')
    plt.savefig(figname, format = 'png', bbox_inches = 'tight', dpi = 150)

    print('\n-------------------------- curve --------------------------')
    print('source: {}\n        {}'.format(png_name, source_png))
    print('data: {}'.format(datname))
    print('plot: {}'.format(figname))

    print("\nDone!\n")
    plt.close('all') 
    
# Function to convert CCDLAB XYFrac and XYInts to X, Y positions in 4k.
def CCDLAB_to_4k(Int, Frac):
    coo_in_4k = ((Int + Frac) - 16) / 4.0
    return coo_in_4k
    
# Function to convert CCDLAB files to a compatible events list.     
def process_ccdlab(output = None,
                   time_list = None, 
                   XY_integers = None, 
                   XY_fractions = None, 
                   flat_list = None, 
                   framecount_per_sec = framecount_per_sec):
                   
    """Generate a Curvit compatible events list from CCDLAB files.

    Parameters
    ----------
    output : file path
        The name of the output events list FITS file.
        
    time_list : file path
        The name of the CCDLAB time list FITS file
        
    XY_integers : file path
        The name of the CCDLAB XY integers FITS file
        
    XY_fractions : file path
        The name of the CCDLAB XY fractions FITS file
        
    flat_list : file path
        The name of the CCDLAB flat list FITS file
        
    framecount_per_sec : float, optional
        The framerate of the observation, with a default value of 28.7185
        frames per second for 512 x 512 window mode. 
        The most accurate way to get the framerate would be to take the value 
        of (``1 / INT_TIME``). 
        ``INT_TIME`` value can be found from the corresponding image header. 
        Approximate values of framerate for different window modes of UVIT 
        are given in the table below.

        +---------------+---------------------+
        | window mode   | frames per second   |
        +===============+=====================+
        | 512 x 512     | 28.7                |
        +---------------+---------------------+
        | 350 x 350     | 61                  |
        +---------------+---------------------+
        | 300 x 300     | 82                  |
        +---------------+---------------------+
        | 250 x 250     | 115                 |
        +---------------+---------------------+
        | 200 x 200     | 180                 |
        +---------------+---------------------+
        | 150 x 150     | 300                 |
        +---------------+---------------------+
        | 100 x 100     | 640                 |
        +---------------+---------------------+ 
        

    Note
    ---- 
    It is essential to set the correct value of the framerate. 
    Most UVIT observations are carried out in 512 x 512 window mode.
            
    Warning
    -------
    This function is new; please report if you find any bugs.
        
    Example
    --------
    >>> import curvit
    >>> process_ccdlab(output = 'output_events_list.fits',
    ...                time_list = 'sample_TimeList.fits', 
    ...                XY_integers = 'sample_XYInts_List.fits',
    ...                XY_fractions = 'sample_XYFrac_List.fits',
    ...                flat_list = 'sample_FlatList.fits',
    ...                framecount_per_sec = 28.7185)
    
    The above script will generate a FITS table called ``output_events_list.fits``.
    You may then use it as input to ``curve`` or ``makecurves``. 
    """
    
    time = fits.open(time_list)[0].data / 1000
    XYFrac = fits.open(XY_fractions)[0].data
    XYInts = fits.open(XY_integers)[0].data
    weight = fits.open(flat_list)[0].data 
    photons = weight * framecount_per_sec
    fx = CCDLAB_to_4k(XYInts[:,0], XYFrac[:,0])
    fy = CCDLAB_to_4k(XYInts[:,1], XYFrac[:,1])
    
    col1 = fits.Column(name = 'MJD_L2', format = 'D', array = time)
    col2 = fits.Column(name = 'Fx', format = 'D', array = fx)
    col3 = fits.Column(name = 'Fy', format = 'D', array = fy)
    col4 = fits.Column(name = 'EFFECTIVE_NUM_PHOTONS', format = 'D', array = photons)
    col5 = fits.Column(name = 'BAD FLAG', format = 'D', array = np.ones(len(time)))

    cols = fits.ColDefs([col1, col2, col3, col4, col5])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(output, overwrite = True) 
    print('The events list: {}'.format(output))  

def makefits(events_list = events_list, 
             framecount_per_sec = framecount_per_sec):

    """Create a FITS image from the input events list.

    Parameters
    ----------
    events_list : file path
        The name of the events list FITS file.
        
    framecount_per_sec : float, optional
        The framerate of the observation, with a default value of 28.7185
        frames per second for 512 x 512 window mode. 
        The most accurate way to get the framerate would be to take the value 
        of (``1 / INT_TIME``). 
        ``INT_TIME`` value can be found from the corresponding image header. 
        Approximate values of framerate for different window modes of UVIT 
        are given in the table below.

        +---------------+---------------------+
        | window mode   | frames per second   |
        +===============+=====================+
        | 512 x 512     | 28.7                |
        +---------------+---------------------+
        | 350 x 350     | 61                  |
        +---------------+---------------------+
        | 300 x 300     | 82                  |
        +---------------+---------------------+
        | 250 x 250     | 115                 |
        +---------------+---------------------+
        | 200 x 200     | 180                 |
        +---------------+---------------------+
        | 150 x 150     | 300                 |
        +---------------+---------------------+
        | 100 x 100     | 640                 |
        +---------------+---------------------+ 

    Warning
    -------
    If you plan to use the generated FITS image for science,
    make sure to give the proper framerate value. 
        
    Example
    --------
    >>> import curvit
    >>> curvit.makefits('test_events_list.fits', 28.7185)
    
    The above script will generate a FITS image called ``test_events_list_image.fits``.
    You may open it in software such as DS9 to view the image. 
    """

    time, fx, fy, photons = read_columns(events_list)
    weights = photons / framecount_per_sec
    bins = np.arange(0, 4801)
    ndarray, yedges, xedges = np.histogram2d(fy, fx, bins = (bins, bins), weights = weights)  
    fits_name = events_list.replace('.fits', '_image.fits')
    hdu = fits.PrimaryHDU(data = ndarray)
    try:
        hdu.header['RA_PNT'] = fits.open(events_list)[1].header['RA_PNT_RADEC']
        hdu.header['DEC_PNT'] = fits.open(events_list)[1].header['DEC_PNT_RADEC']
    except KeyError:
        pass
    
    try:
        hdu.header['RA_PNT'] = fits.open(events_list)[0].header['RA_PNT']
        hdu.header['DEC_PNT'] = fits.open(events_list)[0].header['DEC_PNT']
    except KeyError:
        pass
        
    try:
        hdu.header['EXPTIME'] = fits.open(events_list)[0].header['EXPTIME']
    except KeyError:
        pass
    hdu.writeto(fits_name, overwrite = True)      
    print('The image: {}'.format(fits_name))    
    
def rebin(arr, bin_factor):
    shape = (int(arr.shape[0] / bin_factor), bin_factor,
             int(arr.shape[1] / bin_factor), bin_factor)
    binned_arr = arr.reshape(shape).mean(-1).mean(1)
    return binned_arr

def get_image_data(fx, fy, photons, framecount_per_sec):
    weights = photons / framecount_per_sec
    bins = np.arange(0, 4801)
    data, _, _ = np.histogram2d(fy, fx, bins = (bins, bins), weights = weights)
    return data

def daofind_on_image_data(data, threshold):
    kernel = Gaussian2DKernel(x_stddev=1.5)
    binned_data = rebin(data, 2)
    smoothed_data = convolve(binned_data, kernel)
    zoomed_data = zoom(smoothed_data, zoom = 2, order = 0)
    
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = MedianBackground()
    bkg = Background2D(zoomed_data, (50, 50), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

    daofind = DAOStarFinder(fwhm = 4, 
                            threshold = threshold * bkg.background_rms_median, 
                            exclude_border = True)
    
    sources = daofind(zoomed_data - bkg.background)    
    
    sources.sort('mag')
    uA = np.array([sources['xcentroid'].data, sources['ycentroid'].data]).T
    return uA 

def new_detect_sources_daofind(fx, fy, photons, threshold, framecount_per_sec):
    data = get_image_data(fx, fy, photons, framecount_per_sec)
    uA = daofind_on_image_data(data, threshold)
    return uA

def combine_events_lists(events_lists_paths = None, 
                         shift_algorithm = shift_algorithm, 
                         min_exptime = min_exptime,
                         framecount_per_sec = framecount_per_sec):
                         
    """Align and combine the events lists from different orbits.

    Parameters
    ----------
    events_lists_paths : list
        The list of events list FITS file paths.
        
    shift_algorithm : {'single_star', 'multiple_star'}, optional
        The parameter to choose between available aligning methods. 
        
        * ``'single_star'``: To use a single star for aligning orbits. 
            Single_star option is useful only when there is only one star 
            in the field and no rotation between frames.
        * ``'multiple_star'``: To use multiple stars for aligning orbits. 
            This is recommended and the default method.
                               
    min_exptime : float, optional
        Orbits having exposure time below this limit will be ignored. 
        the default value is 100 seconds. 
                
    framecount_per_sec : float, optional
        The framerate of the observation.
        
    Warning
    -------
    While combining events lists, do not mix data from ``RAS_VIS`` 
    and ``RAS_NUV``. In most cases, ``RAS_VIS`` data should be preferred. 
    
    Note
    ---- 
    Please cite the Astroalign package if you are using this function.     
    
    Example
    --------
    >>> import curvit
    >>> from glob import glob
    >>> file_path_list = glob('*/*/*/RAS_VIS/*/F*/*F1*ce*')
    >>> curvit.combine_events_lists(file_path_list)
    """
    
    print('\nWorking on the following events lists:')
    print(*events_lists_paths, sep='\n')
    print('')

    exptimes = []
    detected_sources = []
    select_events_lists_paths = []    
    for path in events_lists_paths:
        time, fx, fy, photons = read_columns(path)
        nbin_check = (time.max() - time.min()) / min_exptime
        if nbin_check < 1:
            print('\n{} ignored.\nExposure time below {} seconds.'.format(path, min_exptime))
            continue
            
        uA = []    
        if shift_algorithm == 'single_star':
            number_of_sources = 2
            i = 6
            while len(uA) <= number_of_sources:
                uA = new_detect_sources_daofind(fx, fy, photons, i, framecount_per_sec)
                i = i - 0.5
                if i <= 1:
                    print('If you see this, please contact Curvit developer.')
                    break
        else:
            number_of_sources = 3
            i = 5
            while len(uA) <= number_of_sources:
                uA = new_detect_sources_daofind(fx, fy, photons, i, framecount_per_sec)
                if len(uA) > 200:
                    uA = new_detect_sources_daofind(fx, fy, photons, i + 5, framecount_per_sec)
                i = i - 0.5
                if i <= 1:
                    print('If you see this, please contact Curvit developer.')
                    break  

        print('{} sources detected in {}'.format(len(uA), path))

        exptimes.append(time.max() - time.min())
        detected_sources.append(uA)
        select_events_lists_paths.append(path)
        
    print('')
      
    framerate_from_header = []

    if shift_algorithm == 'multiple_star':    
        _, eventslist_name = ntpath.split(select_events_lists_paths[0])
        eventslist_name = modify_string(eventslist_name)
        combined_eventslist_name = eventslist_name + '_all_orbits.fits'
        reference = detected_sources[np.argmax(exptimes)]
        i = 0
        for elist, path in zip(detected_sources, select_events_lists_paths):
            try:
                transf, (source_list, target_list) = aa.find_transform(elist, reference, max_control_points = 100)  

                img_path = ntpath.split(path)[0] + '/*I_l2img*'
                if len(glob(img_path)) == 1:
                    img_hdu = fits.open(glob(img_path)[0])
                    framerate_from_header.append(1 / img_hdu[0].header['INT_TIME'])
                    RA_pointing = img_hdu[0].header['RA_PNT']
                    DEC_pointing = img_hdu[0].header['DEC_PNT']
                else:
                    print('Requires exactly one file matching *I_l2img* in directory {}!'.format(ntpath.split(path)[0]))
                    sys.exit()

                if i == 0:
                    hdu_base = fits.open(path)

                    photons = hdu_base[1].data['EFFECTIVE_NUM_PHOTONS']
                    mask = photons > 0 

                    nrows_base = hdu_base[1].data[mask].shape[0]

                    X_and_Y = np.array([hdu_base[1].data['Fx'], hdu_base[1].data['Fy']]).T
                    X_centroid, Y_centroid = aa.matrix_transform(X_and_Y, transf.params).T        

                    hdu_base[1].data['Fx'] = X_centroid
                    hdu_base[1].data['Fy'] = Y_centroid

                    hdu = fits.BinTableHDU.from_columns(hdu_base[1].columns, nrows = nrows_base, fill = True)
                    for colname in hdu_base[1].columns.names:
                        hdu.data[colname][:nrows_base] = hdu_base[1].data[colname][mask]

                    hduP = fits.PrimaryHDU()
                    hduP.header = hdu_base[0].header
                    hdu.name = hdu_base[1].name
                    hdu_base = fits.HDUList([hduP, hdu])
                    hdu_base.writeto(combined_eventslist_name, overwrite = True)

                else:
                    hdu_add = fits.open(path)
                    hdu_base = fits.open(combined_eventslist_name)

                    photons = hdu_add[1].data['EFFECTIVE_NUM_PHOTONS']
                    mask = photons > 0 

                    nrows_base = hdu_base[1].data.shape[0]
                    nrows_add = hdu_add[1].data[mask].shape[0]
                    nrows = nrows_base + nrows_add

                    X_and_Y = np.array([hdu_add[1].data['Fx'], hdu_add[1].data['Fy']]).T
                    X_centroid, Y_centroid = aa.matrix_transform(X_and_Y, transf.params).T 

                    hdu_add[1].data['Fx'] = X_centroid
                    hdu_add[1].data['Fy'] = Y_centroid                

                    hdu = fits.BinTableHDU.from_columns(hdu_base[1].columns, nrows = nrows, fill = True)
                    for colname in hdu_base[1].columns.names:
                        hdu.data[colname][:nrows_base] = hdu_base[1].data[colname]
                        hdu.data[colname][nrows_base:] = hdu_add[1].data[colname][mask]

                    hdu.name = hdu_base[1].name
                    hdu_base = fits.HDUList([hduP, hdu])
                    hdu_base.writeto(combined_eventslist_name, overwrite = True)

                i = i + 1
                print('coordinate matching successful for {}'.format(path))

            except aa.MaxIterError:
                print('coordinate matching failed for {}'.format(path))


    if shift_algorithm == 'single_star':
        new_detected_sources = np.array([d[0] for d in detected_sources])
        shifts = new_detected_sources - new_detected_sources[0]

        _, eventslist_name = ntpath.split(select_events_lists_paths[0])
        eventslist_name = modify_string(eventslist_name)
        combined_eventslist_name = eventslist_name + '_all_orbits.fits'
        i = 0
        for shift, path in zip(shifts, select_events_lists_paths):

            img_path = ntpath.split(path)[0] + '/*I_l2img*'
            if len(glob(img_path)) == 1:
                img_hdu = fits.open(glob(img_path)[0])
                framerate_from_header.append(1 / img_hdu[0].header['INT_TIME'])
                RA_pointing = img_hdu[0].header['RA_PNT']
                DEC_pointing = img_hdu[0].header['DEC_PNT']
            else:
                print('Requires exactly one file matching *I_l2img* in directory {}!'.format(ntpath.split(path)[0]))
                sys.exit()

            if i == 0:
                hdu_base = fits.open(path)

                photons = hdu_base[1].data['EFFECTIVE_NUM_PHOTONS']
                mask = photons > 0 

                nrows_base = hdu_base[1].data[mask].shape[0]

                hdu = fits.BinTableHDU.from_columns(hdu_base[1].columns, nrows = nrows_base, fill = True)
                for colname in hdu_base[1].columns.names:
                    hdu.data[colname][:nrows_base] = hdu_base[1].data[colname][mask]

                hduP = fits.PrimaryHDU()
                hduP.header = hdu_base[0].header
                hdu.name = hdu_base[1].name
                hdu_base = fits.HDUList([hduP, hdu])
                hdu_base.writeto(combined_eventslist_name, overwrite = True)


            else:
                hdu_add = fits.open(path)
                hdu_base = fits.open(combined_eventslist_name)


                photons = hdu_add[1].data['EFFECTIVE_NUM_PHOTONS']
                mask = photons > 0 

                nrows_base = hdu_base[1].data.shape[0]
                nrows_add = hdu_add[1].data[mask].shape[0]
                nrows = nrows_base + nrows_add

                hdu_add[1].data['Fx'] = hdu_add[1].data['Fx'] - shift[0]
                hdu_add[1].data['Fy'] = hdu_add[1].data['Fy'] - shift[1]

                hdu = fits.BinTableHDU.from_columns(hdu_base[1].columns, nrows = nrows, fill = True)
                for colname in hdu_base[1].columns.names:
                    hdu.data[colname][:nrows_base] = hdu_base[1].data[colname]
                    hdu.data[colname][nrows_base:] = hdu_add[1].data[colname][mask]

                hdu.name = hdu_base[1].name
                hdu_base = fits.HDUList([hduP, hdu])
                hdu_base.writeto(combined_eventslist_name, overwrite = True)

            i = i + 1

    AVGFRMRT = np.mean(framerate_from_header)
    STDFRMRT = np.std(framerate_from_header)

    print('\nAverage frame rate = {:.6f}'.format(AVGFRMRT))
    print('frame rate standard deviation = {:.6f}'.format(STDFRMRT))

    hdu_base = fits.open(combined_eventslist_name, mode = 'update')
    hdu_base[0].header['AVGFRMRT'] = (AVGFRMRT, 'Added by events lists combine software')
    hdu_base[0].header['STDFRMRT'] = (STDFRMRT, 'Added by events lists combine software')
    hdu_base[0].header['RA_PNT'] = RA_pointing
    hdu_base[0].header['DEC_PNT'] = DEC_pointing
    time = hdu_base[1].data['MJD_L2']
    photons = hdu_base[1].data['EFFECTIVE_NUM_PHOTONS']
    bad_flag = hdu_base[1].data['BAD FLAG']
    mask = photons > 0 
    mask = np.logical_and(mask, bad_flag)
    time = time[mask]
    unique_frames = len(np.unique(time))
    exp_time = unique_frames / AVGFRMRT 
    print('Total combined exposure time = {:.4f} seconds'.format(exp_time))
    hdu_base[0].header['EXPTIME'] = exp_time
    hdu_base.flush()
    print('The combined events list: {}'.format(combined_eventslist_name))    
    print("\nDone!\n")
    
def image_astrometry(UV_image = None, threshold = 3, API_key = AstrometryNet_API_key):

    """Carry out astrometry on a UVIT image using Astrometry.net.

    Parameters
    ----------
    UV_image : file path
        The name of the UVIT FITS image.
                               
    threshold : float, optional
        The threshold parameter associated with the source detection method. 
        The default value is 3.
                
    API_key : string, optional
        The Astrometry.net API key. Ideally, you should get your API key
        from Astrometry.net and use it. 

    Warning
    -------
    Astrometry should be successful on most fields. 
    However, failures found during tests on some crowded fields. 
    Please try changing the source detection threshold in such cases.      
        
    Note
    ---- 
    Please cite Astrometry.net if you are using this function.        
        
    Example
    --------
    >>> import curvit
    >>> curvit.image_astrometry('test.fits')        
    """    
    hdu = fits.open(UV_image)
    sources = daofind_on_image_data(hdu[0].data, threshold)
    sources = sources + 1
    print('Number of detected sources = {}'.format(len(sources)))
    ra_pnt = hdu[0].header['RA_PNT']
    dec_pnt = hdu[0].header['DEC_PNT']
    
    ast = AstrometryNet()
    ast.api_key = API_key

    try:
        wcs_header = ast.solve_from_source_list(sources[:, 0], sources[:, 1],
                                    hdu[0].header['NAXIS1'], hdu[0].header['NAXIS2'],
                                    solve_timeout = 1200, 
                                    center_ra = ra_pnt, 
                                    center_dec = dec_pnt,
                                    radius = 0.33,
                                    crpix_center = True,
                                    scale_type='ev',
                                    scale_units = 'arcsecperpix',
                                    scale_est = 0.41625,
                                    scale_err = 10)
    except TimeoutError:
        print('\nTimed out!')

    if wcs_header:
        print('\nAstrometry.net solve success!')
        
        # For fixing the UVIT L2 pipeline image header.
        try: 
            hdu[0].header.remove('CTYPE1')
            hdu[0].header.remove('CUNIT1') 
            hdu[0].header.remove('CRPIX1') 
            hdu[0].header.remove('CDELT1') 
            hdu[0].header.remove('CRVAL1') 
            hdu[0].header.remove('CTYPE2') 
            hdu[0].header.remove('CUNIT2') 
            hdu[0].header.remove('CRPIX2') 
            hdu[0].header.remove('CDELT2') 
            hdu[0].header.remove('CRVAL2') 
            hdu[0].header.remove('CROTA2') 
            hdu[0].header.remove('CROTA1')
        except KeyError:
            pass 

        try: 
            hdu[0].header.remove('CTYPE1')
            hdu[0].header.remove('CUNIT1') 
            hdu[0].header.remove('CRPIX1') 
            hdu[0].header.remove('CDELT1') 
            hdu[0].header.remove('CRVAL1') 
            hdu[0].header.remove('CTYPE2') 
            hdu[0].header.remove('CUNIT2') 
            hdu[0].header.remove('CRPIX2') 
            hdu[0].header.remove('CDELT2') 
            hdu[0].header.remove('CRVAL2') 
            hdu[0].header.remove('CROTA2') 
            hdu[0].header.remove('CROTA1')
        except KeyError:
            pass 
            
        hdu[0].header.update(wcs_header)
        hdu.writeto(UV_image, overwrite = True)
        print('Image header updated with WCS.')
    else:
        print('\nAstrometry.net solve FAILED!!!\n')        
