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
import matplotlib.pyplot as plt

from glob import glob
from astropy.io import fits
from collections import Counter
from scipy.spatial import KDTree
from scipy.interpolate import interp1d
from matplotlib.colors import LogNorm
from photutils import DAOStarFinder, CircularAperture
from astropy.stats import sigma_clipped_stats, gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel, convolve

 
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
#######################################################################


def read_columns(events_list):
    # Reading few columns.
    f = fits.open(events_list)
    time = f[1].data['MJD_L2']
    fx = f[1].data['Fx']
    fy = f[1].data['Fy']
    photons = f[1].data['EFFECTIVE_NUM_PHOTONS']
    mask = photons > 0 
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
def auto_bg(fx, fy, photons, framecount_per_sec): 
    weights = photons / framecount_per_sec    
    
    lowres_counts, lowres_xedges, \
    lowres_yedges, lowres_Image = plt.hist2d(fx, fy, 
                                             bins = 256,
                                             weights = weights,
                                             norm = LogNorm())

    lowres_xcentres = (lowres_xedges[:-1] + lowres_xedges[1:]) / 2.
    lowres_ycentres = (lowres_yedges[:-1] + lowres_yedges[1:]) / 2.
    flat_counts = lowres_counts.flatten()
    x_mesh, y_mesh = np.meshgrid(lowres_ycentres, lowres_xcentres) #Notice the swap.
    x_mesh = x_mesh.flatten() 
    y_mesh = y_mesh.flatten()
    # To avoid the edges. 
    mask = (x_mesh - 2400) ** 2 + (y_mesh - 2400) ** 2 < 1800 ** 2
    array = np.array([flat_counts, x_mesh, y_mesh]).T 
    polished_array = array[mask]
    sorted_counts = np.sort(polished_array, axis= 0)   

    five_percent = int(0.05 * len(sorted_counts))
    r_count, x_bg, y_bg = random.choice(sorted_counts[:five_percent])
    return lowres_Image, r_count, x_bg, y_bg   

# To estimate background CPS.
def bg_estimate(fx, fy, time, photons, framecount_per_sec, x_bg, y_bg, sky_radius):

    weights = photons / framecount_per_sec    
    T_W = [(t_b, w) for xx_b, yy_b, t_b, w in zip(fx, fy, time, weights)
                     if ((xx_b - x_bg)**2 + (yy_b - y_bg)**2) <= sky_radius**2]

    T, W = np.array(T_W).T

    if len(T) != 0:
        scaled_events = (np.sum(W) * radius**2) / float(sky_radius**2)
        scaled_events_e = (np.sqrt(len(T)) * radius**2) / float(sky_radius**2)
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

    obj_xy = [(obj_x, obj_y) for obj_x, obj_y 
              in zip(fx, fy) 
              if pos_x - sub_size <= obj_x <= pos_x + sub_size
                  and pos_y - sub_size <= obj_y <= pos_y + sub_size]

    obj_fx = np.array(obj_xy)[:,0]
    obj_fy = np.array(obj_xy)[:,1]
    obj_circle = plt.Circle((pos_x, pos_y), cir_rad, 
                            color = 'k', fill = False)

    plt.hist2d(obj_fx, obj_fy, bins = sub_size*2, norm = LogNorm())
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

    mean, median, std = sigma_clipped_stats(data, sigma = 5., maxiters = 1)
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
    T_W = [(t, w) for xx, yy, t, w in zip(fx, fy, time, weights) 
                   if ((xx - xp)**2 +(yy - yp)**2) <= radius**2]

    T, W = np.array(T_W).T 

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
               threshold = threshold,
               how_many = how_many,
               bwidth = bwidth,
               framecount_per_sec = framecount_per_sec,
               background = background,
               sky_radius = sky_radius,
               x_bg = x_bg,
               y_bg = y_bg,
               detection_method = detection_method,
               aperture_correction = aperture_correction,
               saturation_correction = saturation_correction,
               whole_figure_resolution = whole_figure_resolution,
               sub_fig_size = sub_fig_size, 
               fontsize = fontsize):

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
        print('No sources, try increasing the "how_many" parameter.')
        return
    
    coo_file = os.path.join(path_to_events_list, 'sources_' + events_list +'.coo')
    np.savetxt(coo_file, uA, fmt = '%4.f\t%4.f')
    print('\nDetected source coordinates saved in file:\n* {}'.format(coo_file))

    # To automatically choose background region.
    plt.figure(figsize = (10.5, 10))
    if background == 'auto':
        lowres_Image, r_count, x_bg, y_bg = auto_bg(fx, fy, photons, framecount_per_sec)
        plt.gcf().gca().add_artist(lowres_Image)
       
    # To create a quick look figure marking sources and background.
    if whole_figure_resolution != 256 or background != 'auto':  # To avoid doing this twice.
        plt.hist2d(fx, fy, bins = whole_figure_resolution,
                   weights = weights, 
                   norm = LogNorm())

    for u in uA:
        plt.annotate('Source', u, 
                     size = 13, color = 'black', fontweight = 'bold')

        obj_circle = plt.Circle(u, 100, color = 'k', fill = False)
        plt.gcf().gca().add_artist(obj_circle)

    if background != None:
        plt.annotate('Background', (x_bg, y_bg),
                     size = 13, color = 'black', fontweight = 'bold')

        bg_circle = plt.Circle((x_bg, y_bg), 100, color = 'k', fill = False)
        plt.gcf().gca().add_artist(bg_circle)
        
    png_name = os.path.join(path_to_events_list, 'sources_' + events_list + '.png')
    plt.savefig( png_name, format = 'png', bbox_inches = 'tight')
    plt.clf()

    print('Detected sources are plotted in the image:\n* {}'.format(png_name))   
 
    if background != None:
        bg_png = create_sub_image(x_bg, y_bg, 
                                  sub_fig_size, 
                                  sky_radius, 
                                  'background_',
                                  fx, fy,
                                  path_to_events_list,
                                  events_list)

        # To estimate Background CPS.
        bg_CPS, bg_CPS_e = bg_estimate(fx, fy, time, photons, 
                                       framecount_per_sec, x_bg, y_bg, sky_radius)
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
        T_W = np.array([(t, w) for xx, yy, t, w 
                                in zip(fx, fy, time, weights)
                                if ((xx - xp)**2 +(yy - yp)**2) <= radius**2])

        T, W = np.array(T_W).T 
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
            hisdata = [(mb, mc, mwc, uc) for mb, mwc, mc, uc 
                       in zip(bin_centres, weighted_counts, counts, u_counts) 
                       if mc != 0]
        else:
            print('\nThis happens when bwidth is too small\n')
            return

        if len(hisdata) != 0:
            mcentres, weighted_mcounts, mcounts, frames_in_bin = np.array(hisdata).T
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
        lowres_Image, r_count, x_bg, y_bg = auto_bg(fx, fy, photons, framecount_per_sec)
        plt.gcf().gca().add_artist(lowres_Image)
        
    # To create a quick look figure marking source and background.
    if whole_figure_resolution != 256 or background != 'auto':  # To avoid doing this twice.
        plt.hist2d(fx, fy, bins = whole_figure_resolution, 
                   weights = weights,
                   norm = LogNorm())

    plt.annotate("Source", (xp, yp), 
                 size = 13, color = 'black', fontweight = 'bold')

    obj_circle = plt.Circle((xp, yp), 100, color = 'k', fill = False)
    plt.gcf().gca().add_artist(obj_circle)

    if background != None:
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
        bg_png = create_sub_image(x_bg, y_bg, 
                                  sub_fig_size, 
                                  sky_radius, 
                                  'background_',
                                  fx, fy,
                                  path_to_events_list,
                                  events_list)

        # To estimate Background CPS.
        bg_CPS, bg_CPS_e = bg_estimate(fx, fy, time, photons, 
                                       framecount_per_sec, x_bg, y_bg, sky_radius)
        print('\nThe estimated background CPS = {:.5f} +/-{:.5f}'.format(bg_CPS, bg_CPS_e))
        print('Region selected for background estimate:\n* {}'.format(bg_png))
    else:
        bg_CPS, bg_CPS_e = 0, 0 
        
    # selecting events within a circular region.
    T_W = np.array([(t, w) for xx, yy, t, w 
                            in zip(fx, fy, time, weights)
                            if ((xx - xp)**2 +(yy - yp)**2) <= radius**2])

    T, W = np.array(T_W).T 

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
            hisdata = [(mb, mc, mwc, uc) for mb, mwc, mc, uc 
                       in zip(bin_centres, weighted_counts, counts, u_counts) 
                       if mc != 0]
    else:
        print('\nThis happens when bwidth is too small\n')
        return

    if len(hisdata) != 0:
        mcentres, weighted_mcounts, mcounts, frames_in_bin = np.array(hisdata).T
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

    cols = fits.ColDefs([col1, col2, col3, col4])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(output, overwrite = True)   
    return   

def makefits(events_list):
    time, fx, fy, photons = read_columns(events_list)
    weights = photons / framecount_per_sec
    bins = np.arange(0, 4801)
    ndarray, yedges, xedges = np.histogram2d(fy, fx, bins = (bins, bins), weights = weights)  
    fits_name = events_list.replace('.fits', '_quick_look.fits')
    hdu = fits.PrimaryHDU(data = ndarray)
    hdu.writeto(fits_name, overwrite = True)        
