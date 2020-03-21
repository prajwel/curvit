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
from matplotlib.colors import LogNorm


#def get_eventsfile():
#    if len(glob('*ce.fits')) == 1:
#        events_list = glob('*ce.fits')[0]
#    elif len(glob('*ce.fits.gz')) == 1:
#        events_list = glob('*ce.fits.gz')[0] #events file
#    else:
#        events_list = ''
#    return events_list
 
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

# The how_many is only required for makecurves. 
how_many = 4 # number of objects to be auto-detected.

# The coordinates are only required for curves.
xp = 2000
yp = 2000

'''The following parameter determines how the background is determined.
If you prefer to manually specify a region where there are no stars, 
then give 'no' as the value. Otherwise, you can provide 'yes' and
background region will be automatically selected.'''
background_auto = 'yes'  # 'yes' or 'no'.

# If 'no', PLEASE FILL the following.
x_bg = 2100 # background X-coordinate.
y_bg = 1812 # background Y-coordinate.

# Following parameters need not be changed (unless you want to).
whole_figure_resolution = 256 # resolution of full figure.
sub_fig_size = 40 # size of sub figure.
fontsize = 9 #fontsize for plots
#######################################################################


def modify_string(events_list):
    if events_list[-5:] == '.fits':
        events_list = events_list[:-5]
    if events_list[-8:] == '.fits.gz':
        events_list = events_list[:-8]
    return events_list

def tobe_or_notobe(time, FrameCount, bwidth, framecount_per_sec):
    fc_time_start = time.min()
    fc_time_width = (FrameCount.max() - FrameCount.min()) / float(framecount_per_sec)
    fc_time_end = fc_time_width + fc_time_start
    nbin = (fc_time_end - fc_time_start) / bwidth
    return int(nbin) 

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

# To automatically estimate background region.
def auto_bg(fx, fy): 
    lowres_counts, lowres_xedges, \
    lowres_yedges, lowres_Image = plt.hist2d(fx, fy, 
                                             bins = 256, 
                                             norm = LogNorm())

    lowres_xcentres = (lowres_xedges[:-1] + lowres_xedges[1:]) / 2.
    lowres_ycentres = (lowres_yedges[:-1] + lowres_yedges[1:]) / 2.
    flat_counts = lowres_counts.flatten()
    x_mesh, y_mesh = np.meshgrid(lowres_ycentres, lowres_xcentres) #Notice the swap.
    x_mesh = x_mesh.flatten() 
    y_mesh = y_mesh.flatten()
    # To avoid the edges. 
    polished_array = [(lr_count, xm_p, ym_p) for lr_count, xm_p, ym_p
                      in zip(lowres_counts.flatten(), x_mesh, y_mesh)
                      if ((xm_p - 2400) ** 2 + (ym_p - 2400) ** 2) <= 1800 ** 2]

    sorted_counts = sorted(polished_array)
    five_percent = int(0.05 * len(sorted_counts))
    r_count, x_bg, y_bg = random.choice(sorted_counts[:five_percent])
    return lowres_Image, r_count, x_bg, y_bg   
 
def makecurves(events_list = events_list,
               radius = radius,
               how_many = how_many,
               bwidth = bwidth,
               framecount_per_sec = framecount_per_sec,
               background_auto = background_auto,
               sky_radius = sky_radius,
               x_bg = x_bg,
               y_bg = y_bg,
               whole_figure_resolution = whole_figure_resolution,
               sub_fig_size = sub_fig_size, 
               fontsize = fontsize):

    # Reading few columns.
    f = fits.open(events_list)
    time = f[1].data['MJD_L2']
    FrameCount = f[1].data['FrameCount']
    fx = f[1].data['Fx']
    fy = f[1].data['Fy']

    nbin_check = tobe_or_notobe(time = time, 
                                FrameCount = FrameCount, 
                                bwidth = bwidth,
                                framecount_per_sec = framecount_per_sec)

    if nbin_check < 1:
        print('\nThe events list contain little data OR check bwidth parameter.\n')
        return

    original_input = events_list
    path_to_events_list, events_list = ntpath.split(events_list)
    events_list = modify_string(events_list)

    # To find positions of interest (positions with maximum events).
    fxi = [int(round(s)) for s in fx]
    fyi = [int(round(s)) for s in fy]
    # Counting stuff to know who all are popular. 
    counter = Counter(zip(fxi, fyi))
    A = np.array(list(zip(*counter.most_common(how_many * 100)))[0])
    # Sieving out the duplicates. 
    uA = []
    while len(A) != 0:
        uA.append(A[0]) 
        A = np.array([x for x in A 
                      if x not in A[KDTree(A).query_ball_point(uA[-1], 15)]])

    uA = np.array(uA)[:how_many]

    if len(uA) == 0:
        print('No sources, try increasing the "how_many" parameter.')
        return

    # To avoid sources which have coordinates 0 or 4800.
    mask = np.isin(uA, [0, 4800], invert = True)
    mask = mask[:, 0] * mask[:, 1]
    uA = uA[mask]

    if len(uA) == 0:
        print('No sources, try increasing the "how_many" parameter.')
        return
    
    coo_file = os.path.join(path_to_events_list, 'sources_' + events_list +'.coo')
    np.savetxt(coo_file, uA, fmt = '%4.f\t%4.f')
    print('\nDetected source coordinates saved in file:\n* {}'.format(coo_file))

    # To automatically estimate background region.
    plt.figure(figsize = (10.5, 10))
    if background_auto == 'yes':
        lowres_Image, r_count, x_bg, y_bg = auto_bg(fx, fy)
        plt.gcf().gca().add_artist(lowres_Image)
       
    # To create a quick look figure marking sources and background.
    if whole_figure_resolution != 256 or background_auto == 'no':  # To avoid doing this twice.
        plt.hist2d(fx, fy, bins = whole_figure_resolution, 
                   norm = LogNorm())

    for u in uA:
        plt.annotate('Source', u, 
                     size = 13, color = 'black', fontweight = 'bold')

        obj_circle = plt.Circle(u, 100, color = 'k', fill = False)
        plt.gcf().gca().add_artist(obj_circle)

    plt.annotate('Background', (x_bg, y_bg),
                 size = 13, color = 'black', fontweight = 'bold')

    bg_circle = plt.Circle((x_bg, y_bg), 100, color = 'k', fill = False)
    plt.gcf().gca().add_artist(bg_circle)
    png_name = os.path.join(path_to_events_list, 'sources_' + events_list + '.png')
    plt.savefig( png_name, format = 'png', bbox_inches = 'tight')
    plt.clf()

    print('Detected sources are plotted in the image:\n* {}'.format(png_name))   
        
    bg_png = create_sub_image(x_bg, y_bg, 
                              sub_fig_size, 
                              sky_radius, 
                              'background_',
                              fx, fy,
                              path_to_events_list,
                              events_list)

    # For estimating background counts.
    F_b = [Fn_b for xx_b, yy_b, Fn_b in zip(fx, fy, FrameCount)
               if ((xx_b - x_bg)**2 + (yy_b - y_bg)**2) <= sky_radius**2]

    if len(F_b) != 0:
        scaled_events = (len(F_b) * radius**2) / float(sky_radius**2)
        scaled_events_e = (np.sqrt(len(F_b)) * radius**2) / float(sky_radius**2)
    else:
        scaled_events = 0
        scaled_events_e = 0

    # Converting FrameCounts array (the ones inside aperture) to time array.
    fc_time_start = time.min()
    fc_time_width = (FrameCount.max() - FrameCount.min()) / float(framecount_per_sec)
    fc_time_end = fc_time_width + fc_time_start
    # Converting FrameCounts array (note that np.unique is used) to time array.
    framecount_array = (np.unique(FrameCount) - FrameCount.min()) / float(framecount_per_sec)  
    framecount_time = framecount_array + fc_time_start

    # To estimate Background CPS.
    unique_FrameCount = np.unique(FrameCount)
    Number_of_frames = float(len(unique_FrameCount))
    bg_CPS = (scaled_events * framecount_per_sec) / Number_of_frames
    bg_CPS_e = (scaled_events_e * framecount_per_sec) / Number_of_frames
    print('\nThe estimated background CPS = {:.5f} +/-{:.5f}'.format(bg_CPS, bg_CPS_e))
    print('Region selected for background estimate:\n* {}'.format(bg_png))

    # Calculating number of bins.
    nbin = (fc_time_end - fc_time_start) / bwidth
    nbin = int(nbin)

    # End time of the bin range.
    till_here = fc_time_start + (bwidth * nbin)

    # Changing mission elapsed time in seconds to modified julian date. 
    jan2010 = 55197.0  # 2010.0(UTC) expressed with MJD format and scale UTC.
    till_here = (till_here / 86400.0) + jan2010
    fc_time_start = (fc_time_start / 86400.0) + jan2010
    framecount_time = (np.array(framecount_time) / 86400.0) + jan2010

    # Getting the number of unique frames within a bin.
    u_counts, u_bin_edges = np.histogram(framecount_time, bins = nbin,
                                         range = (fc_time_start, till_here),
                                         density = None)

    # selecting events within a circular region.
    print('\n---------------------- light curves ----------------------')
    plt.figure(figsize = (8, 5))
    for uaxy in uA:
        xp, yp = uaxy
        F = [Fn for xx, yy, Fn 
             in zip(fx, fy, FrameCount) 
             if ((xx - xp)**2 +(yy - yp)**2) <= radius**2]
        
        # Converting FrameCounts array to time array
        fc_time = (np.array(F) - FrameCount.min()) / float(framecount_per_sec)
        fc_time = fc_time + time.min()

        # Changing mission elapsed time in seconds to modified julian date. 
        fc_time = (np.array(fc_time) / 86400.0) + jan2010  # 1 julian day = 86400 seconds

        plt.title("X = %s, Y = %s, bin = %ss, radius = %spx" \
                  %(xp, yp, bwidth, radius), fontsize = fontsize)

        plt.xlabel("Time (Julian Date)", fontsize = fontsize)
        plt.ylabel("Counts per second", fontsize = fontsize)
        plt.tick_params(axis = 'both', labelsize = fontsize)
        counts,bin_edges = np.histogram(fc_time, bins = nbin, 
                                        range = (fc_time_start, till_here),
                                        density = None)    

        bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2.

        if np.array_equal(bin_edges, u_bin_edges) == True:
            mbmcuc = [(mb, mc, float(uc)) for mb, mc, uc 
                      in zip(bin_centres, counts, u_counts) if mc != 0]
        else:
            print('\nWhoa! What happened?\n')
            return

        if len(mbmcuc) != 0:
            mcentres, mcounts, frames_in_bin = zip(*mbmcuc)
        else:
            print('No counts for the source at %s' %uaxy)
            continue

        mcentres = np.array(mcentres)

        CPF = np.array(mcounts) / frames_in_bin
        CPF_err = np.sqrt(mcounts) / frames_in_bin
        CPS = CPF * framecount_per_sec
        CPS_err = CPF_err * framecount_per_sec

        #Background subtraction.
        CPS = CPS - bg_CPS    
        CPS_err = np.sqrt(CPS_err**2 + bg_CPS_e**2)

        plt.scatter(mcentres, CPS)
        plt.errorbar(mcentres, CPS, yerr = CPS_err, linestyle = "None")
        empty_space = (till_here - fc_time_start) / 25.0
        plt.xlim(fc_time_start - empty_space, till_here + empty_space)

#        median_time = np.median(framecount_time)
#        figname = 'makecurves_' + str(xp) + '_' + str(yp) + '_' + events_list + '_' + str(median_time) + ".png"
        figname = os.path.join(path_to_events_list, 'makecurves_' + str(xp) + '_' + str(yp) + '_' + events_list + '.png')
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
          background_auto = background_auto,
          sky_radius = sky_radius,
          x_bg = x_bg,
          y_bg = y_bg,
          whole_figure_resolution = whole_figure_resolution,
          sub_fig_size = sub_fig_size,
          fontsize = fontsize):

    # Reading few columns.
    f = fits.open(events_list)
    time = f[1].data['MJD_L2']
    FrameCount = f[1].data['FrameCount']
    fx = f[1].data['Fx']
    fy = f[1].data['Fy']

    nbin_check = tobe_or_notobe(time = time, 
                                FrameCount = FrameCount, 
                                bwidth = bwidth,
                                framecount_per_sec = framecount_per_sec)

    if nbin_check < 1:
        print('\nThe events list contain little data OR check bwidth parameter.\n')
        return

    path_to_events_list, events_list = ntpath.split(events_list)
    events_list = modify_string(events_list)

    # To automatically estimate background region.
    plt.figure(figsize = (10.5, 10))
    if background_auto == 'yes':
        lowres_Image, r_count, x_bg, y_bg = auto_bg(fx, fy)
        plt.gcf().gca().add_artist(lowres_Image)
        
    # To create a quick look figure marking source and background.
    if whole_figure_resolution != 256 or background_auto == 'no':  #Just to avoid doing this twice.
        plt.hist2d(fx, fy, bins = whole_figure_resolution, 
                   norm = LogNorm())

    plt.annotate("Source", (xp, yp), 
                 size = 13, color = 'black', fontweight = 'bold')

    obj_circle = plt.Circle((xp, yp), 100, color = 'k', fill = False)
    plt.gcf().gca().add_artist(obj_circle)

    plt.annotate("Background", (x_bg, y_bg),
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

    bg_png = create_sub_image(x_bg, y_bg, 
                              sub_fig_size, 
                              sky_radius, 
                              'background_',
                              fx, fy,
                              path_to_events_list,
                              events_list)

    # For estimating background counts.
    F_b = [Fn_b for xx_b, yy_b, Fn_b in zip(fx, fy, FrameCount)
               if ((xx_b - x_bg)**2 + (yy_b - y_bg)**2) <= sky_radius**2]

    if len(F_b) != 0:
        scaled_events = (len(F_b) * radius**2) / float(sky_radius**2)
        scaled_events_e = (np.sqrt(len(F_b)) * radius**2) / float(sky_radius**2)
    else:
        scaled_events = 0
        scaled_events_e = 0

    # selecting events within a circular region.
    F = [Fn for xx, yy, Fn 
            in zip(fx, fy, FrameCount)
            if ((xx - xp)**2 +(yy - yp)**2) <= radius**2]
            
    # Converting FrameCounts array (the ones inside aperture) to time array.
    fc_time_start = time.min()
    fc_time_array = (np.array(F) - FrameCount.min()) / float(framecount_per_sec)
    fc_time = fc_time_array + fc_time_start
    fc_time_width = (FrameCount.max() - FrameCount.min()) / float(framecount_per_sec)
    fc_time_end = fc_time_width + fc_time_start
    # Converting FrameCounts array (note that np.unique is used) to time array.
    framecount_array = (np.unique(FrameCount) - FrameCount.min()) / float(framecount_per_sec)  
    framecount_time = framecount_array + fc_time_start

    # To estimate Background CPS.
    unique_FrameCount = np.unique(FrameCount)
    Number_of_frames = float(len(unique_FrameCount))
    bg_CPS = (scaled_events * framecount_per_sec) / Number_of_frames
    bg_CPS_e = (scaled_events_e * framecount_per_sec) / Number_of_frames
    print("\nThe estimated background CPS = {:.5f} +/-{:.5f}".format(bg_CPS, bg_CPS_e))
    print('Region selected for background estimate:\n* {}'.format(bg_png))

    # Calculating number of bins.
    nbin = fc_time_width / bwidth
    nbin = int(nbin)

    # End time of the bin range.
    till_here = fc_time_start + (bwidth * nbin)

    # Changing mission elapsed time in seconds to modified julian date.
    jan2010 = 55197.0  # 2010.0(UTC) expressed with MJD format and scale UTC.
    fc_time = (np.array(fc_time) / 86400.0) + jan2010  # 1 julian day = 86400 seconds.
    till_here = (till_here / 86400.0) + jan2010
    fc_time_start = (fc_time_start / 86400.0) + jan2010
    framecount_time = (np.array(framecount_time) / 86400.0) + jan2010

    # Binning stuff, plotting stuff.
    plt.figure(figsize = (8, 5))
    plt.title('bin = %ss, radius = %spx' %(bwidth, radius), fontsize = fontsize)
    plt.xlabel("Time (Julian Date)", fontsize = fontsize)
    plt.ylabel("Counts per second", fontsize = fontsize)
    plt.tick_params(axis = 'both', labelsize = fontsize)
    counts, bin_edges = np.histogram(fc_time, bins = nbin,
                                     range = (fc_time_start, till_here),
                                     density = None)

    u_counts, u_bin_edges = np.histogram(framecount_time, bins = nbin,
                                         range = (fc_time_start, till_here),
                                         density = None)

    bin_centres = (bin_edges[:-1] + bin_edges[1:]) / 2.

    if np.array_equal(bin_edges, u_bin_edges) == True:
        mbmcuc = [(mb, mc, float(uc)) for mb, mc, uc
                  in zip(bin_centres, counts, u_counts) 
                  if mc != 0]
    else:
        print('\nWhoa! What happened?\n')
        return

    if len(mbmcuc) != 0:
        mcentres, mcounts, frames_in_bin = zip(*mbmcuc)
    else:
        print('No counts inside the aperture!')

    mcentres = np.array(mcentres)
    CPF = np.array(mcounts) / frames_in_bin
    CPF_err = np.sqrt(mcounts) / frames_in_bin
    CPS = CPF * framecount_per_sec
    CPS_err = CPF_err * framecount_per_sec

    # Background subtraction.
    CPS = CPS - bg_CPS
    CPS_err = np.sqrt(CPS_err**2 + bg_CPS_e**2)

    # Let us get on with the plot. 
    plt.scatter(mcentres, CPS)
    plt.errorbar(mcentres, CPS, yerr = CPS_err, linestyle = "None")
    # To make the plot look good.
    empty_space = (till_here - fc_time_start) / 25.0
    plt.xlim(fc_time_start - empty_space, till_here + empty_space)

    #To write the array to output.
    data_to_output = list(zip(mcentres, CPS, CPS_err))
    datname = os.path.join(path_to_events_list, 'curve_' + str(xp) + '_' + str(yp) + '_' + events_list + '.dat')
    np.savetxt(datname, data_to_output,
               fmt = '%10.11f\t%.5e\t%.5e',
               header = 'MJD\t\t\tCPS (bin=%ss)\tCPS_error' %bwidth)

    figname = os.path.join(path_to_events_list, 'curve_' + str(xp) + '_' + str(yp) + '_' + events_list + '.png')
    plt.savefig(figname, format = 'png', bbox_inches = 'tight', dpi = 150)

    print('\n-------------------------- curve --------------------------')
    print('source: {}\n        {}'.format(png_name, source_png))
    print('data: {}'.format(datname))
    print('plot: {}'.format(figname))

    print("\nDone!\n")
    plt.close('all') 

