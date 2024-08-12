#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 18:20:06 2024

@author: jbeckwith
"""

from src import DAB_Analysis_Functions

DAB_A = DAB_Analysis_Functions.DAB()

overall_directory = "example_directory"
imtype = ".svs"
multi_tiff_pixel_size = 4000
pixel_size = 0.2528e-6
NA = 0.75
n_samples_forthresh = 50
percentile = 90
giveNN = False
save_figs = False

"""chop_and_analyse_svs function finds svs files,
chops them up, analyses them, and then deletes the chopped files

Args:
    overall_directory (string): directory to search for svs files in sub-directories
    multi_tiff_pixel_size (int): size of images to chop to
    pixel_size (float): pixel size in m
    NA (float): NA of objective used to generate data
    n_samples_forthresh (int): number of randomly chosen tiffs to use
                    to determine thresholds for overall sample set
    percentile (float): percentile of thresholds to use for overall data
    giveNN (boolean): If true, also calculates Nearest Neighbors.
                        Default False, as this is quite slow.
    save_figs (boolean): If True, saves figures showing the segmentation
                        for each image. Useful as a "sanity check"
"""

DAB_A.chop_and_analyse_svs(
    overall_directory,
    imtype=imtype,
    multi_tiff_pixel_size=multi_tiff_pixel_size,
    pixel_size=pixel_size,
    NA=NA,
    n_samples_forthresh=n_samples_forthresh,
    percentile=percentile,
    giveNN=giveNN,
    save_figs=save_figs,
)