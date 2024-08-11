#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 15:21:34 2023

Class related to analysis of DAB stained images

@author: jbeckwith
"""
import numpy as np
import skimage as ski
import fnmatch
from skimage.measure import label, regionprops_table
from skimage.color import (
    separate_stains,
    combine_stains,
    hdx_from_rgb,
    rgb_from_hdx,
    rgb2gray,
)
import shutil
from skimage.morphology import reconstruction
from skimage.filters import threshold_yen
import polars as pl
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.spatial import cKDTree
import sys
import os
import json
import time
import pathos
from pathos.pools import ThreadPool as Pool

module_dir = os.path.dirname(__file__)
sys.path.append(module_dir)
import IOFunctions

IO = IOFunctions.IO_Functions()
cpu_number = int(pathos.helpers.cpu_count() * 0.8)


class DAB:
    def __init__(self):
        self = self
        return

    def file_search(self, folder, string1):
        """
        Search for files with 'string1' as their end string within 'folder'
        Args:
            folder (str): The directory to search for files.
            string1 (str): The first string to search for in the filenames.

        Returns:
            file_list (list): A sorted list of file paths matching the search criteria.
        """
        # Get a list of all files containing 'string1' in their names within 'folder'
        file_list = [
            os.path.join(dirpath, f)
            for dirpath, dirnames, files in os.walk(folder)
            for f in fnmatch.filter(files, "*" + string1)
        ]
        return file_list

    def chop_and_analyse_svs(
        self,
        overall_directory,
        imtype=".svs",
        multi_tiff_pixel_size=4000,
        pixel_size=0.2528e-6,
        NA=0.75,
        n_samples_forthresh=50,
        percentile=90,
        giveNN=False,
        save_figs=False,
    ):
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
        svs_files = self.file_search(
            overall_directory, imtype
        )  # first get all files in any subfolders
        ws = " "

        for file in svs_files:
            slice_name = os.path.split(file.split(".svs")[0])[-1]
            IO.chop_svs_file(
                os.path.split(file)[0],
                slice_name,
                multi_tiff_pixel_size=multi_tiff_pixel_size,
                pixel_size=pixel_size,
            )
            print(
                ws * 80,
                end="\r",
                flush=True,
            )

            folder = os.path.split(file)[0]
            temp_folder = os.path.join(folder, "temp")
            
            if save_figs == True:
                fig_folder = os.path.join(folder, "figures")
                if not os.path.isdir(fig_folder):
                    os.mkdir(fig_folder)

            tiffs = self.file_search(temp_folder, ".tiff")
            tiffs = np.sort([e for e in tiffs if "mask" not in e])

            thresh_asyn_variance = np.zeros(n_samples_forthresh)
            thresh_nuclei_variance = np.zeros(n_samples_forthresh)
            if n_samples_forthresh < len(tiffs):
                random_tiffs = np.random.choice(
                    tiffs, size=n_samples_forthresh, replace=False
                )
            else:
                random_tiffs = tiffs
                n_samples_forthresh = len(tiffs)

            start = time.time()

            for i in np.arange(n_samples_forthresh):
                data_tothresh = IO.imread(random_tiffs[i])
                maskname = os.path.split(random_tiffs[i])[-1]
                maskname = maskname.split("_")
                maskname.insert(1, "mask")
                maskname = "_".join(maskname)
                maskname = os.path.join(temp_folder, maskname)
                mask = IO.imread(maskname)
                (
                    thresh_asyn_variance[i],
                    thresh_nuclei_variance[i],
                ) = self.analyse_DAB_and_cells(
                    data_tothresh,
                    mask=mask,
                    filename=os.path.split(random_tiffs[i])[-1],
                    give_NN=False,
                    justthresh=True,
                )
                print(
                    "Thresholding; tiff {}/{}   Time elapsed: {:.2f} minutes".format(
                        i + 1,
                        n_samples_forthresh,
                        (time.time() - start) / 60,
                    ),
                    end="\r",
                    flush=True,
                )
            print(
                ws * 80,
                end="\r",
                flush=True,
            )
            thresh_asyn = np.percentile(thresh_asyn_variance, percentile)
            thresh_nuclei = np.percentile(thresh_nuclei_variance, percentile)

            to_save = {"thresh_protein": thresh_asyn, "thresh_nuclei": thresh_nuclei}
            with open(
                os.path.join(folder, slice_name + "_thresholds.json"),
                "w",
                encoding="utf-8",
            ) as f:
                json.dump(to_save, f, ensure_ascii=False, indent=4)

            pixel_size_um = pixel_size * 1e6

            start = time.time()

            for i in np.arange(len(tiffs)):
                data_toanalyse = IO.imread(tiffs[i])
                xpos = np.multiply(
                    pixel_size_um, int(os.path.split(tiffs[i])[-1].split("_")[1])
                )
                ypos = np.multiply(
                    pixel_size_um,
                    int(os.path.split(tiffs[i])[-1].split("_")[-1].split(".tiff")[0]),
                )
                maskname = os.path.split(tiffs[i])[-1]
                maskname = maskname.split("_")
                maskname.insert(1, "mask")
                maskname = "_".join(maskname)
                maskname = os.path.join(temp_folder, maskname)
                mask = IO.imread(maskname)
                (
                    image_mask_asyn,
                    table_asyn_temp,
                    image_mask_nuclei,
                    table_nuclei_temp,
                    thresh,
                    thresh_nuclei,
                ) = self.analyse_DAB_and_cells(
                    data_toanalyse,
                    mask=mask,
                    filename=os.path.split(tiffs[i])[-1],
                    give_NN=giveNN,
                    thresh=thresh_asyn,
                    thresh_nuclei=thresh_nuclei,
                    justthresh=False,
                    pixel_size=pixel_size_um,
                )

                xc = pl.Series(
                    "centroid-x", table_asyn_temp["centroid-0"].to_numpy() + xpos
                )
                yc = pl.Series(
                    "centroid-y", table_asyn_temp["centroid-1"].to_numpy() + ypos
                )
                table_asyn_temp = table_asyn_temp.replace_column(0, xc)
                table_asyn_temp = table_asyn_temp.replace_column(1, yc)

                xnc = pl.Series(
                    "centroid-x", table_nuclei_temp["centroid-0"].to_numpy() + xpos
                )
                ync = pl.Series(
                    "centroid-y", table_nuclei_temp["centroid-1"].to_numpy() + ypos
                )
                table_nuclei_temp = table_nuclei_temp.replace_column(0, xnc)
                table_nuclei_temp = table_nuclei_temp.replace_column(1, ync)

                if save_figs == True:
                    import matplotlib

                    matplotlib.use("Agg")
                    fig, axs = self.plot_masks(
                        IO.imread(tiffs[i]),
                        np.dstack([image_mask_asyn, image_mask_nuclei]),
                    )
                    figname = (
                        str(int(os.path.split(tiffs[i])[-1].split("_")[1]))
                        + "_"
                        + str(
                            int(
                                os.path.split(tiffs[i])[-1]
                                .split("_")[-1]
                                .split(".tiff")[0]
                            )
                        )
                        + "_segmentation_example.svg"
                    )
                    plt.savefig(
                        os.path.join(fig_folder, figname), format="svg", dpi=600
                    )

                if i == 0:
                    table_asyn = table_asyn_temp
                    table_nuclei = table_nuclei_temp
                else:
                    table_asyn = pl.concat(
                        [table_asyn, table_nuclei_temp], rechunk=True
                    )
                    table_nuclei = pl.concat(
                        [table_nuclei, table_nuclei_temp], rechunk=True
                    )

                print(
                    "Analysing tiffs; tiff {}/{}   Time elapsed: {:.2f} minutes".format(
                        i + 1,
                        len(tiffs),
                        (time.time() - start) / 60,
                    ),
                    end="\r",
                    flush=True,
                )
            print(
                ws * 80,
                end="\r",
                flush=True,
            )
            savename_asyn = os.path.join(
                folder,
                slice_name
                + "_pixel_size_"
                + str(pixel_size).replace(".", "p")
                + "_NA_"
                + str(NA).replace(".", "p")
                + "_protein_stain_analysis.csv",
            )
            savename_nuclei = os.path.join(
                folder,
                slice_name
                + "_pixel_size_"
                + str(pixel_size).replace(".", "p")
                + "_NA_"
                + str(NA).replace(".", "p")
                + "_nuclei_stain_analysis.csv",
            )
            table_asyn.write_csv(savename_asyn)
            table_nuclei.write_csv(savename_nuclei)

            shutil.rmtree(temp_folder)
        return

    def im2double(self, img):
        """im2double function
        takes image and normalises to double

        Args:
            img (np.2darray): image object

        Returns:
            imdouble (np.2darray): numpy array converted to double"""
        info = np.iinfo(img.dtype)  # Get the data type of the input image
        imdouble = (
            img.astype(np.float32) / info.max
        )  # Divide all values by the largest possible value in the datatype
        return imdouble

    def bincalculator(self, data):
        """bincalculator function
        reads in data and generates bins according to Freedman-Diaconis rule

        Args:
            data (np.1darray): data to calculate bins

        Returns:
            bins (np.1darray): bins for histogram according to Freedman-Diaconis rule"""
        N = len(data)
        sigma = np.std(data)

        binwidth = np.multiply(np.multiply(np.power(N, np.divide(-1, 3)), sigma), 3.5)
        bins = np.linspace(
            np.min(data),
            np.max(data),
            int((np.max(data) - np.min(data)) / binwidth) + 1,
        )
        return bins

    def pseudo_circularity(self, MajorAxisLength, MinorAxisLength):
        """pseudo_circularity function
        takes major and minor axis length and computes pseudo-circularity

        Args:
            MajorAxisLength (np.float): major axis length in pixels
            MinorAxisLength (np.float): minor axis length in pixels

        Returns:
            p_circ (np.float64): pseudo-circularity (runs from 0--1)"""
        p_circ = np.divide(
            np.multiply(2, MinorAxisLength), np.add(MinorAxisLength, MajorAxisLength)
        )
        return p_circ

    def d_calculation(self, wavelength=0.6, NA=0.75, pixel_size=0.2528):
        """d_calculation function
        takes wavelength, NA and pixel size and calculates diffraction limit

        Args:
            wavelength (float): wavelength in same units as pixel_size
            NA (float): numerical aperture
            pixel_size (float): size of pixel

        Returns:
            cleaned_mask (np.2darray): cleaned up mask"""

        d = int(np.around((0.6 / (2 * NA)) / pixel_size))
        if d < 1:
            d = 1
        return d

    def clean_nuclear_mask(self, image_mask, mask=[], pixel_size=0.2528, NA=0.75):
        """clean_nuclear_mask function
        takes image_mask, and cleans up nuclei
        removes <60 area (i.e. diffraction limited) objects
        clears border, connects larger aggregates if small "holes" inside, etc

        Args:
            image_mask (np.2darray): logical array of image mask

        Returns:
            cleaned_mask (np.2darray): cleaned up mask"""
        d = self.d_calculation(NA=NA, pixel_size=pixel_size)

        mask_disk = 1 * ski.morphology.binary_closing(
            image_mask, ski.morphology.disk(d)
        )
        seed = np.copy(mask_disk)
        seed[1:-1, 1:-1] = mask_disk.max()

        mask_filled = ski.morphology.reconstruction(seed, mask_disk, method="erosion")

        if len(mask) > 0:
            cleaned_mask = mask * ski.segmentation.clear_border(mask_filled)
        else:
            cleaned_mask = ski.segmentation.clear_border(mask_filled)

        label_img = label(cleaned_mask)
        props = regionprops_table(label_img, properties=("area", "axis_minor_length"))
        Area = props["area"]
        indices_toremove = np.unique(np.unique(label_img)[1:] * (Area < 30))[1:]
        mask = np.isin(label_img, indices_toremove)
        cleaned_mask[mask] = 0
        return cleaned_mask

    def clean_protein_mask(self, image_mask, mask=[], pixel_size=0.2528, NA=0.75):
        """clean_protein_mask function
        takes image_mask, and cleans up protein aggregates
        removes 3*3 (i.e. diffraction limited) objects
        clears border, connects larger aggregates if small "holes" inside, etc

        Args:
            image_mask (np.2darray): logical array of image mask

        Returns:
            cleaned_mask (np.2darray): cleaned up mask"""
        d = self.d_calculation(NA=NA, pixel_size=pixel_size)

        mask_disk = 1 * ski.morphology.binary_closing(
            image_mask, ski.morphology.disk(d)
        )
        seed = np.copy(mask_disk)
        seed[1:-1, 1:-1] = mask_disk.max()

        mask_filled = ski.morphology.reconstruction(seed, mask_disk, method="erosion")
        cleaned_mask = ski.segmentation.clear_border(mask_filled)

        if len(mask) > 0:
            cleaned_mask = mask * ski.segmentation.clear_border(mask_filled)
        else:
            cleaned_mask = ski.segmentation.clear_border(mask_filled)

        label_img = label(cleaned_mask)
        props = regionprops_table(label_img, properties=("area", "axis_minor_length"))
        minorA = props["axis_minor_length"]
        Area = props["area"]
        indices_toremove = np.unique(
            np.hstack(
                [
                    np.unique(label_img)[1:] * (minorA < 3),
                    np.unique(label_img)[1:] * (Area < 9),
                ]
            )
        )[1:]
        image_mask = np.isin(label_img, indices_toremove)
        cleaned_mask[image_mask] = 0
        return cleaned_mask

    # TODO: speed up
    def neighbourdistance_self(
        self, table_coords, NN_level=1, pixel_size=0.2528, NA=0.75
    ):
        """N nearest neighbours of single table of coordinates

        Args:
            table_coords (pl.Series): polars series of table coordinates
            NN_level (int): number of neighbours to count
            pixel_size (float): pixel size in um

        Returns:
            NN_array (np.ndarray): NN_level by len(table_coords) array"""

        NN_array = np.zeros([len(table_coords), NN_level])
        for i in np.arange(len(table_coords)):
            temp_NN = np.zeros(len(table_coords) - (i + 1))
            jrange = np.arange(len(table_coords))
            jrange = jrange[jrange != i]

            def distancecalc(kj):
                k = kj[0]
                j = kj[1]
                min_dists, min_dists_idx = cKDTree(table_coords[int(j)]).query(
                    table_coords[int(i)], NN_level
                )
                temp_NN[k] = np.min(min_dists)

            pool = Pool(nodes=cpu_number)
            pool.restart()
            pool.map(distancecalc, zip(np.arange(len(jrange)), jrange))
            pool.close()
            pool.terminate()

            NN_array[i, :] = np.sort(temp_NN)[:NN_level]
        return np.multiply(pixel_size, NN_array)

    # TODO: speed up
    def neighbourdistance_AtoB(
        self, table_coords_A, table_coords_B, NN_level=1, pixel_size=0.2528, NA=0.75
    ):
        """N nearest neighbours between two table of coordinates

        Args:
            table_coords_A (pl.Series): polars series of table coordinates, 1st array
            table_coords_B (pl.Series): polars series of table coordinates, 2nd array
            NN_level (int): number of neighbours to count
            pixel_size (float): pixel size in um

        Returns:
            NN_array (np.ndarray): NN_level by len(table_coords_A) array"""

        NN_array = np.zeros([len(table_coords_A), NN_level])
        for i in np.arange(len(table_coords_A)):
            temp_NN = np.zeros(len(table_coords_B))
            jrange = np.arange(len(table_coords_B))

            def distancecalc(kj):
                k = kj[0]
                j = kj[1]
                min_dists, min_dists_idx = cKDTree(table_coords_B[int(j)]).query(
                    table_coords_A[int(i)], NN_level
                )
                temp_NN[k] = np.min(min_dists)

            pool = Pool(nodes=cpu_number)
            pool.restart()
            pool.map(distancecalc, zip(np.arange(len(jrange)), jrange))
            pool.close()
            pool.terminate()

            NN_array[i, :] = np.sort(temp_NN)[:NN_level]
        return np.multiply(pixel_size, NN_array)

    def yen_filtering(self, image):
        """yen threshold a single colour image

        Args:
            image (np.2darray): single grayscale image
        Returns:
            mask (np.2darray): single boolean image
            thresh (float): single boolean mask image"""

        thresh = threshold_yen(image)
        mask = image <= thresh
        return mask, thresh

    def binary_filtering(self, image, threshold):
        """threshold a single colour image

        Args:
            image (np.2darray): single grayscale image
            threshold (float): specified threshold
        Returns:
            mask (np.2darray): single boolean image"""

        seed = np.copy(image)
        seed[1:-1, 1:-1] = image.max()
        mask = image

        filled = reconstruction(seed, mask, method="erosion")
        holes = np.abs(image - filled)
        holes = holes / np.nanmax(holes)
        mask = holes > threshold
        return mask

    def correct_props_pixel_size(self, props, pixel_size=0.2528):
        """correct property matrix with pixel size

        Args:
            props (pl.Series): properties
            pixel_size (float): pixel size in um
        Returns:
            props (pl.Series): properties"""

        props["area"] = np.multiply(props["area"], np.square(pixel_size))
        props["centroid-0"] = np.multiply(props["centroid-0"], pixel_size)
        props["centroid-0"] = np.multiply(props["centroid-0"], pixel_size)
        props["axis_major_length"] = np.multiply(props["axis_major_length"], pixel_size)
        props["axis_minor_length"] = np.multiply(props["axis_minor_length"], pixel_size)
        return props

    def NN_single_loop(self, label_img_object, table_object, pixel_size, NN_level):
        """NN_single_loop collects loop of nearest neighbour computation

        Args:
            label_img_asyn (np.2darray): labelled asyn image
            table_object (pl.Series): table of asyn info
            pixel_size (float): pixel size in um
            NN_level (int): nearest neighbour level
        Returns:
            table_object (pl.Series): table of asyn info"""

        coords_A = regionprops_table(label_img_object, properties=["coords"])
        table_coords = coords_A["coords"]
        if len(table_coords) > 0:
            A_NN = self.neighbourdistance_self(
                table_coords, NN_level=NN_level, pixel_size=pixel_size
            )
        if len(table_coords) > 0:
            if NN_level == 1:
                table_object = table_object.with_columns(
                    NN_tosametype_distance=np.squeeze(A_NN)
                )
            else:
                for i in np.arange(NN_level, dtype=int):
                    table_object = table_object.with_columns(
                        NN_tosametype_distance=np.squeeze(A_NN[:, i])
                    ).rename(
                        {
                            "NN_tosametype_distance": "NN_tosametype_distance_level"
                            + (i + 1)
                        }
                    )
        return table_object

    def NN_multi_loop(
        self,
        label_img_object_A,
        table_object_A,
        label_img_object_B,
        pixel_size,
        NN_level,
    ):
        """NN_multi_loop collects loop of nearest neighbour computation
        from A to B

        Args:
            label_img_asyn (np.2darray): labelled asyn image
            table_object (pl.Series): table of asyn info
            pixel_size (float): pixel size in um
            NN_level (int): nearest neighbour level
        Returns:
            table_object (pl.Series): table of asyn info"""

        coords_A = regionprops_table(label_img_object_A, properties=["coords"])
        table_coords_A = coords_A["coords"]
        coords_B = regionprops_table(label_img_object_B, properties=["coords"])
        table_coords_B = coords_B["coords"]

        if (len(table_coords_A) > 0) & (len(table_coords_B) > 0):
            A_NN = self.neighbourdistance_AtoB(
                table_coords_A, table_coords_B, NN_level=NN_level, pixel_size=pixel_size
            )

        if (len(table_coords_A) > 0) & (len(table_coords_B) > 0):
            if NN_level == 1:
                table_object_A = table_object_A.with_columns(
                    NN_todifferenttype_distance=np.squeeze(A_NN)
                )
            else:
                for i in np.arange(NN_level, dtype=int):
                    table_object_A = table_object_A.with_columns(
                        NN_todifferenttype_distance=np.squeeze(A_NN[:, i])
                    ).rename(
                        {
                            "NN_todifferenttype_distance": "NN_todifferenttype_distance_level"
                            + (i + 1)
                        }
                    )
        elif (len(table_coords_A) > 0) & (len(table_coords_B) == 0):
            if NN_level == 1:
                table_object_A = table_object_A.with_columns(
                    NN_todifferenttype_distance=np.full(len(table_coords_A), np.NAN)
                )
            else:
                for i in np.arange(NN_level, dtype=int):
                    table_object_A = table_object_A.with_columns(
                        NN_todifferenttype_distance=np.full(len(table_coords_A), np.NAN)
                    ).rename(
                        {
                            "NN_todifferenttype_distance": "NN_todifferenttype_distance_level"
                            + (i + 1)
                        }
                    )
        return table_object_A

    def analyse_DAB(
        self,
        img,
        filename,
        pixel_size=0.2528,
        NA=0.75,
        mask=[],
        give_NN=True,
        NN_level=1,
        thresh=None,
        justthresh=False,
    ):
        """analyse_DAB function
        takes file, and uses initial parameters and rate to separate out
        coloured objects
        then returns table with object information

        Args:
            img (np.ndarray): image data
            filename (str): filename string
            pixel_size (float): pixel size in um
            mask (np.2darray): image mask to ignore some areas in analysis
            give_NN (bool): if True, will calculate nearest neighbours of objects
            NN_level (int): how many neighbours to calculate

        Returns:
            image_mask_asyn (np.2darray): mask of protein
            table_asyn (pd.DataArray): pandas array of asyn data
            thresh (float): float used in thresholding"""

        ihc_hdx = ski.color.separate_stains(self.im2double(img), hdx_from_rgb)
        # Create an RGB image for each of the stains
        null = np.zeros_like(ihc_hdx[:, :, 0])
        ihc_d = combine_stains(
            np.stack((null, ihc_hdx[:, :, 1], null), axis=-1), rgb_from_hdx
        )

        if len(mask) > 0:
            dab_image = ma.masked_array(rgb2gray(ihc_d), 1 - mask)
        else:
            dab_image = rgb2gray(ihc_d)

        if isinstance(thresh, type(None)):
            image_mask_asyn_raw, thresh = self.yen_filtering(dab_image)
            if justthresh == True:
                return thresh
        else:
            image_mask_asyn_raw = dab_image <= thresh

        if len(mask) > 0:
            image_mask_asyn = self.clean_protein_mask(
                image_mask_asyn_raw, mask=mask, pixel_size=pixel_size, NA=NA
            )
        else:
            image_mask_asyn = self.clean_protein_mask(
                image_mask_asyn_raw, pixel_size=pixel_size, NA=NA
            )

        label_img_asyn = label(image_mask_asyn)
        props_asyn = regionprops_table(
            label_img_asyn,
            properties=("area", "centroid", "axis_major_length", "axis_minor_length"),
        )

        # update to be real areas and locations etc
        props_asyn = self.correct_props_pixel_size(props_asyn, pixel_size=pixel_size)

        table_asyn = pl.DataFrame(props_asyn)
        table_asyn = table_asyn.with_columns(
            pseudo_circularity=self.pseudo_circularity(
                props_asyn["axis_major_length"], props_asyn["axis_minor_length"]
            )
        )

        if give_NN == True:
            table_asyn = self.NN_single_loop(
                label_img_asyn, table_asyn, pixel_size, NN_level
            )

        table_asyn = table_asyn.with_columns(
            filename=np.full_like(
                props_asyn["axis_minor_length"], filename, dtype="object"
            )
        )

        return image_mask_asyn, table_asyn, thresh

    def analyse_DAB_and_cells(
        self,
        img,
        filename,
        pixel_size=0.2528,
        NA=0.75,
        mask=[],
        give_NN=True,
        NN_level=1,
        thresh=None,
        thresh_nuclei=None,
        justthresh=False,
    ):
        """analyse_DAB function
        takes file, and uses initial parameters and rate to separate out
        coloured objects
        then returns table with object information

        Args:
            img (np.ndarray): image data
            filename (str): filename string
            pixel_size (float): pixel size in um
            mask (np.2darray): image mask to ignore some areas in analysis
            give_NN (bool): if True, will calculate nearest neighbours of objects
            NN_level (int): how many neighbours to calculate

        Returns:
            image_mask_asyn (np.2darray): mask of protein
            table_asyn (pd.DataArray): pandas array of asyn data
            image_mask_nuclei (np.2darray): mask of nuclei
            table_nuclei (pd.DataArray): pandas array of nuclei data"""

        ihc_hdx = ski.color.separate_stains(self.im2double(img), hdx_from_rgb)
        # Create an RGB image for each of the stains
        null = np.zeros_like(ihc_hdx[:, :, 0])
        ihc_d = combine_stains(
            np.stack((null, ihc_hdx[:, :, 1], null), axis=-1), rgb_from_hdx
        )

        if len(mask) > 0:
            dab_image = ma.masked_array(rgb2gray(ihc_d), 1 - mask)
        else:
            dab_image = rgb2gray(ihc_d)

        if isinstance(thresh, type(None)):
            image_mask_asyn_raw, thresh = self.yen_filtering(dab_image)
        else:
            image_mask_asyn_raw = dab_image <= thresh

        bin_blur = ski.filters.gaussian(image_mask_asyn_raw, sigma=2)
        bin_blur[bin_blur > 0] = 1
        bin_blur = np.asarray(bin_blur, dtype=bool)

        if len(mask) > 0:
            image_mask_asyn = self.clean_protein_mask(
                image_mask_asyn_raw, mask, pixel_size=pixel_size, NA=NA
            )
        else:
            image_mask_asyn = self.clean_protein_mask(
                image_mask_asyn_raw, pixel_size=pixel_size, NA=NA
            )

        data_1 = img[:, :, 0]
        data_2 = img[:, :, 1]
        data_3 = img[:, :, 2]

        data_1[bin_blur] = np.median(data_1)
        data_2[bin_blur] = np.median(data_2)
        data_3[bin_blur] = np.median(data_3)
        data_redone = np.dstack([data_1, data_2, data_3])
        data_sep = separate_stains(self.im2double(data_redone), hdx_from_rgb)

        null = np.zeros_like(data_sep[:, :, 0])
        ihc_h = combine_stains(
            np.stack((data_sep[:, :, 0], null, null), axis=-1), rgb_from_hdx
        )

        if len(mask) > 0:
            nuclei_image = ma.masked_array(rgb2gray(ihc_h), 1 - mask)
        else:
            nuclei_image = rgb2gray(ihc_h)

        if isinstance(thresh_nuclei, type(None)):
            image_mask_nuclei_raw, thresh_nuclei = self.yen_filtering(nuclei_image)
            if justthresh == True:
                return thresh, thresh_nuclei
        else:
            image_mask_nuclei_raw = nuclei_image <= thresh_nuclei

        if len(mask) > 0:
            image_mask_nuclei = self.clean_nuclear_mask(
                image_mask_nuclei_raw, mask, pixel_size=pixel_size, NA=NA
            )
        else:
            image_mask_nuclei = self.clean_nuclear_mask(
                image_mask_nuclei_raw, pixel_size=pixel_size, NA=NA
            )

        label_img_asyn = label(image_mask_asyn)
        props_asyn = regionprops_table(
            label_img_asyn,
            properties=("area", "centroid", "axis_major_length", "axis_minor_length"),
        )

        # update to be real areas and locations etc
        props_asyn = self.correct_props_pixel_size(props_asyn, pixel_size=pixel_size)

        table_asyn = pl.DataFrame(props_asyn)
        table_asyn = table_asyn.with_columns(
            pseudo_circularity=self.pseudo_circularity(
                props_asyn["axis_major_length"], props_asyn["axis_minor_length"]
            )
        )

        if give_NN == True:
            table_asyn = self.NN_single_loop(
                label_img_asyn, table_asyn, pixel_size, NN_level
            )

        label_img_nucl = label(image_mask_nuclei)
        props_nuclei = regionprops_table(
            label_img_nucl,
            properties=("area", "centroid", "axis_major_length", "axis_minor_length"),
        )

        # update to be real areas and locations etc
        props_nuclei = self.correct_props_pixel_size(
            props_nuclei, pixel_size=pixel_size
        )

        table_nuclei = pl.DataFrame(props_nuclei)

        if give_NN == True:
            table_nuclei = self.NN_single_loop(
                label_img_nucl, table_nuclei, pixel_size, NN_level
            )
            table_asyn = self.NN_multi_loop(
                label_img_asyn, table_asyn, label_img_nucl, pixel_size, NN_level
            )
            table_nuclei = self.NN_multi_loop(
                label_img_nucl, table_nuclei, label_img_asyn, pixel_size, NN_level
            )

        table_asyn = table_asyn.with_columns(
            filename=np.full_like(
                props_asyn["axis_minor_length"], filename, dtype="object"
            )
        )

        table_nuclei = table_nuclei.with_columns(
            pseudo_circularity=self.pseudo_circularity(
                props_nuclei["axis_major_length"], props_nuclei["axis_minor_length"]
            )
        )
        return (
            image_mask_asyn,
            table_asyn,
            image_mask_nuclei,
            table_nuclei,
            thresh,
            thresh_nuclei,
        )

    def plot_masks_and_protein_histogram(
        self,
        img,
        masks,
        table_asyn,
        histcolor="gray",
        xaxislabel=r"object area (pixels$^2$)",
        alpha=1,
        pixel_size=0.2528,
        lw=0.25,
    ):
        """plot_masks function
        takes image, and optional masks, and plots them together

        Args:
            img (np.ndarray): image data
            masks (np.ndarry): mask data
            table_asyn (pd.DataFrame): object data
        Returns:
            fig (object): figure object
            axes (object): axis object"""

        fig, axes = plt.subplots(1, 2)

        axes[0].imshow(img)
        if len(masks.shape) > 2:  # if multiple masks
            colors = ["darkred", "darkblue"]
            for i in np.arange(masks.shape[2]):  # plot multiple masks
                axes[0].contour(masks[:, :, i], [0.5], linewidths=lw, colors=colors[i])
        else:
            axes[0].contour(masks, [0.5], linewidths=lw, colors="darkred")

        axes[0].axis("off")

        areas = table_asyn["area"].to_numpy() * np.square(pixel_size)
        bins = self.bincalculator(areas)

        axes[1].hist(
            areas,
            bins=bins,
            density=False,
            color=histcolor,
            alpha=alpha,
        )
        axes[1].grid(True, which="both", ls="--", c="gray", lw=0.25, alpha=0.25)
        axes[1].set_ylabel("frequency", fontsize=8)
        axes[1].set_xlabel(r"area ($\mu$m$^2$)", fontsize=8)
        return fig, axes

    def plot_masks_and_both_histograms(
        self,
        img,
        masks,
        table_asyn,
        table_nuclei,
        xaxislabel=r"object area (pixels$^2$)",
        alpha=0.5,
        pixel_size=0.2528,
        lw=0.25,
    ):
        """plot_masks function
        takes image, and optional masks, and plots them together

        Args:
            img (np.ndarray): image data
            masks (np.ndarry): mask data
            table_asyn (pd.DataFrame): object data
        Returns:
            fig (object): figure object
            axes (object): axis object"""

        fig, axes = plt.subplots(1, 2)

        axes[0].imshow(img)
        if len(masks.shape) > 2:  # if multiple masks
            colors = ["darkred", "darkblue"]
            for i in np.arange(masks.shape[2]):  # plot multiple masks
                axes[0].contour(masks[:, :, i], [0.5], linewidths=lw, colors=colors[i])
        else:
            axes[0].contour(masks, [0.5], linewidths=lw, colors="darkred")

        axes[0].axis("off")

        areas = table_asyn["area"].to_numpy() * np.square(pixel_size)
        bins = self.bincalculator(areas)

        axes[1].hist(
            areas,
            bins=bins,
            density=False,
            color="darkred",
            alpha=alpha,
        )

        areas = table_nuclei["area"].to_numpy() * np.square(pixel_size)
        bins = self.bincalculator(areas)

        axes[1].hist(
            areas,
            bins=bins,
            density=False,
            color="blue",
            alpha=alpha,
        )

        axes[1].grid(True, which="both", ls="--", c="gray", lw=0.25, alpha=0.25)
        axes[1].set_ylabel("frequency", fontsize=8)
        axes[1].set_xlabel(r"area ($\mu$m$^2$)", fontsize=8)
        return fig, axes

    def plot_masks(self, img, masks=None, lw=0.25):
        """plot_masks function
        takes image, and optional masks, and plots them together

        Args:
            img (np.ndarray): image data
            masks (np.ndarry): mask data

        Returns:
            fig (object): figure object
            axes (object): axis object"""

        if isinstance(masks, type(None)):
            fig, axes = plt.subplots(1, 1)
            axes.imshow(img)
            axes.axis("off")
        else:
            fig, axes = plt.subplots(1, 1)

            axes.imshow(img)
            if len(masks.shape) > 2:  # if multiple masks
                colors = ["darkred", "darkblue"]
                for i in np.arange(masks.shape[2]):  # plot multiple masks
                    axes.contour(masks[:, :, i], [0.5], linewidths=lw, colors=colors[i])
            else:
                axes.contour(masks, [0.5], linewidths=lw, colors="darkred")

            axes.axis("off")

        return fig, axes
