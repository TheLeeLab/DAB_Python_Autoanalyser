#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 15:21:34 2023

Class related to analysis of DAB stained images

@author: jbeckwith
"""
import numpy as np
import skimage as ski
import cv2
from skimage.measure import label, regionprops_table
from skimage.color import (
    separate_stains,
    combine_stains,
    hdx_from_rgb,
    rgb_from_hdx,
    rgb2gray,
    hed2rgb,
    rgb2hed,
)
from skimage.morphology import reconstruction
from skimage.filters import threshold_otsu
import polars as pl
import matplotlib.pyplot as plt


class DAB:
    def __init__(self):
        self = self
        return

    def imread(self, file):
        """imread function
        takes RGB image and corrects openCV's ordering

        Args:
            file (str): file path

        Returns:
            img (np.2darray): array"""
        img = cv2.imread(file)
        img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
        return img

    def imdict_read(self, files):
        """imdict_read function
        takes list of files and makes image dict

        Args:
            files (array): array of file path

        Returns:
            img_dict (dict): dict of images"""

        img_dict = {}
        for file in files:
            img = cv2.imread(file)
            img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
            img_dict[file] = img
        return img_dict

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

    def clean_nuclear_mask(self, mask):
        """clean_nuclear_mask function
        takes image_mask, and cleans up nuclei
        removes <60 area (i.e. diffraction limited) objects
        clears border, connects larger aggregates if small "holes" inside, etc

        Args:
            image_mask (np.2darray): logical array of image mask

        Returns:
            cleaned_mask (np.2darray): cleaned up mask"""
        mask_disk = 1 * ski.morphology.binary_closing(mask, ski.morphology.disk(1))
        seed = np.copy(mask_disk)
        seed[1:-1, 1:-1] = mask_disk.max()

        mask_filled = ski.morphology.reconstruction(seed, mask_disk, method="erosion")
        cleaned_mask = ski.segmentation.clear_border(mask_filled)

        label_img = label(cleaned_mask)
        props = regionprops_table(label_img, properties=("area", "axis_minor_length"))
        Area = props["area"]
        indices_toremove = np.unique(np.unique(label_img)[1:] * (Area < 30))[1:]
        mask = np.isin(label_img, indices_toremove)
        cleaned_mask[mask] = 0
        return cleaned_mask

    def clean_protein_mask(self, image_mask):
        """clean_protein_mask function
        takes image_mask, and cleans up protein aggregates
        removes 3*3 (i.e. diffraction limited) objects
        clears border, connects larger aggregates if small "holes" inside, etc

        Args:
            image_mask (np.2darray): logical array of image mask

        Returns:
            cleaned_mask (np.2darray): cleaned up mask"""
        mask_disk = 1 * ski.morphology.binary_closing(
            image_mask, ski.morphology.disk(1)
        )
        seed = np.copy(mask_disk)
        seed[1:-1, 1:-1] = mask_disk.max()

        mask_filled = ski.morphology.reconstruction(seed, mask_disk, method="erosion")
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

    def otsu_filtering(self, image):
        """otsu threshold a single colour image

        Args:
            image (np.2darray): single grayscale image
        Returns:
            mask (np.2darray): single boolean image
            thresh (float): single boolean mask image"""

        seed = np.copy(image)
        seed[1:-1, 1:-1] = image.max()
        mask = image

        filled = reconstruction(seed, mask, method="erosion")
        holes = np.abs(image - filled)
        holes = holes / np.nanmax(holes)
        thresh = threshold_otsu(holes)
        mask = holes > thresh
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

    def otsu_filtering_multiimage(self, image_dict):
        """otsu threshold a single colour image

        Args:
            image_dict (dict): dict of images
        Returns:
            mask_dict (dict): dict of masks
            thresh (float): threshold value"""

        holes_dict = {}
        mask_dict = {}
        for i, key in enumerate(image_dict.keys()):
            image = image_dict[key]
            seed = np.copy(image)
            seed[1:-1, 1:-1] = image.max()
            mask = image

            filled = reconstruction(seed, mask, method="erosion")
            holes_dict[key] = np.abs(image - filled)
            if i == 0:
                holearray = holes_dict[key].ravel()
            else:
                holearray = np.hstack([holearray, holes_dict[key].ravel()])

        thresh = threshold_otsu(holearray)
        for i, key in enumerate(image_dict.keys()):
            mask_dict[key] = holes_dict[key] > thresh
        return mask_dict, thresh

    def analyse_DAB_multiimage(self, img_dict):
        """analyse_DAB function
        takes dictionary of images, and uses otsu's method to separate out
        coloured objects
        then returns table with object information

        Args:
            img (dict): dict of images image data; keys are filenames

        Returns:
            image_mask_asyn (dict): dict of masks
            table_asyn (pd.DataArray): pandas array of asyn data"""

        dab_image_dict = {}
        for key in img_dict.keys():
            ihc_hed = rgb2hed(self.im2double(img_dict[key]))
            # Create an RGB image for each of the stains
            null = np.zeros_like(ihc_hed[:, :, 0])
            ihc_d = hed2rgb(np.stack((null, null, ihc_hed[:, :, 2]), axis=-1))

            dab_image_dict[key] = rgb2gray(ihc_d)

        image_mask_asyn, thresh = self.otsu_filtering_multiimage(dab_image_dict)
        for i, key in enumerate(image_mask_asyn.keys()):
            image_mask_asyn[key] = self.clean_protein_mask(image_mask_asyn[key])

            label_img_asyn = label(image_mask_asyn[key])
            props_asyn = regionprops_table(
                label_img_asyn,
                properties=(
                    "area",
                    "centroid",
                    "axis_major_length",
                    "axis_minor_length",
                ),
            )

            table_asyn = pl.DataFrame(props_asyn)
            table_asyn = table_asyn.with_columns(
                pseudo_circularity=self.pseudo_circularity(
                    props_asyn["axis_major_length"], props_asyn["axis_minor_length"]
                )
            )
            table_asyn = table_asyn.with_columns(
                filename=np.full_like(
                    props_asyn["axis_minor_length"], str(key), dtype="object"
                )
            )
            yield key, image_mask_asyn[key], table_asyn

    def analyse_DAB(self, img, filename):
        """analyse_DAB function
        takes file, and uses initial parameters and rate to separate out
        coloured objects
        then returns table with object information

        Args:
            img (np.ndarray): image data
            filename (str): filename string
            asyn_params (np.1darry): initial default Lmean, aMean, bMean and threshold parameters
            nuclei_params (np.1darray): initial default Lmean, aMean, bMean and threshold parameters

        Returns:
            image_mask_asyn (np.2darray): mask of protein
            table_asyn (pd.DataArray): pandas array of asyn data"""

        ihc_hdx = ski.color.separate_stains(self.im2double(img), hdx_from_rgb)
        # Create an RGB image for each of the stains
        null = np.zeros_like(ihc_hdx[:, :, 0])
        ihc_d = combine_stains(
            np.stack((null, ihc_hdx[:, :, 1], null), axis=-1), rgb_from_hdx
        )

        dab_image = rgb2gray(ihc_d)

        image_mask_asyn_raw, thresh = self.otsu_filtering(dab_image)
        image_mask_asyn = self.clean_protein_mask(image_mask_asyn_raw)

        label_img_asyn = label(image_mask_asyn)
        props_asyn = regionprops_table(
            label_img_asyn,
            properties=("area", "centroid", "axis_major_length", "axis_minor_length"),
        )

        table_asyn = pl.DataFrame(props_asyn)
        table_asyn = table_asyn.with_columns(
            pseudo_circularity=self.pseudo_circularity(
                props_asyn["axis_major_length"], props_asyn["axis_minor_length"]
            )
        )
        table_asyn = table_asyn.with_columns(
            filename=np.full_like(
                props_asyn["axis_minor_length"], filename, dtype="object"
            )
        )

        return image_mask_asyn, table_asyn, thresh

    def analyse_DAB_and_cells(self, img, filename):
        """analyse_DAB function
        takes file, and uses initial parameters and rate to separate out
        coloured objects
        then returns table with object information

        Args:
            img (np.ndarray): image data
            filename (str): filename string

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

        dab_image = rgb2gray(ihc_d)

        image_mask_asyn_raw, thresh = self.otsu_filtering(dab_image)

        bin_blur = ski.filters.gaussian(image_mask_asyn_raw, sigma=2)
        bin_blur[bin_blur > 0] = 1
        bin_blur = np.asarray(bin_blur, dtype=bool)

        image_mask_asyn = self.clean_protein_mask(image_mask_asyn_raw)

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
        nuclei_image = rgb2gray(ihc_h)

        image_mask_nuclei_raw = self.binary_filtering(nuclei_image, thresh)
        image_mask_nuclei = self.clean_nuclear_mask(image_mask_nuclei_raw)

        label_img_asyn = label(image_mask_asyn)
        props_asyn = regionprops_table(
            label_img_asyn,
            properties=("area", "centroid", "axis_major_length", "axis_minor_length"),
        )

        table_asyn = pl.DataFrame(props_asyn)
        table_asyn = table_asyn.with_columns(
            pseudo_circularity=self.pseudo_circularity(
                props_asyn["axis_major_length"], props_asyn["axis_minor_length"]
            )
        )
        table_asyn = table_asyn.with_columns(
            filename=np.full_like(
                props_asyn["axis_minor_length"], filename, dtype="object"
            )
        )

        label_img_nucl = label(image_mask_nuclei)
        props_nuclei = regionprops_table(
            label_img_nucl,
            properties=("area", "centroid", "axis_major_length", "axis_minor_length"),
        )
        table_nuclei = pl.DataFrame(props_nuclei)
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
        )

    def plot_masks_and_histogram(
        self,
        img,
        masks,
        table_asyn,
        histcolor="gray",
        xaxislabel=r"object area (pixels$^2$)",
        alpha=1,
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
                axes[0].contour(masks[:, :, i], [0.5], linewidths=0.5, colors=colors[i])
        else:
            axes[0].contour(masks, [0.5], linewidths=0.5, colors="darkred")

        axes[0].axis("off")

        areas = table_asyn["area"].values
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
        return fig, axes

    def plot_masks(self, img, masks=None):
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
                    axes.contour(
                        masks[:, :, i], [0.5], linewidths=0.5, colors=colors[i]
                    )
            else:
                axes.contour(masks, [0.5], linewidths=0.5, colors="darkred")

            axes.axis("off")

        return fig, axes
