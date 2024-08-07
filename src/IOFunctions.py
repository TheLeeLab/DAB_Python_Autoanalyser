# -*- coding: utf-8 -*-
"""
This class contains functions pertaining to IO of files based for the RASP code.
jsb92, 2024/01/02
"""
import os
import numpy as np
import polars as pl

import slideio
import re
from shapely.geometry import Polygon
from PIL import Image
import matplotlib.path as mpl_path
import time
from copy import copy

import sys


module_dir = os.path.dirname(__file__)
sys.path.append(module_dir)


class IO_Functions:
    def __init__(self):
        self = self
        tiff = Image.open("src/.example.ome.tif")
        self.example_tags = tiff.tag_v2
        return

    def chop_svs_file(
        self,
        overall_directory,
        slice_name,
        multi_tiff_pixel_size=4000,
        pixel_size=0.2528e-6,
    ):
        slide = slideio.open_slide(
            os.path.join(overall_directory, slice_name + ".svs"), "SVS"
        )
        scene = slide.get_scene(0)

        temp_directory = os.path.join(overall_directory, "temp")
        if not os.path.isdir(temp_directory):
            os.mkdir(temp_directory)

        def get_identical_element_indices(arr):
            indices = []
            for i in range(len(arr)):
                for j in range(i + 1, len(arr)):
                    if arr[i] == arr[j]:
                        if (i + 1) != j:
                            indices.append([i, j])
            return indices

        with open(
            os.path.join(overall_directory, slice_name + ".svs_stroma.txt"), "r"
        ) as f:
            # Read the entire file as a string
            file_str = f.read()

        # Use regular expressions to extract the numerical values from each line
        pattern = r"Point:\s*(-?\d+\.\d+)\s*,\s*(-?\d+\.\d+)"
        matches = re.findall(pattern, file_str)

        # Convert the matches to a list of tuples
        polygon_vertices = np.array(
            [[float(x), float(y)] for x, y in matches]
        )  ## unit: pixel
        exterior = [(x, y) for x, y in polygon_vertices]
        identical_element_indices = get_identical_element_indices(exterior)
        if len(identical_element_indices) == 0:
            polygon_vertices = np.vstack((polygon_vertices[:], polygon_vertices[0]))
            exterior = [(x, y) for x, y in polygon_vertices]
            identical_element_indices = get_identical_element_indices(exterior)
        polygons = []

        for index in identical_element_indices:
            sub_points = polygon_vertices[index[0] : (index[1] + 1)]
            polygon = Polygon([(x, y) for x, y in sub_points])
            polygons.append(polygon)
        interior_polys = []
        exterior_polys = []

        # Check each polygon against the others
        for i, poly1 in enumerate(polygons):
            is_interior = False
            for j, poly2 in enumerate(polygons):
                if i != j and poly1.within(poly2):
                    is_interior = True
                    break
            if is_interior:
                interior_polys.append(poly1)
            else:
                exterior_polys.append(poly1)

        start = time.time()

        polygon_paths = {}
        for i, polygon in enumerate(exterior_polys):
            polygon_paths[i] = mpl_path.Path(polygon.exterior.coords)

        j = 0
        k = 0
        l = 0
        maxval = 0

        for polygon in exterior_polys:
            coords_array = np.array(polygon.exterior.coords)
            min_x = int(np.min(coords_array[:, 0]))
            min_y = int(np.min(coords_array[:, 1]))
            max_x = int(np.max(coords_array[:, 0]))
            max_y = int(np.max(coords_array[:, 1]))
            maxval += len(np.arange(min_x, max_x, multi_tiff_pixel_size)) * len(
                np.arange(min_y, max_y, multi_tiff_pixel_size)
            )

        for polygon in exterior_polys:

            coords_array = np.array(polygon.exterior.coords)
            min_x = int(np.min(coords_array[:, 0]))
            min_y = int(np.min(coords_array[:, 1]))
            max_x = int(np.max(coords_array[:, 0]))
            max_y = int(np.max(coords_array[:, 1]))
            xrange = np.arange(min_x, max_x, multi_tiff_pixel_size)
            yrange = np.arange(min_y, max_y, multi_tiff_pixel_size)

            for x_cor in xrange:
                end_x_cor = min(x_cor + multi_tiff_pixel_size, max_x)
                for y_cor in yrange:
                    end_y_cor = min(y_cor + multi_tiff_pixel_size, max_y)
                    img = scene.read_block(
                        (x_cor, y_cor, multi_tiff_pixel_size, multi_tiff_pixel_size)
                    )

                    ### create corrdinate for each multi-tiff
                    x_cor_list = np.linspace(
                        x_cor, end_x_cor - 1, multi_tiff_pixel_size
                    )
                    y_cor_list = np.linspace(
                        y_cor, end_y_cor - 1, multi_tiff_pixel_size
                    )

                    # Create a grid of coordinates using meshgrid and reshape the result
                    x_grid, y_grid = np.meshgrid(x_cor_list, y_cor_list)
                    corr = np.vstack([x_grid.ravel(), y_grid.ravel()]).T
                    mask_total = np.zeros_like(corr[:, 0], dtype=bool)
                    for polygon in exterior_polys:
                        polygon_path = mpl_path.Path(polygon.exterior.coords)
                        mask = polygon_path.contains_points(corr)
                        mask_total = np.logical_or(mask_total, mask)
                    mask_total = np.reshape(
                        mask_total, (multi_tiff_pixel_size, multi_tiff_pixel_size)
                    )
                    if np.nansum(mask_total.ravel()) > 0:
                        new_image_path = os.path.join(
                            temp_directory, f"image_{x_cor}_{y_cor}.tiff"
                        )
                        img_pil = Image.fromarray(img, mode="RGB")
                        tags_tosave = copy(self.example_tags)
                        tags_tosave[296] = 3
                        tags_tosave[282] = 1e-3 / pixel_size
                        tags_tosave[283] = 1e-3 / pixel_size
                        tags_tosave[270] = "Python=3.10.12f\nunit=micron\n"
                        tags_tosave[256] = img_pil.width
                        tags_tosave[257] = img_pil.height
                        tags_tosave[278] = img_pil.height
                        tags_tosave[279] = img_pil.height * img_pil.width * 3
                        img_pil.save(
                            new_image_path, format="TIFF", tiffinfo=tags_tosave
                        )

                        mask_image_path = os.path.join(
                            temp_directory, f"image_mask_{x_cor}_{y_cor}.tiff"
                        )
                        mask_pil = Image.fromarray(mask_total)
                        masktags_tosave = copy(self.example_tags)
                        masktags_tosave[296] = 3
                        del masktags_tosave[258]
                        masktags_tosave[259] = 1
                        masktags_tosave[262] = 1
                        del masktags_tosave[273]
                        del masktags_tosave[50838]
                        del masktags_tosave[50839]
                        masktags_tosave[282] = 1e-3 / pixel_size
                        masktags_tosave[283] = 1e-3 / pixel_size
                        masktags_tosave[270] = "Python=3.10.12f\nunit=micron\n"
                        masktags_tosave[277] = 1
                        del masktags_tosave[256]
                        del masktags_tosave[257]
                        del masktags_tosave[278]
                        del masktags_tosave[279]
                        mask_pil.save(
                            mask_image_path, format="TIFF", tiffinfo=masktags_tosave
                        )
                    l += 1
                    print(
                        "Iterating over sliced DAB   subsection {}/{}   Time elapsed: {:.3f} s".format(
                            j + k + l,
                            maxval,
                            time.time() - start,
                        ),
                        end="\r",
                        flush=True,
                    )
                k += 1
            j += 1
        return
