import numpy as np
import matplotlib.pyplot as plt
import slideio
import re 
from shapely.geometry import Polygon
from PIL import Image
import matplotlib.path as mpl_path

slice_name = '771791'
slide = slideio.open_slide('data_SH/'+slice_name+'.svs','SVS')
num_scenes = slide.num_scenes
scene = slide.get_scene(0)
print('pixel numbers: ',scene.rect)
raw_string = slide.raw_metadata
print('slide metadata: ',raw_string)
### Set the multi-tiff size in pixel
multi_tiff_pixel_size = 4000

def get_identical_element_indices(arr):
    indices = []
    for i in range(len(arr)):
        for j in range(i+1, len(arr)):
            if arr[i] == arr[j]:
                if (i+1)!=j:
                    indices.append([i, j])
    return indices
# pixel_size = 0.2528
with open('data_SH/'+slice_name+'.svs_stroma.txt', 'r') as f:
    # Read the entire file as a string
    file_str = f.read()       

# Use regular expressions to extract the numerical values from each line
pattern = r'Point:\s*(-?\d+\.\d+)\s*,\s*(-?\d+\.\d+)'
matches = re.findall(pattern, file_str)

# Convert the matches to a list of tuples
polygon_vertices = np.array([[float(x), float(y)] for x, y in matches]) ## unit: pixel
exterior = [(x, y) for x, y in polygon_vertices]
identical_element_indices = get_identical_element_indices(exterior)
if len(identical_element_indices) == 0:
    polygon_vertices = np.vstack((polygon_vertices[:],polygon_vertices[0]))
    exterior = [(x, y) for x, y in polygon_vertices]
    identical_element_indices = get_identical_element_indices(exterior)
polygons = []

for index in identical_element_indices:
    sub_points = polygon_vertices[index[0]:(index[1]+1)]
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
print('Number of interior polygons:', len(interior_polys))
print('Number of exterior polygons:', len(exterior_polys))
for polygon in exterior_polys:
    coords_array = np.array(polygon.exterior.coords)
    min_x = int(np.min(coords_array[:, 0]))
    min_y = int(np.min(coords_array[:, 1]))
    max_x = int(np.max(coords_array[:, 0]))
    max_y = int(np.max(coords_array[:, 1]))
    for x_cor in range(min_x, max_x, multi_tiff_pixel_size):
        end_x_cor = min(x_cor + multi_tiff_pixel_size, max_x)
        for y_cor in range(min_y, max_y, multi_tiff_pixel_size):
            end_y_cor = min(y_cor + multi_tiff_pixel_size, max_y)
            print('x_cor: ', x_cor, 'end_x_cor: ', end_x_cor, 'y_cor: ', y_cor, 'end_y_cor: ', end_y_cor)
            img = scene.read_block((x_cor, y_cor, multi_tiff_pixel_size, multi_tiff_pixel_size))
            print('img shape: ', img.shape)
            tiff_path = f'image{x_cor}_{y_cor}.tiff'  # Update this with your desired output path
            img_pil = Image.fromarray(img)
            img_pil.save(tiff_path, format='TIFF')
            ### create corrdinate for each multi-tiff
            x_cor_list = np.linspace(x_cor, end_x_cor-1, multi_tiff_pixel_size)
            y_cor_list = np.linspace(y_cor, end_y_cor-1, multi_tiff_pixel_size)
            # Create a grid of coordinates using meshgrid and reshape the result
            x_grid, y_grid = np.meshgrid(x_cor_list, y_cor_list)
            corr = np.vstack([x_grid.ravel(), y_grid.ravel()]).T
            mask_total = np.zeros_like(corr[:,0], dtype=bool)
            for polygon in exterior_polys:
                polygon_path = mpl_path.Path(polygon.exterior.coords)
                mask = polygon_path.contains_points(corr)
                mask_total = np.logical_or(mask_total,mask)
            corr = corr[mask_total]
            mask_total = np.reshape(mask_total, (multi_tiff_pixel_size, multi_tiff_pixel_size))
            mask_total = np.repeat(mask_total[:, :, np.newaxis], 3, axis=2)
            new_image = mask_total*img
            new_image_path = f'masked_image{x_cor}_{y_cor}.tiff'
            img_pil = Image.fromarray(new_image)
            img_pil.save(new_image_path, format='TIFF')