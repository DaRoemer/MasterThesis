#-------------------------------------------------------
# Moprhology analysis tools
# By: Felix Romer
# Last updated: 26.04.2025
# Description:
# This script contains functions for analyzing and visualizing the morphology of segmented cells form microscopy images.
# It includes functions for counting neighboring regions, calculating actin fibers, and generating property images.
# The script also provides functions for processing and initializing images, performing region analysis, and saving results.
# The code is designed to work with labeled images and outlines, and it uses various libraries such as skimage, scipy, and matplotlib for image processing and visualization.
#-------------------------------------------------------
from scipy.ndimage import binary_dilation
import feret
from skimage.measure import profile_line
from skimage.filters import gaussian
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from skimage.measure import regionprops_table
from skimage import io, color
import tifffile as tiff
import os

#-------------------------------------------------------
# Analysis functions
#-------------------------------------------------------

def count_neighbors(region: dict, label_image: np.ndarray, labels: list) -> int:
    """
    Count the number of unique neighboring regions for a given region in a labeled image.

    Parameters:
        region (dict): The region properties dictionary.
        label_image (np.ndarray): The labeled image.
        labels (list): List of labels present in the image.

    Returns:
        int: Number of unique neighboring regions.
    """
    if region['label'] == 0:  # Skip the background label
        return None

    # Get the coordinates of the region (assuming coords is a list of tuples)
    coords = np.array(region['coords'])

    # Create a binary mask for the current region
    mask = np.zeros(label_image.shape, dtype=bool)
    mask[coords[:, 0], coords[:, 1]] = True

    # Dilate the mask by the specified number of iterations
    dilated_mask = binary_dilation(mask, iterations=1)

    # Extract the surrounding pixels from the original labeled image
    surrounding_pixels = label_image[dilated_mask]

    # Get unique neighbor labels and filter only those present in 'label' column
    neighbors = np.unique(surrounding_pixels)
    neighbors = [x for x in neighbors if x in labels]

    # Return the number of unique neighbors
    return len(neighbors)

def find_intersection(x1: float, y1: float, x2: float, y2: float, minf_angle: float, minf_t: float, x0: float, y0: float) -> list:
    """
    Find the intersection point between a line segment and a line defined by the minimum Feret angle and distance.

    Parameters:
        x1, y1, x2, y2 (float): Coordinates of the line segment.
        minf_angle (float): Minimum Feret angle.
        minf_t (float): Minimum Feret distance.
        x0, y0 (float): Coordinates of the point on the minimum Feret line.

    Returns:
        list: Coordinates of the intersection point.
    """
    # Line equation for the edge: y = m_edge * x + c_edge
    if x2 != x1:  # Avoid division by zero
        m_edge = (y2 - y1) / (x2 - x1)
        c_edge = y1 - m_edge * x1
    else:
        m_edge = None  # Vertical line
        c_edge = x1

    # Line equation for the min feret line: y = m_line * x + c_line
    m_line = np.tan(minf_angle)
    c_line = minf_t

    # If minf_angle is 90 degrees, the min feret line is vertical
    if minf_angle == np.pi / 2:
        x_intersect = x0
        y_intersect = m_edge * x_intersect + c_edge
        return [y_intersect, x_intersect]
    
    # If minf_angle is 0 degrees, the min feret line is horizontal
    if minf_angle == 0:
        y_intersect = y0
        if m_edge is None:
            x_intersect = x1
        else:
            x_intersect = (y_intersect - c_edge) / (m_line - m_edge)
        return [y_intersect, x_intersect]

    

    # Find intersection
    if m_edge is not None:
        # Solve for x: m_line * x + c_line = m_edge * x + c_edge
        x_intersect = (c_edge - c_line) / (m_line - m_edge)
        y_intersect = m_line * x_intersect + c_line
    else:
        # Vertical line case: x = c_edge
        y_intersect = c_edge
        x_intersect = m_line * x_intersect + c_line

    # Check if the intersection is within the segment bounds
    #if min(x1, x2) <= x_intersect <= max(x1, x2) and min(y1, y2) <= y_intersect <= max(y1, y2):
    return[y_intersect, x_intersect]



def calculate_actin_fibers(region_label: int, actin_image: np.ndarray, label_image: np.ndarray, fcal: float) -> tuple:
    """
    Calculate the number of actin fibers in a given region of an actin image.

    Parameters:
        region_label (int): The label of the region.
        actin_image (np.ndarray): The actin image.
        label_image (np.ndarray): The labeled image.
        fcal (float): Calibration factor.

    Returns:
        tuple: 
            - int: Number of actin fibers.
            - np.ndarray: Profile line.
            - np.ndarray: Smoothed line.
            - list: Peaks.
            - list: Start and end points.
            - float: Minimum Feret distance.
    """
    im = np.where(label_image == region_label, 1, 0)
    y_max, x_max = im.shape
    feret_results = feret.calc(im, edge = True)
    if feret_results is None:
        return None, None, None, None, None
    
    coords = feret_results.minf_coords

    
    x1, y1 = coords.T[1][1], coords.T[0][1]
    x2, y2 = coords.T[1][0], coords.T[0][0]
    y3, x3 = coords.T[1][2], coords.T[0][2]


    intersections = find_intersection(x1, y1, x2, y2, feret_results.minf_angle, feret_results.minf_t, y3, x3)

    start = [x3, y3]
    end = intersections
    #print(start, end)
    profile            = profile_line(actin_image, start,end, linewidth=1)
    smoothed_line      = gaussian(profile, sigma=1)
    pks, _             = find_peaks(smoothed_line, prominence=50)
    num_actin_fibers   = len(pks)
    min_feret_distance = feret_results.minf * fcal
    return num_actin_fibers, profile, smoothed_line, pks, [start, end], min_feret_distance

#-------------------------------------------------------
# Data processing and visulaization
#-------------------------------------------------------
# Processing and initialization
def process_and_initialize(file: str, actin_file: str, image_dir: str, actin_dir: str, outdir: str, flow_direction_map: dict, flow_direction: str, params: dict, npy: bool = True) -> tuple:
    """
    Process and initialize the labeled image and outlines from a file.

    Parameters:
        file (str): The file name.
        actin_file (str): The actin file name.
        image_dir (str): The directory containing the image.
        actin_dir (str): The directory containing the actin image.
        outdir (str): The output directory.
        flow_direction_map (dict): Mapping of flow directions.
        flow_direction (str): The flow direction.
        params (dict): Analysis parameters.
        npy (bool): Whether the file is a numpy file.

    Returns:
        tuple: 
            - np.ndarray: Labeled image.
            - np.ndarray: Outlines.
            - np.ndarray: Actin image.
            - str: Name of the file.
            - str: Filename for saving.
            - str: Output directory.
    """
    name_region, display_orientation, crop_image, crop_heigth, actin = map(params.get, 
                                                                    ['name_region', 
                                                                     'display_orientation', 
                                                                     'crop_image', 
                                                                     'crop_heigth',
                                                                     'calculate_actin_fibers_number'])    
    name = file[:-8]  # Remove the '_seg.npy' part
    dirname = image_dir + file
    filename = outdir + name + name_region

    # Read and rotate image
    if npy:
        label_image = np.load(dirname, allow_pickle=True).item().get('masks')
        outlines    = np.load(dirname, allow_pickle=True).item().get('outlines')
    else:
        label_image = io.imread(dirname)
        label_image = color.rgb2lab(label_image)
        outlines    = np.zeros(label_image.shape)
    # Load actin image
    if actin == True:
        actin_image = io.imread(image_dir + actin_file)
    else:
        actin_image = np.zeros(label_image.shape)
    rotations = 0
    if display_orientation != flow_direction:
        rotations = (flow_direction_map[display_orientation] - flow_direction_map[flow_direction]) % 4
        label_image = np.rot90(label_image, k=rotations)
        outlines = np.rot90(outlines, k=rotations)
        actin_image = np.rot90(actin_image, k=rotations)

    # Crop image if needed
    if crop_image:
        y_min, y_max = np.where(label_image != 0)[0].min(), np.where(label_image != 0)[0].max()
        if y_max - y_min >= crop_heigth:
            crop_min = y_min
            crop_max = y_max
        else:
            crop_min = y_min - (crop_heigth - (y_max - y_min)) // 2
            crop_max = crop_min + crop_heigth
        label_image = label_image[crop_min:crop_max, :]
        outlines = outlines[crop_min:crop_max, :]
        actin_image = actin_image[crop_min:crop_max, :]
    
    # make outdir
    outdir = outdir + name + '/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    return label_image, outlines, actin_image, name, filename, outdir

# Region analysis
def region_analysis(label_image: np.ndarray, params: dict) -> pd.DataFrame:
    """
    Analyze regions in a labeled image and filter them based on area constraints.

    Parameters:
        label_image (np.ndarray): The labeled image.
        params (dict): Analysis parameters.

    Returns:
        pd.DataFrame: DataFrame containing the filtered region properties.
    """
    min_area, max_area, fcal = map(params.get, ['min_area', 'max_area', 'fcal', 'calculate_actin_fibers_number'])
    stats = pd.DataFrame(regionprops_table(label_image,
                                            properties=('label', 'area', 'orientation', 'major_axis_length', 'minor_axis_length', 'perimeter', 'coords', 'centroid')))
    
    stats_filtered = stats[(stats['area'] > min_area) & (stats['area'] < max_area)].reset_index(drop=True)

    # add rows 
    stats_filtered['aspect_ratio']       = stats_filtered['major_axis_length'] / stats_filtered['minor_axis_length']
    stats_filtered['area_um']            = stats_filtered['area'] * fcal**2
    stats_filtered['tortuosity']         = stats_filtered['perimeter'] / np.sqrt(stats_filtered['area'])
    stats_filtered['orientation_degree'] = 90-np.degrees(stats_filtered['orientation']).abs()

    # get number of nieghbors
    lables = stats_filtered['label'].values
    stats_filtered['num_neighbors']      = stats_filtered.apply(count_neighbors, axis=1, label_image=label_image, labels=lables)

    # Read the actin image 
    stats_filtered['num_actin_fibers']    = 0
    stats_filtered['actin_fibers_per_um'] = 0
    stats_filtered['intersections']       = 0
    return stats_filtered

# Visualize
def create_propertie_image(stats_filtered: pd.DataFrame, label_image: np.ndarray, outlines: np.ndarray, current_outdir: str, name: str, analysis_params: dict) -> list:
    """
    Create property images for visualization of region properties.

    Parameters:
        stats_filtered (pd.DataFrame): Filtered region properties.
        label_image (np.ndarray): The labeled image.
        outlines (np.ndarray): Outlines of the regions.
        current_outdir (str): Current output directory.
        name (str): Name of the file.
        analysis_params (dict): Analysis parameters.

    Returns:
        list: List of property images.
    """
    # Extract cmaps & normalizations
    cmap_area, cmap_orientation, cmap_asp_ratio, cmap_tortuosity, cmap_num_neighbors, cmap_actin_fibers = analysis_params.get('cmaps')
    norm_area, norm_orientation, norm_asp_ratio, norm_tortuosity, norm_num_neighbors, norm_actin_fibers = analysis_params.get('normalizations')
    save_plots = analysis_params.get('save_plots')
    name_orientation = analysis_params.get('name_orientation')
    name_area = analysis_params.get('name_area')
    name_asp_ratio = analysis_params.get('name_asp_ratio')
    name_tortuosity = analysis_params.get('name_tortuosity')
    name_num_neighbors = analysis_params.get('name_num_neighbors')
    name_actin_fibers = analysis_params.get('name_actin_fibers')
    # Initialize empty images for visualization
    im_orientation   = np.zeros((label_image.shape[0], label_image.shape[1], 3), dtype=np.uint8)
    im_area          = np.zeros((label_image.shape[0], label_image.shape[1], 3), dtype=np.uint8)
    im_asp_ratio     = np.zeros((label_image.shape[0], label_image.shape[1], 3), dtype=np.uint8)
    im_tortuosity    = np.zeros((label_image.shape[0], label_image.shape[1], 3), dtype=np.uint8)
    im_num_neighbors = np.zeros((label_image.shape[0], label_image.shape[1], 3), dtype=np.uint8)
    im_actin_fibers  = np.zeros((label_image.shape[0], label_image.shape[1], 3), dtype=np.uint8)
    for ind, region in stats_filtered.iterrows():
        if region['label'] == 0:  # Skip the background
            continue
        
        area_um           = region['area_um']
        orientation       = region['orientation_degree']
        aspect_ratio      = region['aspect_ratio']
        tortuosity        = region['tortuosity']
        num_neighbors     = region['num_neighbors']
        num_actin_fibers  = region['actin_fibers_per_um']

        # Normalize the area value for color mapping
        rgb_area          = (np.array(cmap_area         (norm_area         (area_um))         [:3]) * 255).astype(np.uint8)
        rgb_orientation   = (np.array(cmap_orientation  (norm_orientation  (orientation))     [:3]) * 255).astype(np.uint8)
        rgb_asp_ratio     = (np.array(cmap_asp_ratio    (norm_asp_ratio    (aspect_ratio))    [:3]) * 255).astype(np.uint8)
        rgb_tortuosity    = (np.array(cmap_tortuosity   (norm_tortuosity   (tortuosity))      [:3]) * 255).astype(np.uint8)
        rgb_num_neighbors = (np.array(cmap_num_neighbors(norm_num_neighbors(num_neighbors))   [:3]) * 255).astype(np.uint8)
        rgb_actin_fibers  = (np.array(cmap_actin_fibers (norm_actin_fibers (num_actin_fibers))[:3]) * 255).astype(np.uint8)
    
        # Apply the color to all pixels in the region
        for coordinates in region['coords']:
            im_area[         coordinates[0], coordinates[1]] = rgb_area
            im_orientation[  coordinates[0], coordinates[1]] = rgb_orientation
            im_asp_ratio[    coordinates[0], coordinates[1]] = rgb_asp_ratio
            im_tortuosity[   coordinates[0], coordinates[1]] = rgb_tortuosity
            im_num_neighbors[coordinates[0], coordinates[1]] = rgb_num_neighbors
            im_actin_fibers[ coordinates[0], coordinates[1]] = rgb_actin_fibers
        
    outlines_expanded = np.expand_dims(outlines, axis=-1)

    im_area         = np.where(outlines_expanded != 0, 0, im_area)
    im_orientation  = np.where(outlines_expanded != 0, 0, im_orientation)
    im_asp_ratio    = np.where(outlines_expanded != 0, 0, im_asp_ratio)
    im_tortuosity   = np.where(outlines_expanded != 0, 0, im_tortuosity)
    im_num_neighbors= np.where(outlines_expanded != 0, 0, im_num_neighbors)
    im_actin_fibers = np.where(outlines_expanded != 0, 0, im_actin_fibers)

    # Save images to disk
    if save_plots:
        tiff.imwrite(current_outdir + name + name_orientation, im_orientation)
        tiff.imwrite(current_outdir + name + name_area, im_area)
        tiff.imwrite(current_outdir + name + name_asp_ratio, im_asp_ratio)
        tiff.imwrite(current_outdir + name + name_tortuosity, im_tortuosity)
        tiff.imwrite(current_outdir + name + name_num_neighbors, im_num_neighbors)
        tiff.imwrite(current_outdir + name + name_actin_fibers, im_actin_fibers)
    return [im_orientation, im_area, im_asp_ratio, im_tortuosity, im_num_neighbors, im_actin_fibers]

def save_as_figure(prop_images: list, stats_filtered: pd.DataFrame, current_outdir: str, name: str, analysis_params: dict):
    """
    Save property images as figures with additional annotations.

    Parameters:
        prop_images (list): List of property images.
        stats_filtered (pd.DataFrame): Filtered region properties.
        current_outdir (str): Current output directory.
        name (str): Name of the file.
        analysis_params (dict): Analysis parameters.
    """
    # Extract images
    im_orientation, im_area, im_asp_ratio, im_tortuosity, im_num_neighbors, im_actin_fibers = prop_images
    cmap_area, cmap_orientation, cmap_asp_ratio, cmap_tortuosity, cmap_num_neighbors, cmap_actin_fibers = analysis_params.get('cmaps')
    norm_area, norm_orientation, norm_asp_ratio, norm_tortuosity, norm_num_neighbors, norm_actin_fibers = analysis_params.get('normalizations')
    name_region_prop = analysis_params.get('name_region_prop')
    name_region_prop_numbers = analysis_params.get('name_region_prop_numbers')
    # Plot and save the visualization
    fig, axes = plt.subplots(2, 3, figsize=(20, 5))

    # Orientation plot
    im0 = axes[0, 0].imshow(im_orientation, cmap=cmap_orientation, norm=norm_orientation)
    axes[0, 0].set_title('Orientation')
    axes[0, 0].axis('off')  # Remove axis ticks and labels
    cb0 = fig.colorbar(im0, ax=axes[0, 0], orientation='horizontal')
    cb0.set_label('Angle (degrees)')

    # Area plot
    im1 = axes[0, 1].imshow(im_area, cmap=cmap_area, norm=norm_area)
    axes[0, 1].set_title('Area')
    axes[0, 1].axis('off')  # Remove axis ticks and labels
    cb1 = fig.colorbar(im1, ax=axes[0, 1], orientation='horizontal')
    cb1.set_label('Area (\u03bcm\u00b2)')  # Unicode for µm²

    # Aspect Ratio plot
    im2 = axes[0, 2].imshow(im_asp_ratio, cmap=cmap_asp_ratio, norm=norm_asp_ratio)
    axes[0, 2].set_title('Aspect Ratio')
    axes[0, 2].axis('off')  # Remove axis ticks and labels
    cb2 = fig.colorbar(im2, ax=axes[0, 2], orientation='horizontal')
    cb2.set_label('Aspect Ratio')

    # Tortuosity plot
    im3 = axes[1, 0].imshow(im_tortuosity, cmap=cmap_tortuosity, norm=norm_tortuosity)
    axes[1, 0].set_title('Shape Factor')
    axes[1, 0].axis('off')  # Remove axis ticks and labels
    cb3 = fig.colorbar(im3, ax=axes[1, 0], orientation='horizontal')
    cb3.set_label('Shape Factor')

    # Number of Neighbors plot
    im4 = axes[1, 1].imshow(im_num_neighbors, cmap=cmap_num_neighbors, norm=norm_num_neighbors)
    axes[1, 1].set_title('Number of Neighbors')
    axes[1, 1].axis('off')  # Remove axis ticks and labels
    cb4 = fig.colorbar(im4, ax=axes[1, 1], orientation='horizontal')
    cb4.set_label('Number of Neighbors')
    

    # Actin Fibers plot
    im5 = axes[1, 2].imshow(im_actin_fibers, cmap=cmap_actin_fibers, norm=norm_actin_fibers)
    axes[1, 2].set_title('Actin Fibers')
    axes[1, 2].axis('off')  # Remove axis ticks and labels
    cb5 = fig.colorbar(im5, ax=axes[1, 2], orientation='horizontal')
    cb5.set_label('Number of Actin Fibers per \u03bcm')  # Unicode for µm
    

    plt.tight_layout()


    plot_filename = current_outdir + name + name_region_prop        
    fig.savefig(plot_filename)

    # Plot the number of neighbors of each region at the position of their centroid
    for ind, region in stats_filtered.iterrows():
        axes[1, 1].text(region['centroid-1'], region['centroid-0'], str(region['num_neighbors']),
                        color='white', fontsize=8, ha='center', va='center')
        
    # Plot the number of actin fibers of each region at the position of their centroid
    for ind, region in stats_filtered.iterrows():
        axes[1, 2].text(region['centroid-1'], region['centroid-0'], str(region['actin_fibers_per_um']),
                        color='white', fontsize=8, ha='center', va='center')

    plot_filename_numbers = current_outdir + name + name_region_prop_numbers
    fig.savefig(plot_filename_numbers)

