#-------------------------------------------------------
# AFT_tools.py
# By Marcotti et al. 2021
# DOI: 10.3389/fcomp.2021.745831
# Adapted by: Felix Romer
# Last update: 26.04.2025
# Description:
# This script contains functions for analyzing and visualizing the local alignment of images.
# It includes functions for periodic decomposition, least moment calculation, local order calculation,
# image normalization, and parameter search for window size and neighborhood radius.
# The script also includes functions for comparing parameters between two datasets and plotting the results.
#-------------------------------------------------------

#-------------------------------------------------------
# Import necessary libraries
#-------------------------------------------------------
from typing import List, Tuple, Union, Optional
import numpy as np
from scipy.fft import fft2, fftshift, ifft2, ifftshift, fft, ifft, fftfreq
import skimage.io as io
import cv2                                                     # for filtering vector fields
from skimage.morphology import disk        # morphology operations
import numpy.matlib as matlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu
import os

#-------------------------------------------------------
# Utility functions
#-------------------------------------------------------

def image_norm(im: np.ndarray) -> np.ndarray:
    """
    Calculate the norm of an image.

    Args:
        im (np.ndarray): Input image.

    Returns:
        np.ndarray: Normalized image.
    """
    im_norm = np.sqrt(np.real(im * np.conj(im)))
    return im_norm

def pad_to_square(image: np.ndarray) -> np.ndarray:
    """
    Pad an image to make it square.

    Args:
        image (np.ndarray): Input image.

    Returns:
        np.ndarray: Padded square image.
    """
    a, b = image.shape[:2]  # Get original dimensions
    size = max(a, b)  # The target size to make it square
    
    # Compute padding (balanced on both sides)
    pad_top = (size - a) // 2
    pad_bottom = size - a - pad_top
    pad_left = (size - b) // 2
    pad_right = size - b - pad_left
    
    # Create padding format dynamically
    pad_width = [(pad_top, pad_bottom), (pad_left, pad_right)]
    
    if image.ndim == 3:  # If the image has channels, don't pad them
        pad_width.append((0, 0))
    
    # Apply padding
    padded_image = np.pad(image, pad_width, mode='constant', constant_values=0)
    
    return padded_image

def rotate_and_crop(
    label_image: np.ndarray,
    flow_direction: str,
    display_orientation: str,
    flow_direction_map: dict,
    crop_image: bool
) -> np.ndarray:
    """
    Rotate and crop an image based on flow direction and display orientation.

    Args:
        label_image (np.ndarray): Input label image.
        flow_direction (str): Flow direction of the image.
        display_orientation (str): Desired display orientation.
        flow_direction_map (dict): Mapping of flow directions to rotations.
        crop_image (bool): Whether to crop the image.

    Returns:
        np.ndarray: Rotated and cropped image.
    """
    # Read and rotate image
    rotations = 0
    if display_orientation != flow_direction:
        rotations = (flow_direction_map[display_orientation] - flow_direction_map[flow_direction]) % 4
        label_image = np.rot90(label_image, k=rotations)

    if crop_image == False:
        return label_image
    if (flow_direction == 'to_right') | (flow_direction == 'to_left'):
        dimension = 0
    else:
        dimension = 1
    crop_min, crop_max = np.where(label_image != 0)[dimension].min(), np.where(label_image != 0)[dimension].max()
    if dimension == 0:
        croped_image = label_image[crop_min:crop_max, :]
    else:
        croped_image = label_image[:, crop_min:crop_max]

    return croped_image

#-------------------------------------------------------
# Image processing functions
#-------------------------------------------------------

def periodic_decomposition(im: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform periodic decomposition of an image.

    Args:
        im (np.ndarray): Input image.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Periodic component and smooth component of the image.
    """
    im = im.astype('float32')
    # find the number of rows and cols
    N_rows, N_cols = im.shape
    # create an zero matrix the size of the image
    v = np.zeros((N_rows,N_cols))
    # fill the edges of V with the difference between the opposite edge of the real image
    v[0,:] = im[0,:] - im[-1,:]
    v[-1,:] = -v[0,:]
    v[:,0] = v[:,0] + im[:,0] - im[:,-1]
    v[:,-1] = v[:,-1] - im[:,0] + im[:,-1]
    # calculate the frequencies of the image
    fx = matlib.repmat(np.cos(2 * np.pi * np.arange(0,N_cols) / N_cols),N_rows,1)
    fy = matlib.repmat(np.cos(2 * np.pi * np.arange(0,N_rows) / N_rows),N_cols,1).T
    # set the fx[0,0] to 0 to avoid division by zero
    fx[0,0] = 0
    # calculate the smoothed image component
    s = np.real(ifft2(fft2(v) * 0.5 / (2 - fx - fy)))
    # If you want to calculate the periodic fft directly
    # p_fft = fftshift(fft2(actin) - fft2(v) * 0.5 / (2 - fx - fy))
    p = im - s

    return p, s

def least_moment(image: np.ndarray, xcoords: Optional[np.ndarray] = [], ycoords: Optional[np.ndarray] = []) -> Tuple[float, float]:
    """
    Calculate the least moment of an image.

    Args:
        image (np.ndarray): Input image.
        xcoords (Optional[np.ndarray]): X-coordinates for the image. Defaults to [].
        ycoords (Optional[np.ndarray]): Y-coordinates for the image. Defaults to [].

    Returns:
        Tuple[float, float]: Angle of axis (theta) and eccentricity.
    """
    # get the image shape
    N_rows, N_cols = image.shape

    # check if xcoords and ycoords are passed in the function
    if len(xcoords) == 0:
        # create coordinates for the image
        xcoords, ycoords = np.meshgrid(np.arange(0,N_cols) , np.arange(0,N_rows))

    #calculate the moments
    M00 = np.sum(image)
    M10 = np.sum(image * xcoords)
    M01 = np.sum(image * ycoords)
    M11 = np.sum(image * xcoords * ycoords)
    M20 = np.sum(image * xcoords * xcoords)
    M02 = np.sum(image * ycoords * ycoords)

    # center of mass
    xave = M10 / M00
    yave = M01 / M00

    # calculate the central moments
    mu20 = M20/M00 - xave**2
    mu02 = M02/M00 - yave**2
    mu11 = M11/M00 - xave*yave

    # angle of axis
    theta = 0.5 * np.arctan2((2 * mu11),(mu20 - mu02))

    # multiply by -1 to correct for origin being in top left corner instead of bottom right
    theta = -1 * theta

    # find eigenvectors
    lambda1 = (0.5 * (mu20 + mu02)) + (0.5 * np.sqrt(4 * mu11**2 + (mu20 - mu02)**2))
    lambda2 = (0.5 * (mu20 + mu02)) - (0.5 * np.sqrt(4 * mu11**2 + (mu20 - mu02)**2))
    
    # calculate the eccentricity (e.g. how oblong it is)
    eccentricity = np.sqrt(1 - lambda2/lambda1)

    return theta, eccentricity

def image_local_order(
    i: int,
    imstack: np.ndarray,
    window_size: int = 33,
    overlap: float = 0.5,
    im_mask: Optional[np.ndarray] = None,
    intensity_thresh: float = 0,
    eccentricity_thresh: float = 0,
    plot_overlay: bool = False,
    plot_angles: bool = False,
    plot_eccentricity: bool = False,
    save_figures: bool = False,
    save_path: str = ''
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate the local order of an image or stack of images.

    Args:
        i (int): Index of the image.
        imstack (np.ndarray): Stack of images.
        window_size (int): Size of the window. Defaults to 33.
        overlap (float): Overlap between windows. Defaults to 0.5.
        im_mask (Optional[np.ndarray]): Mask for the image. Defaults to None.
        intensity_thresh (float): Intensity threshold. Defaults to 0.
        eccentricity_thresh (float): Eccentricity threshold. Defaults to 0.
        plot_overlay (bool): Whether to plot overlay. Defaults to False.
        plot_angles (bool): Whether to plot angles. Defaults to False.
        plot_eccentricity (bool): Whether to plot eccentricity. Defaults to False.
        save_figures (bool): Whether to save figures. Defaults to False.
        save_path (str): Path to save figures. Defaults to ''.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]: 
        X-coordinates, Y-coordinates, U-component, V-component, theta stack, and eccentricity stack.
    """
    
    # check if an output directory is given
    if len(save_path) > 0:
        # if directory doesn't exist, make it
        if os.path.isdir(save_path) == False:
            os.mkdir(save_path)
            
    # check to see if it's a stack of images or a single image
    if len(imstack.shape) == 2:
        imstack = np.expand_dims(imstack, axis=0)
    if im_mask is not None:
        if imstack.shape[0] == 1:
            im_mask = np.expand_dims(im_mask, axis=0)

    # get the image shape
    N_images, N_rows, N_cols = imstack.shape

    # make window size off if it isn't already
    if window_size % 2 == 0:
        window_size += 1
    
    # define the radius of the window
    radius = int(np.floor((window_size) / 2))
    
    # make a list of the r,c positions for the windows
    rpos = np.arange(radius,N_rows-radius,int(window_size * overlap))
    cpos = np.arange(radius,N_cols-radius,int(window_size * overlap))

    # make a structuring element to filter the mask
    bpass_filter = disk(radius * .5)

    # make window mask
    window_mask = np.zeros((window_size, window_size))
    window_mask[int(np.floor(window_size/2)), int(np.floor(window_size/2))] = 1

    # filter the mask with the structuring element to define the ROI
    window_mask = cv2.filter2D(window_mask, -1, bpass_filter)
    window_mask = np.rint(window_mask) == 1

    # check if there is an image mask
    if im_mask is None:
        im_mask = np.ones_like(imstack).astype('bool')
    else:
        # make sure the input mask is a boolean
        im_mask = im_mask.astype('bool')

    # make x and y coordinate matrices
    xcoords, ycoords = np.meshgrid(np.arange(0,window_size) , np.arange(0,window_size))

    # length of orientation vector
    arrow_length = radius / 2

    # make lists to hold for multiple frames
    theta_stack, ecc_stack, u_stack, v_stack = [], [], [], []

    for frame,im in enumerate(imstack):

        # make empty list variables
        im_theta, im_ecc = [], []
        x, y, u, v = [], [], [], []
        
        
        # loop through each position and measure the local orientation
        for r in rpos:
            for c in cpos:
                # store the row and column positions
                x.append(c)
                y.append(r)

                # check to see if point is within image mask
                if im_mask[frame,r,c] == True:
                    # define the window to analyze
                    im_window = im[r-radius:r+radius+1,c-radius:c+radius+1]

                    # check that it's above the intensity threshold
                    if np.mean(im_window) > intensity_thresh:
                        # separate out the periodic and smooth components
                        im_window_periodic, im_window_smooth = periodic_decomposition(im_window)
                        # take the FFT of the periodic component
                        im_window_fft = fftshift(fft2(im_window_periodic))
                        # find the image norm and mulitply by the mask
                        im_window_fft_norm = image_norm(im_window_fft) * window_mask
                        # calculate the angle and eccentricity of orientation based on the FFT moments
                        theta, eccentricity = least_moment(im_window_fft_norm, xcoords, ycoords)

                        # correct for real space
                        theta = theta + np.pi/2

                        # map everything back to between -pi/2 and pi/2
                        if theta > np.pi/2:
                            theta -= np.pi

                        # filter based on eccentricity
                        if eccentricity < eccentricity_thresh:
                            eccentricity = np.nan
                            theta = np.nan

                        # add the values to each list
                        im_theta.append(theta)
                        im_ecc.append(eccentricity)
                        u.append(np.cos(theta) * arrow_length)
                        v.append(np.sin(theta) * arrow_length)
                    else:
                        im_theta.append(np.nan)
                        im_ecc.append(np.nan)
                        u.append(np.nan)
                        v.append(np.nan)
                else:
                    im_theta.append(np.nan)
                    im_ecc.append(np.nan)
                    u.append(np.nan)
                    v.append(np.nan)

        # turn all the lists into arrays for simple indexing
        x = np.array(x)
        y = np.array(y)
        u = np.array(u)
        v = np.array(v)
        im_theta = np.reshape(im_theta,(len(rpos),len(cpos)))
        im_ecc = np.reshape(im_ecc,(len(rpos),len(cpos)))

        if plot_angles:
            plt.figure()
            plt.imshow(im_theta * 180 / np.pi, vmin=-90, vmax=90, cmap='hsv')
            plt.colorbar()
            plt.title('Orientation')
            plt.show()
            if save_figures:
                plt.savefig(save_path + 'angle_map_frame_%03d' % (frame) + f'_{i}.tif', format='png', dpi=300)

        if plot_eccentricity:
            plt.figure()
            plt.imshow(im_ecc, vmin=0, vmax=1)
            plt.colorbar()
            plt.title('Eccentricity')
            plt.show()
            if save_figures:
                plt.savefig(save_path + 'eccentrcitiy_map_frame_%03d' % (frame) + f'_{i}.tif', format='png', dpi=300)

        if plot_overlay:
            plt.figure()
            plt.imshow(im, cmap='Greys_r')
            plt.quiver(x,y,u,v, color='yellow', pivot='mid', scale_units='xy', scale=overlap/2, headaxislength=0, headlength=0, width=0.005)
            plt.title('Overlay')
            plt.show()
            if save_figures:
                plt.savefig(save_path + 'overlay_frame_%03d' % (frame) + f'_{i}.tif', format='png', dpi=300)

        theta_stack.append(im_theta)
        ecc_stack.append(im_ecc)
        u_stack.append(u)
        v_stack.append(v)

    # reduce dimensions if only one frame
    if N_images == 1:
        u_stack = u_stack[0]
        v_stack = v_stack[0]
        theta_stack = theta_stack[0]
        ecc_stack = ecc_stack[0]

    return x, y, u_stack, v_stack, theta_stack, ecc_stack

#-------------------------------------------------------
# Analysis functions
#-------------------------------------------------------

def calculate_order_parameter(im_theta_stack: Union[np.ndarray, List[np.ndarray]], neighborhood_radius: int) -> Union[float, List[float]]:
    """
    Calculate the order parameter for an image or stack of images.

    Args:
        im_theta_stack (Union[np.ndarray, List[np.ndarray]]): Stack of theta matrices.
        neighborhood_radius (int): Radius of the neighborhood.

    Returns:
        Union[float, List[float]]: Order parameter(s).
    """
    # check if it's a list
    if type(im_theta_stack) == np.ndarray:
        im_theta_stack = [im_theta_stack]

    # get the number of images
    N_images = len(im_theta_stack)

    # empty list to hold order parameter values for each image
    im_orderparameter_stack = []

    for im_theta in im_theta_stack:
        order_list = []

        rpos = np.arange(neighborhood_radius,im_theta.shape[0]-neighborhood_radius)
        cpos = np.arange(neighborhood_radius,im_theta.shape[1]-neighborhood_radius)

        for r in rpos:
            for c in cpos:
                search_window = im_theta[r-neighborhood_radius:r+neighborhood_radius+1,c-neighborhood_radius:c+neighborhood_radius+1]
                vector_array = np.ones_like(search_window) * im_theta[r,c]
                order_array = np.cos(search_window - vector_array) ** 2 - 0.5
                order_array = np.delete(order_array, (2*neighborhood_radius + 1)**2 // 2)
                if not np.isnan(order_array).all() == True:
                    if order_array.size > 0:
                        order_list.append(2 * np.nanmean(order_array))

        im_orderparameter_stack.append(np.nanmedian(order_list))

    # Reduce dimensions if only one image
    if N_images == 1:
        im_orderparameter_stack = im_orderparameter_stack[0]

    return im_orderparameter_stack

def parameter_search(
    image_list: List[str],
    min_win_size: int,
    win_size_interval: int,
    overlap: float,
    flow_direction_map: dict,
    flow_direction: str = 'to_right',
    display_orientation: str = 'to_right',
    max_win_size: Optional[int] = None,
    plot_figure: bool = True,
    from_cellpose_segmentation: bool = False,
    crop_image: bool = False,
    pad: bool = False
) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Perform parameter search for window size and neighborhood radius.

    Args:
        image_list (List[str]): List of image file paths.
        min_win_size (int): Minimum window size.
        win_size_interval (int): Interval for window size.
        overlap (float): Overlap between windows.
        flow_direction_map (dict): Mapping of flow directions.
        flow_direction (str): Flow direction of the image. Defaults to 'to_right'.
        display_orientation (str): Desired display orientation. Defaults to 'to_right'.
        max_win_size (Optional[int]): Maximum window size. Defaults to None.
        plot_figure (bool): Whether to plot the results. Defaults to True.
        from_cellpose_segmentation (bool): Whether the input is from Cellpose segmentation. Defaults to False.
        crop_image (bool): Whether to crop the image. Defaults to False.
        pad (bool): Whether to pad the image. Defaults to False.

    Returns:
        Tuple[pd.DataFrame, np.ndarray]: Dataframe of results and window-neighborhood matrix.
    """
    # turn off warning for division by NaN
    np.seterr(divide='ignore', invalid='ignore')
    
    # read in an image to get the image shape
    if from_cellpose_segmentation:
        im_temp = np.load(image_list[0], allow_pickle = True).item().get('outlines')
        im      = np.where(im_temp==0,0,1)
    else:
        im = io.imread(image_list[0])

    

    if flow_direction != display_orientation:
        im = rotate_and_crop(im, flow_direction, display_orientation, flow_direction_map, crop_image)
    
    # Check if a max_win_size is assigned and if it is usable
    if max_win_size is None:
        # set max win size - this is arbitary but would give you 3 windows with an overlap of 50%
        max_win_size = (np.max(im.shape) -1 ) // 3 # -1 accounts for rare case where image is perfectly divisble by 3 into an even number
    elif max_win_size < min_win_size or max_win_size > ((np.max(im.shape) -1 ) // 3):
        # set max win size - this is arbitary but would give you 3 windows with an overlap of 50%
        print('The seelcted max_win_size is not usable. Has to be larger then min_win_size but smaller then (np.max(im.shape) -1 ) // 3. \n' 
              '(np.max(im.shape) -1 ) // 3 is the maximal possible winsize and used instead')
        max_win_size = (np.max(im.shape) -1 ) // 3 # -1 accounts for rare case where image is perfectly divisble by 3 into an even number
    # make a list of window sizes to search through
    win_size_list = np.arange(min_win_size, max_win_size, win_size_interval)
    
    # make them odd if they aren't already
    win_size_list[win_size_list % 2 == 0] += 1

    # make empty lists to hold the data
    win_size_result, image_result, order_parameter_result, neighborhood_result = [], [], [], []
    
    # for loops to go through each image, window size and neighborhood radius
    for idx, image in enumerate(image_list):
        # Read in the image
        if from_cellpose_segmentation:
            im_temp = np.load(image, allow_pickle = True).item().get('outlines')
            im      = np.where(im_temp==0,0,1)
        else:
            im      = io.imread(image).astype('float32')
            
        if pad:
            im = pad_to_square(im)

        if flow_direction != display_orientation:
            im = rotate_and_crop(im, flow_direction, display_orientation, flow_direction_map, crop_image)
        
        # loop through the different window sizes
        for win_size in win_size_list:
            # calculate the theta matric       
            _,_,_,_,im_theta,_ = image_local_order(idx, im, window_size = win_size, overlap = overlap, plot_overlay = False, plot_angles=False, plot_eccentricity=False)
            
            # get the number of windows
            n_windows = np.max(im_theta.shape)
            
            # make a list of neighborhood radii to measure
            neighborhood_list = np.arange(1, ((n_windows - 1) // 2) + 1)
            
            # loop through all neighborhood radiui
            for neighborhood in neighborhood_list:
                # calculate the order parameter
                im_orderparameter = calculate_order_parameter(im_theta, neighborhood_radius=neighborhood)
                
                # add the parameters to the relevant lists
                win_size_result.append(win_size)
                image_result.append(image)
                order_parameter_result.append(im_orderparameter)
                neighborhood_result.append(neighborhood)

    # make a dictionary of our lists
    data_dict = {
        'window_size' : win_size_result,
        'neighborhood_radius' : neighborhood_result,
        'order_parameter' : order_parameter_result,
        'image' : image_result
    }
    
    # convert the dictionary to a dataframe
    Order_dataframe = pd.DataFrame(data_dict)  
    
    # Find the number of different windows and neighborhoods tested 
    N_windows = len(win_size_list)
    neighborhood_list = np.arange(1,np.max(Order_dataframe['neighborhood_radius'])+1)
    N_neighborhood = len(neighborhood_list)

    # make a matrix of window size (rows) and neighborhood size (cols)
    window_neighborhood = np.empty((N_windows, N_neighborhood))
    window_neighborhood[:] = np.nan

    # populate the matrix with the median of the order parameter of all the images
    for i,win_size in enumerate(win_size_list):
        for j, neighborhood in enumerate(neighborhood_list):
            temp = Order_dataframe[(Order_dataframe.window_size==win_size) & (Order_dataframe.neighborhood_radius==neighborhood)]
            if len(temp) > 0:
                window_neighborhood[i,j] = np.median(temp.order_parameter)
    
    if plot_figure:
        # get labels for the plots
        win_size_labels = np.unique(Order_dataframe['window_size'])
        neighborhood_labels = np.unique(Order_dataframe['neighborhood_radius'])

        plt.figure()
        plt.imshow(window_neighborhood, vmin=0, vmax=1)
        plt.colorbar(orientation='horizontal')
        plt.yticks(np.arange(0,len(win_size_labels)),labels=win_size_labels)
        plt.ylabel('Window Size (px)')
        plt.xticks(np.arange(0,len(neighborhood_labels),2),labels=neighborhood_labels[::2], rotation='vertical')
        plt.xlabel('Neighbourhood Radius (px)')
        plt.title('Median order parameter')
        plt.tight_layout()
        plt.show()

    return Order_dataframe, window_neighborhood

def parameter_comparison(
    Order_dataframe1: pd.DataFrame,
    window_neighborhood1: np.ndarray,
    Order_dataframe2: pd.DataFrame,
    window_neighborhood2: np.ndarray,
    save_figures: bool = False,
    save_path: str = ''
) -> Tuple[np.ndarray, np.ndarray, List[int], List[int]]:
    """
    Compare parameters between two datasets.

    Args:
        Order_dataframe1 (pd.DataFrame): Dataframe of the first dataset.
        window_neighborhood1 (np.ndarray): Window-neighborhood matrix of the first dataset.
        Order_dataframe2 (pd.DataFrame): Dataframe of the second dataset.
        window_neighborhood2 (np.ndarray): Window-neighborhood matrix of the second dataset.
        save_figures (bool): Whether to save figures. Defaults to False.
        save_path (str): Path to save figures. Defaults to ''.

    Returns:
        Tuple[np.ndarray, np.ndarray, List[int], List[int]]: Order difference, p-value matrix, window sizes, and neighborhood sizes.
    """
    # Find the number of different windows and neighborhoods tested 
    win_size_list = sorted(np.unique(Order_dataframe1['window_size']))
    N_windows = len(win_size_list)
    neighborhood_list = sorted(np.unique(Order_dataframe1['neighborhood_radius']))
    N_neighborhood = len(neighborhood_list)

    # make a matrix of window size (rows) and neighborhood size (cols)
    p_median = np.empty((N_windows, N_neighborhood))
    p_median[:] = np.nan

    # populate the matrix with the median of the order parameter of all the images
    for i,win_size in enumerate(win_size_list):
        for j, neighborhood in enumerate(neighborhood_list):
            # select the values we need
            temp1 = Order_dataframe1[(Order_dataframe1.window_size==win_size) & (Order_dataframe1.neighborhood_radius==neighborhood)]
            temp2 = Order_dataframe2[(Order_dataframe2.window_size==win_size) & (Order_dataframe2.neighborhood_radius==neighborhood)]
            # make sure the list isn't full of nan
            if (len(temp1) > 0) and (len(temp2) > 0):
                # calculate the p-value between the two groups
                _, p_median[i,j] = mannwhitneyu(temp1.order_parameter, temp2.order_parameter)


    # difference between order parameters
    order_diff = window_neighborhood1-window_neighborhood2
    # plot difference in order parameter
    plt.figure()
    plt.imshow(order_diff, cmap='jet')
    plt.xlabel('Neighborhood (Vectors)')
    plt.xticks(np.arange(0,len(neighborhood_list),2),labels=neighborhood_list[::2], rotation='vertical')
    plt.ylabel('Window Size (px)')
    plt.yticks(np.arange(0,len(win_size_list)),labels=win_size_list)
    plt.colorbar(orientation='horizontal')
    plt.title('Order Parameter Difference (1st sample - 2nd sample)')
    plt.show()
    if save_figures:
        plt.savefig(save_path + 'parameter_search_difference.png', format='png', dpi=300)

    # plot p-value comparison
    plt.figure()
    plt.imshow(p_median, cmap='spring_r')
    plt.xlabel('Neighborhood (Vectors)')
    plt.xticks(np.arange(0,len(neighborhood_list),2),labels=neighborhood_list[::2], rotation='vertical')
    plt.ylabel('Window Size (px)')
    plt.yticks(np.arange(0,len(win_size_list)),labels=win_size_list)
    plt.title('P-value comparison')
    plt.colorbar(orientation='horizontal')
    plt.show()
    if save_figures:
        plt.savefig(save_path + 'parameter_search_p_value.png', format='png', dpi=300)
        
    return order_diff, p_median, win_size_list, neighborhood_list
