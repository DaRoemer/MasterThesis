#-------------------------------------------------------
# AFT_tools.py
# By: Felix Romer
# Last updated: 24.04.2025
# Description:
# This script contains functions to load, process, and analyze data from MATLAB AFT files.
# It includes functions for loading data, performing statistical tests, and plotting results.
# The script is designed to work with data from the AFT (Alignment by fourier transform) analysis.
# It is intended for use in the analysis of cell alignment and order parameters in microscopy images.
#-------------------------------------------------------

# ------------------------------------------------------
# Import necessary libraries
# ------------------------------------------------------
import os
import numpy as np
import pandas as pd
import scipy.stats as ss
import scikit_posthocs as sp
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.io as sio
from typing import Tuple
from typing import Optional

# ------------------------------------------------------
# Data Loading and Processing
# ------------------------------------------------------
def load_and_process_matlab_data(
    main_dr: str,
    folder_names: list,
    windowsize: int,
    factor_px_to_um: float,
    av_cell_size_um2: float,
    remove_files: Optional[list] = None,
    replicates_to_use: Optional[list] = None,
    cell_count_path: Optional[str] = None
) -> Tuple[pd.DataFrame, float]:
    """
    Load and process data from specified folders from matlab AFT to prepare for analysis.

    Parameters:
    - main_dr (str): Main directory containing the data.
    - folder_names (list): List of folder names corresponding to conditions.
    - windowsize (int): The window size in pixels.
    - factor_px_to_um (float): Conversion factor from pixels to micrometers.
    - av_cell_size_um2 (float): Average cell size in square micrometers.
    - remove_files (list, optional): List of file names to exclude.
    - replicates_to_use (list, optional): List of replicates to include, if empty all replicates are included.
    - cell_count_path (str, optional): Path to cell count CSV file.

    Returns:
    - tuple: DataFrame with filtered and processed data and neighbourhood size in micrometers.
    """
    # Read av_ordermat1 CSVs (assuming 6 CSVs, one per cell)
    av_ordermats = []
    Replicates = {}
    FileName = {}

    for folder_name in folder_names:
        folder_path = os.path.join(main_dr, folder_name)
        all_files = os.listdir(folder_path)
        tif_files = [f for f in all_files if f.endswith('.tif')]
        extracted_data = []
        for file_name in tif_files:
            parts = file_name.split('_')
            if len(parts) > 0:
                prefix = parts[0]
                suffix = parts[1][:3]
                if suffix == 'RFP':
                    prefix = prefix + '_RFP'
                extracted_data.append(prefix)
            else:
                extracted_data.append(file_name.split('.')[0])
        Replicates[folder_name] = extracted_data
        FileName[folder_name] = tif_files

        file_path = os.path.join(main_dr, folder_name, 'output_parameter_search', 'median_order_parameter_search.mat')
        av_ordermat = sio.loadmat(file_path)
        av_ordermats.append(av_ordermat)

    # Read the parameters CSV
    parameters_file_path = os.path.join(main_dr, folder_names[0], 'output_parameter_search', 'parameters.txt')
    parameters_df = pd.read_csv(parameters_file_path)

    # Read the cell count CSV if provided
    if cell_count_path:
        cell_count_df = pd.read_csv(cell_count_path)
        cell_count_df['Condition'] = cell_count_df['Condition'].str.capitalize()
        # Map the replicate values
        replicate_mapping = {20230227: 1, 20230322: 2, 20230324: 3}
        cell_count_df['replicate'] = cell_count_df['Replicate'].map(replicate_mapping)
        # Add 4 new rows with replicate 1 and NaN for the rest
        new_rows = pd.DataFrame({'Number of cells': [np.nan]*4, 'Condition': 'Infected', 'Replicate': [np.nan]*4, 'replicate': [1]*4})
        cell_count_df = pd.concat([cell_count_df.iloc[:4], new_rows, cell_count_df.iloc[4:]]).reset_index(drop=True)
    else:
        cell_count_df = None

    # Process the data
    df, neighbourhood_um = process_data(
        parameters_df=parameters_df,
        Replicates=Replicates,
        FileName=FileName,
        windowsize=windowsize,
        av_ordermats=av_ordermats,
        folder_names=folder_names,
        remove_files=remove_files,
        use_replicates=replicates_to_use,
        cell_count_df=cell_count_df,
        factor_px_to_um=factor_px_to_um,
        av_cell_size_um2=av_cell_size_um2
    )

    return df, neighbourhood_um

def process_data(
    parameters_df: pd.DataFrame,
    Replicates: dict,
    FileName: dict,
    windowsize: int,
    av_ordermats: list,
    folder_names: list,
    remove_files: Optional[list] = None,
    use_replicates: Optional[list] = None,
    cell_count_df: Optional[pd.DataFrame] = None,
    factor_px_to_um: float = 1.,
    av_cell_size_um2: float = 1.
) -> Tuple[pd.DataFrame, float]:
    """
    Process data to calculate order parameters, radial distances, and cell counts for each condition.

    Parameters:
    - parameters_df (pd.DataFrame): DataFrame containing parameters for analysis.
    - Replicates (dict): Dictionary of replicates for each condition.
    - FileName (dict): Dictionary of file names for each condition.
    - windowsize (int): The window size in pixels.
    - av_ordermats (list): List of average order matrices.
    - folder_names (list): List of folder names corresponding to conditions.
    - remove_files (list, optional): List of file names to remove.
    - use_replicates (list, optional): List of replicates to include.
    - cell_count_df (pd.DataFrame, optional): DataFrame of cell counts for each condition.
    - factor_px_to_um (float, optional): Conversion factor from pixels to micrometers.
    - av_cell_size_um2 (float, optional): Average cell size in square micrometers.

    Returns:
    - tuple: Combined DataFrame and neighbourhood size in micrometers.
    """
    # Extract order parameter for the requested window size
    array = np.arange(int(parameters_df['min_winsize_px']), int(parameters_df['max_winsize_px']+parameters_df['winsize_int_px']), int(parameters_df['winsize_int_px']))
    ind = np.where(array == windowsize)[0][0]

    df = pd.DataFrame()

    for i, av_order_mat in enumerate(av_ordermats):
        av_order_df = pd.DataFrame(av_order_mat['av_ordermat_output'][ind][0])
        condition = folder_names[i]
        replicates = Replicates[condition]
        file_names = FileName[condition]
        # Change name to 3x, 5x, 7x, 9x, 11x, 13x,...
        av_order_df.columns = [f'{(i+1)*2+1}x' for i in range(av_order_df.shape[1])]

        # Add replicate label
        av_order_df['replicate'] = replicates
        av_order_df['File name'] = file_names
        # Remove rows where File name is in remove_files
        if remove_files is not None:
            av_order_df = av_order_df[~av_order_df['File name'].isin(remove_files)]
        
        # Keep only rows where replicate is in use_replicates
        if use_replicates is not None:
            av_order_df = av_order_df[av_order_df['replicate'].isin(use_replicates)]

        # Add cell count if provided, then melt df
        if cell_count_df is not None:
            temp_cell_count_df = cell_count_df[cell_count_df['Condition'] == condition]
            av_order_df['Number of cells'] = temp_cell_count_df['Number of cells'].values
            melted_df = pd.melt(av_order_df, id_vars=['File name', 'replicate', 'Number of cells'], var_name='Neighbourhood', value_name='Order parameter')
        else:
            melted_df = pd.melt(av_order_df, id_vars=['File name', 'replicate'], var_name='Neighbourhood', value_name='Order parameter')

        # Add condition column
        melted_df['Group'] = condition

        # Append to the combined DataFrame
        df = pd.concat([df, melted_df], ignore_index=True)

    # Calculate radial distance, window size, and number of cells per neighbourhood
    neighbourhood_um = windowsize * factor_px_to_um
    df['Radial Distance (µm)'] = df['Neighbourhood'].apply(lambda x: int(((int(x.split('x')[0]) -1) /2) * neighbourhood_um / 2))
    df['Window size (um2)'] = (df['Radial Distance (µm)'] *2)**2
    df['Cells_per_neighbourhood'] = (df['Window size (um2)'] / av_cell_size_um2)

    # Reset index
    df.reset_index(drop=True, inplace=True)

    return df, neighbourhood_um

def process_python_AFT_data(df, 
                            windowsize:       int, 
                            factor_px_to_um:  float, 
                            av_cell_size_um2: float,
                            use_replicates:   Optional[list] = None,
                            remove_files:     Optional[list] = None
) -> pd.DataFrame:
    """
    Process data to calculate order parameters, radial distances, and cell counts for each condition.

    Parameters:
    - df (pd.DataFrame or str): DataFrame containing parameters for analysis or a path to a CSV file.
    - windowsize (int): The window size in pixels used for filtering data.
    - factor_px_to_um (float): Conversion factor from pixels to micrometers.
    - av_cell_size_um2 (float): Average cell size in square micrometers.
    - use_replicates (list, optional): List of replicate identifiers to include in the analysis. Defaults to an empty list (no filtering).
    - remove_files (list, optional): List of file names to exclude from the analysis. Defaults to an empty list.

    Returns:
    - pd.DataFrame: Processed DataFrame containing the calculated order parameters, radial distances, and cell counts.

    Raises:
    - ValueError: If the input data is neither a DataFrame nor a valid path to a CSV file.
    - ValueError: If the specified window size is not found in the DataFrame.

    Notes:
    - The function extracts metadata such as 'Date', 'Channel', and 'Replicate' from the 'Image' column.
    - Filters data based on the provided `use_replicates` and `remove_files` lists.
    - Adds columns for radial distance, window size, and cell counts based on input parameters.
    """
    
    # Check if df is a DataFrame
    if not isinstance(df, pd.DataFrame) and not isinstance(df, str):
        raise ValueError('The input data must be a DataFrame or a path to a CSV file.')
    elif isinstance(df, str):
        df = pd.read_csv(df)

    # FCheck if windowsize is in the DataFrame
    if windowsize not in df['Window size'].unique():
        raise ValueError(f'Window size {windowsize} not found in the DataFrame.')
    
    # Filter for the requested window size
    df = df[df['Window size'] == windowsize]

    # Extract metadata from file names
    df['Date']      = df['Image'].str.split('_').str[0].str.split('/').str[-1]
    df['Channel']   = df['Image'].str.split('.').str[-2].str.split('_').str[-2]
    df['Replicate'] = df['Image'].str.split('_').str[-1].str.split('.').str[-1]

    # Filter 
    if use_replicates is not None:
        df = df[df['replicate'].isin(use_replicates)]
    if remove_files is not None:
        df = df[~df['Image'].isin(remove_files)]

    # Convert from px to um
    df['Radial Distance (µm)'] = df['Neighborhood radius'] * df['Window size'] * factor_px_to_um / 2
    df['Window size (µm2)'] = (df['Radial Distance (µm)'] * 2) **2
    df['Cells per neighboorhood'] = df['Window size (µm2)'] / av_cell_size_um2
    df.reset_index(drop=True, inplace=True)

    return df

# ------------------------------------------------------
# Statistical Analysis
# ------------------------------------------------------
def perform_kruskal_conover(
    df: pd.DataFrame,
    value_col: str = 'Order parameter',
    group_col: str = 'Group',
    neighbourhood_col: str = 'Neighbourhood'
) -> Tuple[dict, pd.DataFrame]:
    """
    Perform Kruskal-Wallis H-test and Conover post-hoc test for each neighbourhood in the provided DataFrame.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the data for analysis.
    - value_col (str): Column name for the values being compared (e.g., 'Order parameter').
    - group_col (str): Column name for the group labels (e.g., 'Group').
    - neighbourhood_col (str): Column name for the neighbourhoods (e.g., 'Neighbourhood').

    Returns:
    - tuple: Dictionary of Kruskal-Wallis test results and DataFrame with Conover post-hoc test results.
    """
    # Perform Kruskal-Wallis H-test for each neighbourhood
    kruskal_results = {}
    for neighbourhood in df[neighbourhood_col].unique():
        subset = df[df[neighbourhood_col] == neighbourhood]
        data = [subset[subset[group_col] == group][value_col].values for group in subset[group_col].unique()]
        H, p = ss.kruskal(*data)
        kruskal_results[neighbourhood] = {'H-statistic': H, 'p-value': p}

    # Perform Conover post-hoc test for each neighbourhood
    conover_results = {}
    for neighbourhood in df[neighbourhood_col].unique():
        subset = df[df[neighbourhood_col] == neighbourhood]
        result = sp.posthoc_conover(subset, val_col=value_col, group_col=group_col, p_adjust='holm')
        conover_results[neighbourhood] = result

    # Build a DataFrame for Conover post-hoc results
    conover_df = []
    for neighbourhood, result in conover_results.items():
        groups = df[group_col].unique()        
        for i, group1 in enumerate(groups):
            for group2 in groups[i+1:]:  # Avoid duplicate comparisons
                try:
                    p_value = result.loc[group1, group2]
                except KeyError:
                    p_value = np.nan
                significance = 'ns'  # Default to not significant
                if p_value < 0.001:
                    significance = '***'
                elif p_value < 0.01:
                    significance = '**'
                elif p_value < 0.05:
                    significance = '*'
                conover_df.append([neighbourhood, f'{group1} - {group2}', p_value, significance])
    

    conover_df = pd.DataFrame(conover_df, columns=[neighbourhood_col, 'Comparison', 'p-value', 'Significance'])

    return kruskal_results, conover_df

# ------------------------------------------------------
# Plotting Functions
# ------------------------------------------------------
def plot_by_windowsize(
    df: pd.DataFrame,
    windowsize: Optional[int] = None,
    unit: str = 'Neighbourhood',
    group_col: str = 'Group',
    value_col: str = 'Order parameter',
    custom_palette: Optional[dict] = None,
    custom_ticks: Optional[list] = None,
    custom_ticks_cell: Optional[list] = None,
    group_order: Optional[list] = None,
    start_at: int = 0,
    show_raw_data: bool = True,
    secondary_label_col: str = 'Cells_per_neighbourhood',
    second_label_pos: Optional[int] = None,
    av_cell_size_um2: float = 1,
    doge_condtions: bool = False,
    y_axis_lim: Optional[Tuple[float, float]] = None,
    linestyle_map: Optional[dict] = None
) -> plt.Figure:
    """
    Plot the order parameter by neighbourhood and condition with options for customization.

    Parameters:
    - df (pd.DataFrame): The DataFrame containing data to be plotted.
    - windowsize (int): The window size (in pixels) used in the experiment.
    - unit (str): The unit for the x-axis (default is 'Neighbourhood').
    - start_at (int): Starting position of the x-axis (default is 0).
    - second_label_pos (int, optional): Position of the secondary x-axis labels.
    - group_col (str): Column name for the grouping variable (e.g., 'Group').
    - value_col (str): Column name for the y-axis values (e.g., 'Order parameter').
    - secondary_label_col (str): Column name for secondary x-axis labels (e.g., 'Cells_per_neighbourhood').

    Returns:
    - plt.Figure: The resulting plot figure.
    """
    # Determine unique neighborhoods and x-values
    if windowsize is not None:
        df = df[df['Window size'] == windowsize]
    num_neighbourhoods = len(df[unit].unique())

    # If no colour patter is given, automatically create one determine group order dynamically
    if custom_palette is None:
        if group_order is None:
            group_order = df[group_col].unique()
        custom_palette = sns.color_palette("tab10", len(group_order))
    else:
        group_order = list(custom_palette.keys())
    if linestyle_map is None:
        linestyle_map = {group: '-' for group in group_order}
    linestyles = [linestyle_map[group] for group in group_order]

    # Create the boxplot
    fig, ax1 = plt.subplots(figsize=(10, 6))
    sns.pointplot(x=unit, y=value_col, hue=group_col, data=df, 
                  hue_order=group_order,
                  palette=custom_palette,  # Dynamically created color palette
                  dodge=doge_condtions,               # Do not separate points horizontally
                  join=False,                # Do not join the points with lines
                  ci=95,                     # Show 95% confidence interval as error bars
                  markers='o',               # Marker style
                  scale=1,                   # Scale of the markers
                  ax=ax1,
                  errwidth=2,                # Width of the error bars
                  capsize=0.2,               # Add horizontal lines at the end of the bars
                  )               

    if show_raw_data:
        # Overlay raw data points using stripplot
        sns.stripplot(x=unit, y=value_col, hue=group_col, data=df, 
                    hue_order=group_order,
                    jitter=True,      
                    dodge=doge_condtions,       
                    alpha=0.3,        
                    palette=custom_palette,  
                    marker='o', 
                    size=5,
                    ax=ax1, 
                    legend=False)      

    # Set the x-axis ticks and labels
    x_min_pos = start_at - 0.5
    x_max_pos = num_neighbourhoods - 0.5
    min_val = min(df[unit])
    max_val = max(df[unit])
    diff_val = max_val - min_val
    if custom_ticks is None:
        custom_ticks = [0, 50, 100, 150, 200]
    custom_tick_positions = [(num_neighbourhoods - 1) * (x - min_val) / diff_val for x in custom_ticks]
    print(custom_tick_positions)
    ax1.set_xticks(custom_tick_positions)
    ax1.set_xticklabels(custom_ticks)   

    # Add secondary x-axis labels if specified
    if second_label_pos is not None and secondary_label_col in df.columns:
        if custom_ticks_cell is None:
            custom_ticks_cell = [1, 10, 30, 60, 100]
        custom_tick_positions = [(num_neighbourhoods - 1) * (np.sqrt(x*av_cell_size_um2)/2 - min_val) / diff_val for x in custom_ticks_cell]
        ax2 = ax1.secondary_xaxis('bottom')
        ax2.set_xticks(custom_tick_positions)
        ax2.set_xticklabels(custom_ticks_cell) 
        ax2.set_xlabel('Approx. number of cells in one compared region')
        ax2.spines['bottom'].set_position(('outward', second_label_pos))

    # Customize the plot
    plt.title(f'Order Parameter by {unit} and Condition with windowsize of {windowsize} px', fontsize=14, pad=20)
    plt.xlabel(unit)
    plt.ylabel('Order Parameter (-)')
    plt.legend(title='Group', loc='upper right')
    # Set y-axis limits
    if y_axis_lim is None:
        plt.ylim(0, 1)
    else:
        plt.ylim(y_axis_lim)
     
    plt.xlim(x_min_pos, x_max_pos)
    sns.despine(ax=ax1, top=True, right=True)

    # Finalize layout and return the figure
    plt.tight_layout()
    return fig


def add_stat_annotation(
    fig: plt.Figure,
    conover_df: pd.DataFrame,
    start_at: int = 0
) -> plt.Figure:
    """
    Add statistical annotation to a plot based on Conover post-hoc test results.

    Parameters:
    - fig (plt.Figure): The plot figure to annotate.
    - conover_df (pd.DataFrame): DataFrame containing Conover test results.
    - start_at (int): Start annotation at this x-axis tick index.

    Returns:
    - plt.Figure: The annotated figure.
    """
    map = {'Control': -0.3, 'Uninfected': 0, 'Infected': 0.3}
    ax = fig.axes[0]
    y_max = 1.1
    counter = 0
    for tick in ax.get_xticks():
        subset = conover_df.iloc[counter*3:(counter+1)*3]
        counter += 1
        if counter <= start_at:
            continue
        for idx, row in subset.iterrows():
            group1, group2 = row['Comparison'].split(' - ')
            significance = row['Significance']
            x1, x2 = tick + map[group1], tick + map[group2]
            if x1 > x2:
                x1, x2 = x2, x1
            x1 += 0.05
            x2 -= 0.05
            ax.plot([x1, x2], [y_max - 0.1, y_max - 0.1], color='black', lw=0.5)
            ax.text((x1 + x2) / 2, y_max - 0.11, significance, ha='center', va='bottom', color='black')
    return fig

