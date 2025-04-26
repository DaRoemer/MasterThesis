%% AFT Batch Processing Script %%
% This script processes .tif image files for nematic ordering analysis.
% It supports batch processing of files in a specified directory and its subfolders.
%
% Key Features:
% - User-defined input options for directory paths and analysis parameters.
% - Automatic creation of output folders with incremental suffixes.
% - Saves results as .mat and .csv files.
%
% Usage:
% 1. Set the `parent_d` variable to the directory containing your input files.
% 2. Configure analysis parameters in the "User-defined Input Options" section.
% 3. Run the script to process all .tif files in the specified directory.
%
% Marcotti et al. 2021
% DOI: 10.3389/fcomp.2021.745831
% Adapted by: Felix Romer
% Last update: 26.04.2025

%% Initialization %%
clear; clc
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'DefaultFigureVisible', 'off');
warning off

%% User-defined Input Options %%
% Specify the parent directory containing input files
parent_d = '/path/to/your/input/files'; % Example: '/media/felix/Felixdrive/Analysis2025/...'
subfolders = false; % Set to true if input files are in subfolders

% Define analysis parameters
parameters.overlap                  = 0.5;  % Overlap between windows
parameters.min                      = 0;    % Minimum intensity threshold
parameters.max                      = 255;  % Maximum intensity threshold
parameters.figures                  = 1;    % Enable/disable figure generation
parameters.checkpoint               = 0;    % Enable/disable checkpointing
parameters.mask_method              = 0;    % Mask method (0: none, 1: local)
parameters.filter_blank             = 0;    % Filter blank images
parameters.filter_isotropic         = 0;    % Filter isotropic regions
parameters.eccentricity_threshold   = 0;    % Eccentricity threshold

% Define parameter ranges for analysis
radius = [1]; % Example: [1, 2, 3]
winsize = [125]; % Example: [125, 150, 175]

%% Load Input Files %%
% Get all subfolders in the main directory if subfolders are enabled
if subfolders
    allFolders = dir(parent_d);
    allFolders = allFolders([allFolders.isdir]); % Only keep directories
    allFolders = allFolders(~ismember({allFolders.name}, {'.', '..', 'Backup', 'Analysis', 'Figure', 'Python_output'}));
else
    allFolders = [1];
end

%% Process Each Folder %%
for folder = 1:length(allFolders)
    fprintf('\nProcessing Folder %d of %d:', folder, length(allFolders))
    matlab_folder = cd;
    if subfolders
        currentdir = fullfile(parent_d, allFolders(folder).name);
    else
        currentdir = parent_d;
    end
    cd(currentdir)
    listing = dir('*.tif');
    n_files = length(listing);
    
    % Loop Through Parameter Pairs %%
    for parameter = 1:length(radius)
        fprintf('\nProcessing Parameter %d of %d:', parameter, length(radius))
        parameters.st = radius(parameter);
        parameters.winsize = winsize(parameter);

        % Create the output folder with incremental suffixes if needed
        output_folder_name = 'output_AFT_1';
        suffix = 0;
        cd(currentdir)
        while exist(output_folder_name, 'dir')
            suffix = suffix + 1;
            output_folder_name = ['output_AFT_' num2str(suffix)];
        end
        mkdir(output_folder_name);
        output_folder_name = ['/' output_folder_name];
        parameters.output = output_folder_name;
        cd(matlab_folder)
        
        % Analyze Each File %%
        av_ordermat = zeros(n_files, 1);
        strForm = sprintf('%%.%dd', length(num2str(n_files)));
        fprintf('\nAnalyzing file ');
        
        for i = 1:n_files
            % Display progress
            procStr = sprintf(strForm, i);
            fprintf(1, [procStr, ' of ', num2str(n_files)]);
            
            file = listing(i).name;
            directory = listing(i).folder;
            
            % File and directory name for mask (if local)
            if parameters.mask_method == 1
                file_mask = listing_masks(i).name;
                directory_mask = listing_masks(i).folder;
                parameters.mask_name = fullfile(directory_mask, file_mask);
            end
           
            % Call analysis function
            av_ordermat(i, 1) = AFT_function_adapted_V2(file, directory, parameters, true);
        
            if i < n_files
                for jj = 1:(length(procStr) + 4 + length(num2str(n_files)))
                    fprintf(1, '\b');
                end
            end
        end
        
        % Save Results %%
        save(fullfile([currentdir output_folder_name], 'median_order_parameter.mat'), 'av_ordermat');
        T = table(av_ordermat);
        T.Properties.VariableNames = {'median_order_parameter'};
        writetable(T, fullfile([currentdir output_folder_name], 'median_order_parameter.csv'))
    end
end