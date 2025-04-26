%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code and documentation can be found at
% https://github.com/OakesLab/AFT-Alignment_by_Fourier_Transform
%
% This routine uses a vector field of alignment directions using small
% sub-windows in the real space image to calculate an alignment order
% parameter.  Images should be grayscale images.
% Angles determined are oriented as follows:
%
%                            ^  180�
%                            |
%                            |
%                            ------>  90�
%                            |
%                            |
%                            v  0�
%
% Order Parameter value:   0 = Completely random alignment
%                          1 = Perfectly aligned
%
% All rights and permissions belong to:
% Patrick Oakes
% poakes@gmail.com
% 05/15/2015
%
% Citation:
% Cetera et al. Nat Commun 2014; 5:1-12
% http://www.ncbi.nlm.nih.gov/pubmed/25413675
%
% Adapted by Felix Romer
% Last update: 26.04.2025
% This function is a modified version of the original AFT_function.m
% function. It includes additional features such as:
% - Improved handling of image padding
% - Enhanced output options
% - More flexible parameter settings
% - Additional comments and documentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function av_ordermat = AFT_function_adapted_V2(file, directory, parameters, pading)
    arguments
            file
            directory
            parameters
            pading = true
        end
    % DIfference V1 to V2: V1 crop the image with vectorfield, wethere V2 adds
    % a black border to the heatmap
    
    % load image
    im = imread(fullfile(directory, file));
    im = im2double(im);
    
    output_dir = [directory, parameters.output];
    
    % calculate angle vector field
    [anglemat,pc,pr,vc2,ur2] = AFT_anglemat(im, parameters);
    winsize = parameters.winsize;
    maxwin = max(pc)+winsize/2;
    diff_max = size(im,2) - maxwin;
    minwin = min(pc)-winsize/2;
    diff_min = 0 + minwin;
    
    % Scale vc2 and ur2 to always point upwards
    %vc2(ur2 < 0) = -vc2(ur2 < 0);
    %ur2 = abs(ur2);
    
    % Save pc, pr, vc2 and ur2 as .mat
    save([output_dir '/quiver_data_' file(1:end-4) '.mat'], 'pc', 'pr', 'vc2', 'ur2');
    
    % plots     
    if parameters.figures == 1
        
        % vector field
        str = '\nCreate vector image';
        lenghtstr = length(str);
        fprintf(str)
        figure
        blackImage = zeros(size(im), 'like', im);
        se                = strel('disk', 1);
        im = imdilate(im, se);
        imshow(im, [parameters.min/65535 parameters.max/65535]) % divide by 65535 to nrmalize limits
        hold on
        if exist('pc','var') ~= 0
            % Center the lines around the point
            quiver(pc, pr, vc2/2, ur2/2, 0, 'y', 'showarrowhead', 'off', 'linewidth', 1.5) % Changed linewith from 2 to 1
            quiver(pc, pr, -vc2/2, -ur2/2, 0, 'y', 'showarrowhead', 'off', 'linewidth', 1.5) % Added to center the lines
        end
        ax = gca;
        im_out = getframe(ax);
        im_out = im_out.cdata;
        im_out = im_out(2:end, 2:end, :);
        str = '\nSave vector image';
        fprintf(str)
        lenghtstr = lenghtstr + length(str);
        imwrite(im_out, fullfile(output_dir, ['vectors_' file(1:end-4) '.tif']));
        close
        
        % angle heat map
        str = '\nCreate heat map';
        fprintf(str)
        lenghtstr = lenghtstr + length(str);
        figure('Position',[100 100 size(im,2)/3 size(im,1)/3]);
        if pading == true
            scale_up       = round(winsize/10);
            anglemat       = kron(anglemat, ones(scale_up)); % incerease resolution
        
            [numRows, numCols] = size(anglemat);
            heatmap_to_origianl_factor = numRows / maxwin;
            pad_left_top       = round(heatmap_to_origianl_factor * diff_min);
            pad_right_bottom   = round(heatmap_to_origianl_factor * diff_max);
        
            % Add pading
            zerosMatrix        = nan(pad_left_top, numCols);
            anglemat           = [zerosMatrix; anglemat];
        
            [numRows, numCols] = size(anglemat);
            zerosMatrix        = nan(numRows, pad_left_top);
            anglemat           = [zerosMatrix, anglemat];
        
            [numRows, numCols] = size(anglemat);
            zerosMatrix        = nan(numRows, pad_right_bottom);
            anglemat           = [anglemat, zerosMatrix];
            
            [numRows, numCols] = size(anglemat);
            zerosMatrix        = nan(pad_right_bottom, numCols);
            anglemat           = [anglemat; zerosMatrix];
        end
        imagesc(rad2deg(anglemat));
    
        % Colorbar
        hsv_nan = [[0, 0, 0];colormap('hsv')];
        set(gca,'visible','off')
        clim([0,180]);
        colormap(hsv_nan);
        set(gcf, 'InvertHardCopy', 'off');
        set(gcf, 'Color', [0 0 0]);
        set(gca, 'Color', 'none');
    
        colorbar;
    
        % Add black border
        %diff_min = parameters.diff_min;
        %anglemat_with_border = padarray(rad2deg(kron(anglemat, ones(125))), [diff_min diff_min], 0, 'both');
        %imagesc(anglemat_with_border);
    
        im_out = getframe(gcf);
        im_out = im_out.cdata;
        str = '\nSave heat map';
        fprintf(str)
        lenghtstr = lenghtstr + length(str);
        imwrite(im_out, fullfile(output_dir, ['angle_heatmap_colorbar' file(1:end-4) '.tif']));
    
        axis_out = getframe().cdata;
        imwrite(axis_out, fullfile(output_dir, ['angle_heatmap' file(1:end-4) '.tif']));
    
        close
        for s = 1:lenghtstr-4
            fprintf(1,'\b');
        end
    
    end
    
    % calculate order parameter
    av_ordermat = AFT_ordermat(anglemat, parameters);
    
    end