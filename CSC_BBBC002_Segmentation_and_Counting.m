% CSC_BBBC002_Segmentation_and_Counting
% is a script of the CSC package that extracts and counts cells from 
% from the images of the BBBC002 dataset. 
% This script reads tiff images from the directory .\Images\BBBC002.
% Segmented images are written into the directoy .\Output\BBBC002\BinaryMasks.
% Annotated images are written into the directoy .\Output\BBBC002\Annotated.
% This script exploits functions from the CSC package. 
%
% For mor details about this algorithm, please refer to the following paper:
% 
%     A New Unsupervised Approach for Segmenting and Counting Cells in 
%     High-throughput Microscopy Image Sets", IEEE Journal of Biomedical 
%     and Health Informatics, pp. 1-12, 2018.
%     Authors: Daniel Riccio, Maria Frucci, Nadia Brancati, Diego Gragnaniello
%     Matlab Implementation by Daniel Riccio, June 06, 2017.
%     Copyright (C) 2017 Daniel Riccio (daniel.riccio@unina.it)
%     website: https://www.docenti.unina.it/daniel.riccio
%     also see https://github.com/DnlVisionSoft/CSC.git
%
%     This file is part of the Cell Segmentation and Counting based on 
%     Gray Lvel Clustering (CSC) package.
%
% 
%     The CSC package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY. 
%     Notice: permission is hereby granted, free of charge, to any person 
%     obtaining a copy of this software and associated documentation files 
%     (the "Software"), to deal in the Software without restriction, 
%     including without limitation the rights to use, copy, modify, merge, 
%     publish, distribute, sublicense, and/or sell copies of the Software, 
%     and to permit persons to whom the Software is furnished to do so, 
%     subject to the following conditions:
%
%     The above copyright notice and this permission notice shall be 
%     included in all copies or substantial portions of the Software.
%
%     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
%     CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
%     TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
%     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

% This variable select in which modality the method will run between:
%   0 -> Manual Tuning: parameters specifically tuned for this dataset are used
%   1 -> Automatic: general parameters are used
Automatic = 1;


if(exist('.\Output','dir')==0)
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\', 'Output');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'BBBC002');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC002', 'BinaryMasks');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC002', 'Annotated');
end

if(exist('.\Output','dir') && (exist('.\Output\BBBC002','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'BBBC002');
end

if(exist('.\Output\BBBC002','dir') && (exist('.\Output\BBBC002\BinaryMasks','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC002', 'BinaryMasks');
end

if(exist('.\Output\BBBC002','dir') && (exist('.\Output\BBBC002\Annotated','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC002', 'Annotated');
end
    
% We check if the input directory does exist
if(exist('.\Images\BBBC002','dir'))
   
    % we list all tiff files into the directory
    d = dir('Images\BBBC002\*.tif')

    % we open the ground truth file(txt) for cell counting
    fp = fopen('Images\BBBC002\BBBC002_v1_counts.txt', 'r');
    
    % we skip unuseful information
    fscanf(fp, '%s', 12);

    % we initialize the output 
    Statistics = [];


    if(Automatic == 0)
        % segmentation
        opts.patch_size = 256;
        opts.preprocessing = 0;
        opts.patch_step = 224;
        opts.comp_size_th = 10;
        opts.lambda = 0.7;
        opts.glc_th1 = 1;
        opts.glc_th2 = 0.15;
        
        % counting
        opts.clp_th = 20;

        % output
        opts.verbose = 0;
 
    else
        % segmentation
        opts.patch_size = 0;
        opts.preprocessing = 0;
        opts.patch_step = 0;
        opts.comp_size_th = 10;
        opts.lambda = 0.7;
        opts.glc_th1 = 1;
        opts.glc_th2 = 0;
        
        % counting
        opts.clp_th = 20;

        % output
        opts.verbose = 0;
    end

    hw = waitbar(0,'Segmenting and Counting...');
    
    gt1 = 0;
    gt2 = 0;
    cnt = 0;
    while(feof(fp)==0)
        
        nome = fscanf(fp, '%s', 1);
        gt1 = fscanf(fp, '%d', 1);
        gt2 = fscanf(fp, '%d', 1);
        
        if(length(nome)==0)
            break;
        end
        
        nome
        A = imread(sprintf('Images\\BBBC002\\%s.TIF', nome));

        [Results1] = CSC_CellSegmentation(A, opts);

        [Results] = CSC_CellCounting(Results1.Mask>0, opts);
        
        nomesave = sprintf('Output\\BBBC002\\BinaryMasks\\%s_segm.tif', nome);
        imwrite(uint8(255*Results1.Mask), nomesave);
        nomesave = sprintf('Output\\BBBC002\\Annotated\\%s_Annotated.tif', nome);
        imwrite(uint8(255*Results.AnnotatedMask), nomesave);
        
        cnt = cnt + 1;
        Statistics.ncells_seg(cnt) = Results.nCells;
        Statistics.ncells_gt1(cnt) = gt1;
        Statistics.ncells_gt2(cnt) = gt2;
        Statistics.ncells_err1 = sum((Statistics.ncells_seg-Statistics.ncells_gt1))/cnt
        Statistics.ncells_err2 = sum((Statistics.ncells_seg-Statistics.ncells_gt2))/cnt
        Statistics.ncells_err3 = sum((Statistics.ncells_gt1-Statistics.ncells_gt2))/cnt
        Statistics.ncells_std1 = std((Statistics.ncells_seg-Statistics.ncells_gt1))
        Statistics.ncells_std2 = std((Statistics.ncells_seg-Statistics.ncells_gt2))
        Statistics.ncells_std3 = std((Statistics.ncells_gt1-Statistics.ncells_gt2))
        
        waitbar(cnt/size(d,1),hw, sprintf('Segmenting and Counting (%d/%d)', cnt, size(d,1)));
        pause(0.1);
    end
    close(hw);
    fclose(fp);
    Statistics.MSE1 = mean((Statistics.ncells_seg-Statistics.ncells_gt1).^2);
    Statistics.MSE2 = mean((Statistics.ncells_seg-Statistics.ncells_gt2).^2);
    Statistics.MSE3 = mean((Statistics.ncells_gt1-Statistics.ncells_gt1).^2);
    
    save('Output\\BBBC002\\Statistics.mat', 'Statistics');
else
    sprintf('The input directory .\Images\BBBC002 does not exist')
end