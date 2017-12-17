% CSC_BBBC006_Segmentation_and_Counting
% is a script of the CSC package that extracts and counts cells from
% from the images of the BBBC006 dataset.
% This script reads tiff images from the directory .\Images\BBBC006.
% Segmented images are written into the directoy .\Output\BBBC006\BinaryMasks.
% This script exploits functions from the CSC package.
%
% For mor details about this algorithm, please refer to the following paper:
%
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
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'BBBC006');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC006', 'BinaryMasks');
end

if(exist('.\Output','dir') && (exist('.\Output\BBBC006','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'BBBC006');
end

if(exist('.\Output\BBBC006','dir') && (exist('.\Output\BBBC006\BinaryMasks','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC006', 'BinaryMasks');
end


if(Automatic == 0)
    % segmentation
    opts.patch_size = 128;
    opts.preprocessing = 1;
    opts.patch_step = 64;
    opts.comp_size_th = 200;
    opts.lambda = 0.7;
    opts.glc_th1 = 1;
    opts.glc_th2 = 0.5;
    
    % counting
    opts.clp_th = 20;
    
    % output
    opts.verbose = 0;
    
else
    % segmentation
    opts.patch_size = 0;
    opts.preprocessing = 1;
    opts.patch_step = 0;
    opts.comp_size_th = 200;
    opts.lambda = 0.7;
    opts.glc_th1 = 1;
    opts.glc_th2 = 0;
    
    % counting
    opts.clp_th = 20;
    
    % output
    opts.verbose = 0;
end

if(exist('.\Images\BBBC006\Input2','dir'))
    
    hw = waitbar(0,'Segmenting...');
    
    d = dir('.\Images\BBBC006\Input2\*.tif')
     
    cellcnt = 1;
    
    indices = zeros(1,14);
    indices(2) = 1; % Rand Index
    indices(4:10) = 1; 

    
    for i=28:1:size(d,1)
        
        nome = sprintf('Images\\BBBC006\\Input2\\%s', d(i).name);
        
        A = double(imread(nome));
        
        
        [Results1] = CSC_CellSegmentation(A, opts);
        
        [Results2] = CSC_CellCounting(Results1.Mask, opts);
        
        nomesave = sprintf('.\\Output\\BBBC006\\BinaryMasks\\%s', d(i).name);
        imwrite(uint8(255*Results1.Mask), nomesave);
        
        nome2 = d(i).name;
        t = findstr(nome2, 's1');
        if(isempty(t))
            t=findstr(nome2, 's2');
        end
        nome2 = [nome2(1:t+1),'.png'];
        nome2 = sprintf('Images\\BBBC006\\Labels\\%s', nome2);
        nome2 = strrep(nome2, 'original', 'reference');
        Labels = double(imread(nome2));
        
        [Measures] = CSC_EvaluateSegmentation(Labels, Results2.Labels, indices);
        
        Statistics.nCellRef(i) = n;
        Statistics.nCellSeg(i) = m;
        
        Statistics.RI(i) = Measures.RI;
        Statistics.JI(i) = Measures.JI;
        Statistics.HM(i) = Measures.HM;
        Statistics.NSD(i) = Measures.NSD;
        Statistics.Split(i) = Measures.Split;
        Statistics.Merged(i) = Measures.Merged;
        Statistics.Added(i) = Measures.Added;
        Statistics.Missing(i) = Measures.Missing;
        
        mean(Statistics.RI(1:i))
        mean(Statistics.JI(1:i))
        mean(Statistics.HM(1:i))
        mean(Statistics.NSD(1:i))
        mean(Statistics.Split(1:i))
        mean(Statistics.Merged(1:i))
        mean(Statistics.Added(1:i))
        mean(Statistics.Missing(1:i))
        sprintf('--------------')
        
        waitbar(i/size(d,1),hw, sprintf('Segmenting (%d/%d)', i, size(d,1)));
        pause(0.1);
        
        save('Output\BBBC006\Statistics.mat', 'Statistics');
    end
    close(hw);
    
else
    sprintf('The input directory .\Images\BBBC006\Input2 does not exist')
end