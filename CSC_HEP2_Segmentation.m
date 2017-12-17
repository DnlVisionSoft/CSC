% CSC_HEP2_Segmentation
% is a script of the CSC package that extracts and counts cells from
% from the images of the HEP-2 dataset.
% This script reads tiff images from the directory .\Images\HEP-2.
% Segmented images are written into the directoy .\Output\HEP-2.
% This script exploits functions from the CSC package.
%
% For mor details about this algorithm, please refer to the following paper:
%
% [1] Daniel Riccio, Nadia Brancati, Maria Frucci and Diego Gragnaniello,
% "A New Unsupervised Approach for Segmenting and Counting Cells in
% High-throughput Microscopy Image Sets", submitted to IEEE Journal of
% Biomedical and Health Informatics, pp. 1-10, 2017
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
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'HEP-2');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\HEP-2', 'BinaryMasks');
end

if(exist('.\Output','dir') && (exist('.\Output\HEP-2','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'HEP-2');
end

if(exist('.\Output\HEP-2','dir') && (exist('.\Output\HEP-2\BinaryMasks','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\HEP-2\', 'BinaryMasks');
end


% We check if the input directory does exist
if(exist('Images\HEP-2','dir'))
    
    d = dir('Images\HEP-2\*.tif');
    
    
    if(Automatic == 0)
        % segmentation
        opts.patch_size = 256;
        opts.preprocessing = 1;
        opts.patch_step = 64;
        opts.comp_size_th = 200;
        opts.lambda = 0.7;
        opts.glc_th1 = 1;
        opts.glc_th2 = 0.3;
        opts.smooth_size = 3;
        
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
        opts.smooth_size = 3;
        
        % counting
        opts.clp_th = 20;

        % output
        opts.verbose = 0;
    end

    indices = zeros(1,14);
    indices(1:3) = 1;
    indices(12:14) = 1;
    
    hw = waitbar(0,'Segmenting...');

    Statistics = struct;
    Statistics.FMeasure = 0;
    Statistics.JD = 0;
    Statistics.RI = 0;
    Statistics.GCE = 0;
    Statistics.VOI = 0;
    Statistics.BDE = 0;
    
    cnt = 0;
    for i=1:2:size(d,1)
        nome = d(i).name
        A=imread(sprintf('Images\\HEP-2\\%s', nome));
        
        nomegt = d(i+1).name;
        GT=imread(sprintf('Images\\HEP-2\\%s', nomegt));
        
        [Results] = CSC_CellSegmentation(A, opts);
        
        [Measures] = CSC_EvaluateSegmentation(Results.Mask>0, GT>0, indices);
        cnt = cnt+1;
        Statistics.FMeasure(cnt) = Measures.Fmeasure;
        Statistics.JD(cnt) = Measures.JD;
        Statistics.RI(cnt) = Measures.RI;
        Statistics.GCE(cnt) = Measures.GCE;
        Statistics.VOI(cnt) = Measures.VOI;
        Statistics.BDE(cnt) = Measures.BDE;
        
        [mean(Statistics.FMeasure(1:cnt)), mean(Statistics.JD(1:cnt)), mean(Statistics.RI(1:cnt)), mean(Statistics.GCE(1:cnt)), mean(Statistics.VOI(1:cnt)), mean(Statistics.BDE(1:cnt))]
        
        nome = strrep(nome, '.tif', '_Mask.tif');
        imwrite(uint8(255*Results.Mask), sprintf('Output\\HEP-2\\BinaryMasks\\%s',nome));
        
        waitbar(i/size(d,1),hw, sprintf('Segmenting (%d/%d)', (i+1)/2, size(d,1)/2));
        pause(0.1);
    end
    close(hw);
else
    sprintf('The input directory .\Images\HEP-2 does not exist')
end

save('Output\\HEP-2\\Statistics.mat', 'Statistics');
