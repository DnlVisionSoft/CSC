% CSC_BBBC008_Segmentation
% is a script of the CSC package that extracts and counts cells from
% from the images of the BBBC008 dataset.
% This script reads tiff images from the directory .\Images\BBBC008.
% Segmented images are written into the directoy .\Output\BBBC008\BinaryMasks.
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
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'BBBC008');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC008', 'BinaryMasks');
end

if(exist('.\Output','dir') && (exist('.\Output\BBBC008','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'BBBC008');
end

if(exist('.\Output\BBBC008','dir') && (exist('.\Output\BBBC008\BinaryMasks','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC008', 'BinaryMasks');
end


if(Automatic == 0)
    % segmentation
    opts.patch_size = 512;
    opts.preprocessing = 0;
    opts.patch_step = 512;
    opts.comp_size_th = 10;
    opts.lambda = 0.5;
    opts.glc_th1 = 1;
    opts.glc_th2 = 0.05;
    
    % counting
    opts.clp_th = 3;
    
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
    opts.glc_th2 = 0.0;
    
    % counting
    opts.clp_th = 3;
    
    % output
    opts.verbose = 0;
end

indices = zeros(1,14);
indices(1) = 1;
indices(11) = 1;

if(exist('.\Images\BBBC008','dir'))
    
    d = dir('.\Images\BBBC008\Channel1_Input\*.tif')
    
    Accuracy = 0;
    Specificity = 0;
    Sensitivity = 0;
    Precision = 0;
    Fmeasure = 0;
    RER = 0;
    
    hw = waitbar(0,'Segmenting...');
    fp = fopen('.\Images\BBBC008\names.csv', 'w');
    for i=1:size(d,1)
        
        nome = d(i).name;
        
        A = imread(sprintf('.\\Images\\BBBC008\\Channel1_Input\\%s', nome));

        [Results] = CSC_CellSegmentation(A, opts);
        
        nomesave = sprintf('.\\Output\\BBBC008\\BinaryMasks\\%s', nome);
        imwrite(uint8(255*Results.Mask), nomesave);

        fprintf(fp, '%s\n', nome);
        nome = sprintf('.\\Images\\BBBC008\\Channel1_GT\\%s', d(i).name);
        GT = double(imread(nome));
        
        [Measures] = CSC_EvaluateSegmentation(Results.Mask, GT, indices);

        Accuracy = Accuracy + Measures.Accuracy;
        Specificity = Specificity + Measures.Specificity;
        Sensitivity = Sensitivity + Measures.Sensitivity;
        Precision = Precision + Measures.Precision;
        Fmeasure = Fmeasure + Measures.Fmeasure;
        RER = RER + Measures.RER;
        
        [Accuracy/i Specificity/i Sensitivity/i Precision/i, Fmeasure/i RER/i]
        
        waitbar(i/size(d,1),hw, sprintf('Segmenting (%d/%d)', i, size(d,1)));
        pause(0.1);
    end
    fclose(fp);
    close(hw);
    
    Statistics = struct;
    Statistics.Accuracy = Accuracy/i;
    Statistics.Specificity = Specificity/i;
    Statistics.Sensitivity = Sensitivity/i;
    Statistics.Precision = Precision/i; 
    Statistics.Fmeasure = Fmeasure/i; 
    Statistics.RER = RER/i;
    
    save('Output\\BBBC008\\Statistics.mat', 'Statistics');
else
    sprintf('The input directory .\Images\BBBC008 does not exist')
end
