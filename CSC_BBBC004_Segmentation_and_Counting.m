% CSC_BBBC004_Segmentation_and_Counting
% is a script of the CSC package that extracts and counts cells from
% from the images of the BBBC004 dataset.
% This script reads tiff images from the directory .\Images\BBBC004.
% Segmented images are written into the directoy .\Output\BBBC004\BinaryMasks.
% Annotated images are written into the directoy .\Output\BBBC004\Annotated.
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
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'BBBC004');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC004', 'BinaryMasks');
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC004', 'Annotated');
end

if(exist('.\Output','dir') && (exist('.\Output\BBBC004','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output', 'BBBC004');
end

if(exist('.\Output\BBBC004','dir') && (exist('.\Output\BBBC004\BinaryMasks','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC004', 'BinaryMasks');
end

if(exist('.\Output\BBBC004','dir') && (exist('.\Output\BBBC004\Annotated','dir')==0))
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir('.\Output\BBBC004', 'Annotated');
end

if(Automatic == 0)
    % segmentation
    opts.patch_size = 512;
    opts.preprocessing = 0;
    opts.patch_step = 384;
    opts.comp_size_th = 200;
    opts.lambda = 0.5;
    opts.glc_th1 = 1;
    opts.glc_th2 = 0.2;
    
    % counting
    opts.clp_th = 3;
    
    % output
    opts.verbose = 0;
    
else
    % segmentation
    opts.patch_size = 0;
    opts.preprocessing = 0;
    opts.patch_step = 0;
    opts.comp_size_th = 200;
    opts.lambda = 0.7;
    opts.glc_th1 = 1;
    opts.glc_th2 = 0;
    
    % counting
    opts.clp_th = 3;
    
    % output
    opts.verbose = 0;
end

% We process all images according to the degree of overlapping
indices = zeros(1,14);
indices(1) = 1; % we compute only the F-Measure

if(exist('.\Images\BBBC004','dir'))
    % Overlapping 0
    d = dir('Images\BBBC004\synthetic_000_images\*.tif')
    
    hw = waitbar(0,'Segmenting and Counting Overlap(0%%)...');
    Errore = zeros(size(d,1), 5);
    Counts = zeros(size(d,1), 5);
    Accuracy = zeros(size(d,1), 5);
    Sensitivity = zeros(size(d,1), 5);
    Specificity = zeros(size(d,1), 5);
    Fmeasure = zeros(size(d,1), 5);
    MSE = zeros(1, 5);
    for i=1:size(d,1)
        nome = sprintf('Images\\BBBC004\\synthetic_000_images\\%s', d(i).name);
        A = imread(nome);
        
        [Results1] = CSC_CellSegmentation(A, opts);
        [Results2] = CSC_CellCounting(Results1.Mask, opts);
        
        nomesave = sprintf('Output\\BBBC004\\BinaryMasks\\%s', strrep(d(i).name, '.tif', '_overlap0_segm.tif'));
        imwrite(uint8(255*Results1.Mask), nomesave);
        nomesave = sprintf('Output\\BBBC004\\Annotated\\%s', strrep(d(i).name, '.tif', '_overlap0_annoated.tif'));
        imwrite(uint8(255*Results2.AnnotatedMask), nomesave);
        
        Errore(i,1) = (Results2.nCells-300)/300*100;
        Counts(i,1) = Results2.nCells;
        Counts(:,1)
        nome = strrep(nome, 'synthetic_000_images', 'synthetic_000_foreground');
        nome = strrep(nome, 'GRAY', '');
        GT = rgb2gray(imread(nome))>0;
        
        [Measures] = CSC_EvaluateSegmentation(Results1.Mask, GT, indices);
        Accuracy(i,1) = Measures.Accuracy;
        Sensitivity(i,1) = Measures.Sensitivity;
        Specificity(i,1) = Measures.Specificity;
        Fmeasure(i,1) = Measures.Fmeasure;
        mean(Errore(1:i,1))
        std(Errore(1:i,1))
        mean(Accuracy(1:i,1))
        mean(Sensitivity(1:i,1))
        mean(Specificity(1:i,1))
        mean(Fmeasure(1:i,1))
        waitbar(i/size(d,1),hw, sprintf('Segmenting and Counting (0%%) (%d/%d)', i, size(d,1)));
        pause(0.1);
    end
    mean(Errore)
    std(Errore)
    mean(Accuracy)
    mean(Sensitivity)
    mean(Specificity)
    mean(Fmeasure)
    MSE(1) = mean((Counts(:,1)-300).^2);
    
    
    % Overlapping 15
    d = dir('Images\\BBBC004\\synthetic_015_images\\*.tif')
    for i=1:size(d,1)
        nome = sprintf('Images\\BBBC004\\synthetic_015_images\\%s', d(i).name);
        A = imread(nome);
        [Results1] = CSC_CellSegmentation(A, opts);
        [Results2] = CSC_CellCounting(Results1.Mask, opts);
        
        nomesave = sprintf('Output\\BBBC004\\BinaryMasks\\%s', strrep(d(i).name, '.tif', '_overlap15_segm.tif'));
        imwrite(uint8(255*Results1.Mask), nomesave);
        nomesave = sprintf('Output\\BBBC004\\Annotated\\%s', strrep(d(i).name, '.tif', '_overlap15_annoated.tif'));
        imwrite(uint8(255*Results2.AnnotatedMask), nomesave);
        
        Errore(i,2) = (Results2.nCells-300)/300*100;
        Counts(i,2) = Results2.nCells;
        Counts(:,1:2)
        mean(Errore(1:i,:))
        nome = strrep(nome, 'synthetic_015_images', 'synthetic_015_foreground');
        nome = strrep(nome, 'GRAY', '');
        GT = rgb2gray(imread(nome))>0;
        
        [Measures] = CSC_EvaluateSegmentation(Results1.Mask, GT, indices);
        Accuracy(i,2) = Measures.Accuracy;
        Sensitivity(i,2) = Measures.Sensitivity;
        Specificity(i,2) = Measures.Specificity;
        Fmeasure(i,2) = Measures.Fmeasure;
        mean(Errore(1:i,2))
        std(Errore(1:i,2))
        mean(Accuracy(1:i,2))
        mean(Sensitivity(1:i,2))
        mean(Specificity(1:i,2))
        mean(Fmeasure(1:i,2))
        waitbar(i/size(d,1),hw, sprintf('Segmenting and Counting (15%%) (%d/%d)', i, size(d,1)));
        pause(0.1);
    end
    mean(Errore)
    std(Errore)
    mean(Accuracy)
    mean(Sensitivity)
    mean(Specificity)
    mean(Fmeasure)
    MSE(2) = mean((Counts(:,1)-300).^2);
    
    % Overlapping 30
    d = dir('Images\\BBBC004\\synthetic_030_images\\*.tif')
    for i=1:size(d,1)
        nome = sprintf('Images\\BBBC004\\synthetic_030_images\\%s', d(i).name);
        A = imread(nome);
        [Results1] = CSC_CellSegmentation(A, opts);
        [Results2] = CSC_CellCounting(Results1.Mask, opts);
        
        nomesave = sprintf('Output\\BBBC004\\BinaryMasks\\%s', strrep(d(i).name, '.tif', '_overlap30_segm.tif'));
        imwrite(uint8(255*Results1.Mask), nomesave);
        nomesave = sprintf('Output\\BBBC004\\Annotated\\%s', strrep(d(i).name, '.tif', '_overlap30_annoated.tif'));
        imwrite(uint8(255*Results2.AnnotatedMask), nomesave);
        
        Errore(i,3) = (Results2.nCells-300)/300*100;
        Counts(i,3) = Results2.nCells;
        Counts(:,1:3)
        nome = strrep(nome, 'synthetic_030_images', 'synthetic_030_foreground');
        nome = strrep(nome, 'GRAY', '');
        GT = rgb2gray(imread(nome))>0;
        
        [Measures] = CSC_EvaluateSegmentation(Results1.Mask, GT, indices);
        Accuracy(i,3) = Measures.Accuracy;
        Sensitivity(i,3) = Measures.Sensitivity;
        Specificity(i,3) = Measures.Specificity;
        Fmeasure(i,3) = Measures.Fmeasure;
        mean(Errore(1:i,3))
        std(Errore(1:i,3))
        mean(Accuracy(1:i,3))
        mean(Sensitivity(1:i,3))
        mean(Specificity(1:i,3))
        mean(Fmeasure(1:i,3))
        waitbar(i/size(d,1),hw, sprintf('Segmenting and Counting (30%%) (%d/%d)', i, size(d,1)));
        pause(0.1);
    end
    mean(Errore)
    std(Errore)
    mean(Accuracy)
    mean(Sensitivity)
    mean(Specificity)
    mean(Fmeasure)
    MSE(3) = mean((Counts(:,1)-300).^2);
    
    % Overlapping 45
    d = dir('Images\\BBBC004\\synthetic_045_images\\*.tif')
    for i=1:size(d,1)
        nome = sprintf('Images\\BBBC004\\synthetic_045_images\\%s', d(i).name);
        A = imread(nome);
        [Results1] = CSC_CellSegmentation(A, opts);
        [Results2] = CSC_CellCounting(Results1.Mask, opts);
        
        nomesave = sprintf('Output\\BBBC004\\BinaryMasks\\%s', strrep(d(i).name, '.tif', '_overlap45_segm.tif'));
        imwrite(uint8(255*Results1.Mask), nomesave);
        nomesave = sprintf('Output\\BBBC004\\Annotated\\%s', strrep(d(i).name, '.tif', '_overlap45_annoated.tif'));
        imwrite(uint8(255*Results2.AnnotatedMask), nomesave);
        
        Errore(i,4) = (Results2.nCells-300)/300*100;
        Counts(i,4) = Results2.nCells;
        Counts(:,1:4)
        nome = strrep(nome, 'synthetic_045_images', 'synthetic_045_foreground');
        nome = strrep(nome, 'GRAY', '');
        GT = rgb2gray(imread(nome))>0;
        
        [Measures] = CSC_EvaluateSegmentation(Results1.Mask, GT, indices);
        Accuracy(i,4) = Measures.Accuracy;
        Sensitivity(i,4) = Measures.Sensitivity;
        Specificity(i,4) = Measures.Specificity;
        Fmeasure(i,4) = Measures.Fmeasure;
 
        mean(Errore(1:i,4))
        std(Errore(1:i,4))
        mean(Accuracy(1:i,4))
        mean(Sensitivity(1:i,4))
        mean(Specificity(1:i,4))
        mean(Fmeasure(1:i,4))
        waitbar(i/size(d,1),hw, sprintf('Segmenting and Counting (45%%) (%d/%d)', i, size(d,1)));
        pause(0.1);
    end
    mean(Errore)
    std(Errore)
    mean(Accuracy)
    mean(Sensitivity)
    mean(Specificity)
    mean(Fmeasure)
    MSE(4) = mean((Counts(:,1)-300).^2);
    
    % Overlapping 60
    d = dir('Images\\BBBC004\\synthetic_060_images\\*.tif')
    for i=1:size(d,1)
        nome = sprintf('Images\\BBBC004\\synthetic_060_images\\%s', d(i).name);
        A = imread(nome);
        [Results1] = CSC_CellSegmentation(A, opts);
        [Results2] = CSC_CellCounting(Results1.Mask, opts);
        
        nomesave = sprintf('Output\\BBBC004\\BinaryMasks\\%s', strrep(d(i).name, '.tif', '_overlap60_segm.tif'));
        imwrite(uint8(255*Results1.Mask), nomesave);
        nomesave = sprintf('Output\\BBBC004\\Annotated\\%s', strrep(d(i).name, '.tif', '_overlap60_annoated.tif'));
        imwrite(uint8(255*Results2.AnnotatedMask), nomesave);
        
        Errore(i,5) = (Results2.nCells-300)/300*100;
        Counts(i,5) = Results2.nCells;
        Counts(:,1:5)
        mean(Errore(1:i,:))
        nome = strrep(nome, 'synthetic_060_images', 'synthetic_060_foreground');
        nome = strrep(nome, 'GRAY', '');
        GT = rgb2gray(imread(nome))>0;
        
        [Measures] = CSC_EvaluateSegmentation(Results1.Mask, GT, indices);
        Accuracy(i,5) = Measures.Accuracy;
        Sensitivity(i,5) = Measures.Sensitivity;
        Specificity(i,5) = Measures.Specificity;
        Fmeasure(i,5) = Measures.Fmeasure;
        
        mean(Errore(1:i,5))
        std(Errore(1:i,5))
        mean(Accuracy(1:i,5))
        mean(Sensitivity(1:i,5))
        mean(Specificity(1:i,5))
        mean(Fmeasure(1:i,5))
        waitbar(i/size(d,1),hw, sprintf('Segmenting and Counting (60%%) (%d/%d)', i, size(d,1)));
        pause(0.1);
    end
    mean(Errore)
    std(Errore)
    mean(Accuracy)
    mean(Sensitivity)
    mean(Specificity)
    mean(Fmeasure)
    MSE(5) = mean((Counts(:,1)-300).^2);
    
    Statistics = struct;
    Statistics.Errore = Errore;
    Statistics.Counts = Counts;
    Statistics.Accuracy = Accuracy;
    Statistics.Sensitivity = Sensitivity;
    Statistics.Specificity = Specificity;
    Statistics.Fmeasure = Fmeasure;
    Statistics.AvgErrore = mean(Errore);
    Statistics.StdErrore = std(Errore);
    Statistics.AvgAccuracy = mean(Accuracy);
    Statistics.AvgSensitivity = mean(Sensitivity);
    Statistics.AvgSpecificity = mean(Specificity);
    Statistics.AvgFmeasure = mean(Fmeasure);
    Statistics.MSE = MSE;
    
    save('Output\\BBBC004\\Statistics.mat', 'Statistics');
else
    sprintf('The input directory .\Images\BBBC004 does not exist')
end


