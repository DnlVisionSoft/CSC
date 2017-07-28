% CSC_Demo() 
% is a function of the CSC package that shows an example of how to use 
% functions for cell segmentation and cell counting. For mor details 
% about this algorithm, please refer to the following paper:
% 
% [1] Daniel Riccio, Nadia Brancati, Maria Frucci and Diego Gragnaniello, 
% "A New Unsupervised Approach for Segmenting and Counting Cells in 
% High-throughput Microscopy Image Sets", submitted to IEEE Journal of
% Biomedical and Health Informatics, pp. 1-10, 2017.
% 
%
%     Authors: Daniel Riccio, Maria Frucci, Nadia Brancati, Diego Gragnaniello
%     Matlab Implementation by Daniel Riccio, June 06, 2017.
%     Copyright (C) 2017 Daniel Riccio (daniel.riccio@unina.it)
%     also see https://github.com/DnlVisionSoft/CSC.git
%
%     This file is part of the Cell Segmentation and Counting based on Gray Lvel Clustering (CSC) package.
%
% 
%     CSC package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

function CSC_Demo()

      opts.patch_size = 512;
      opts.preprocessing = 0;
      opts.patch_step = 384;
      opts.comp_size_th = 200;
      opts.lambda = 0.5;
      opts.glc_th1 = 1.5;
      opts.glc_th2 = 0.2;
      opts.verbose = 1;

      A=imread('Cells.tif');

      [Results_seg] = CSC_CellSegmentation(A, opts)

      opts.clp_th = 3;

      [Results_cnt] = CSC_CellCounting(Results_seg.Mask, opts)