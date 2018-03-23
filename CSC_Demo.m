% CSC_Demo() 
% is a function of the CSC package that shows an example of how to use 
% functions for cell segmentation and cell counting. For mor details 
% about this algorithm, please refer to the following paper:
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
      