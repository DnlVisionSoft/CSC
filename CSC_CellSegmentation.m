% [Results] = CSC_CellSegmentation(A, opts)
% is a function of the CSC package that extracts cells from a grayscale image. 
% This function implements a segmentation process that is based on 
% gray level clustering. For mor details about this algorithm, please refer
% to the following paper:
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
%
%   Input
%   ------------------ 
%   A       A gray scale image of size HxW.
%   
%   opts    a structure of options with the following optional fields:
%
%       - patch_size [0, min(W,H)] (default 256) - it indicates the
%           size of the square patch that are independently segmented.
%           It corresponds to the parameter "n" in [1].
%
%       - preprocessing [0,1] (default 1) - it specifies if the 
%           preprocessing is applied (1) or not(0) to the input image.
%           It corresponds to the parameter "Preproc" in [1].
%
%       - patch_step (default patch_size/4) - it is the step of the sliding 
%           window used for extracting patches.
%           It corresponds to the parameter "beta" in [1].
%
%       - comp_size_th [0, min(W,H)] (default 200 px) - it is the 
%           threshold used to drop out small connected components corresponding 
%           to noise.
%           It corresponds to the parameter "epsilon" in [1].
%
%       - lambda (default 0.7) - it specifies the weight used to update the
%           threshold applied to binarize the current patch. 
%           It corresponds to the parameter "lambda" in [1].
%
%       - glc_th1 (default 1) - gray levels, whose distance is lower than 
%           this threshold are assigned to the same cluster by the GLC
%           algorithm.
%           It is fixed to 1.5 in [1].
%
%       - glc_th2 (default 0.5) - it is the threshold used to extract the
%           foreground in local patches
%           It corresponds to the parameter "alfa" in [1].
%
%       - verbose (default 0) - it specifies if messages and partial results 
%           must be shown(1) or not(0) during the segmentation process.
%
%
%   Output
%   ------------------ 
%   Results     structure of the segmented vessels with the following
%               fields:
%
%                   - QuantImage (size HxW) a matrix containing the
%                       quantized gray levels produced by the GLC
%                       algorithm.
%                   - Mask (size HxW) binary mask, where the
%                       foreground represents the extracted cells.
%           
%
%   Examples
%   --------
%
%   See CSC_Demo().
%
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

function [Results] = CSC_CellSegmentation(A, opts)

% we measure the execution time of the function in seconds
tic;

%% we parse the options and setup all variables

if(isfield(opts, 'patch_size') == 0)
    opts.patch_size = 256;
end

if(isfield(opts, 'preprocessing') == 0)
    opts.preprocessing = 1;
end

if(isfield(opts, 'patch_step') == 0)
    opts.patch_step = opts.patch_size/4;
end

if(isfield(opts, 'verbose') == 0)
    verbose = 0;
else
    verbose = (opts.verbose>0);
end

if(isfield(opts, 'lambda') == 0)
    lambda = 0.7;
else
    lambda = max(0,min(opts.lambda,1));
end

if(isfield(opts, 'comp_size_th') == 0)
    comp_size_th = 200;
else
    comp_size_th = opts.comp_size_th;
end

if(isfield(opts, 'glc_th1') == 0)
    opts.glc_th1 = 1.5;
else
    opts.glc_th1 = max(0,min(opts.glc_th1,1));
end

if(isfield(opts, 'glc_th2') == 0)
    opts.glc_th2 = 0.5;
else
    opts.glc_th2 = max(0,min(opts.glc_th2,1));
end

if(isfield(opts, 'smooth_size') == 0)
    opts.smooth_size = 0;
else
    opts.smooth_size;
end

% we set the output variable
Results = [];

% we check if the image format is compatible with next processing steps
% gray leve image
if(size(A,3)==3)
    A = rgb2gray(A);
end

% some images present a lot of noisy pixels with gray tones less then 2
A = double(A).*(A>2);

QuantImage = zeros(size(A));
Mask = zeros(size(A));
    

% we apply global preprocessing
if(opts.preprocessing)
    A=128*(1+quasi_sigmoidal(double(A), mean(A(:))));
    [A] = NormCell(uint8(A));
end
A = double(A);

if(opts.patch_size == 0)
%     S = qtdecomp(double(A(1:512, 1:512))/255,.1, [64 512]);
%     s=[8,16,32,64,128,256,512];
%     h=zeros(1,7);
%     for i=1:7
%         t=find(S==s(i));
%         h(i)=length(t);
%     end
%     q=max(find(h==max(h)));
%     q = s(q);
%     
%     opts.patch_size = 4*q;
%     opts.patch_step = round(3*q);
%     [opts.patch_size opts.patch_step]
    [ws, n, beta] = find_patch(A, 256, opts);
    opts.patch_size = n;
    opts.patch_step = beta;
end
        
stp = opts.patch_size;
stp2 = opts.patch_step;

M=0;
hw = waitbar(0,'Segmenting square patches...');

% we split the input image in square patches and process each of them
% independently according to a zig-zag path
jp = [1:stp2:size(A,2)];
jp = fliplr(jp);

for i=1:stp2:size(A,1)
    jp = fliplr(jp);

    for j=1:length(jp)

        xstart = i;
        ystart = jp(j);
        xstop  = min(i+stp, size(A,1));
        ystop  = min(jp(j)+stp, size(A,2));
        
        c = A(xstart:xstop, ystart:ystop);
                  
        [Mask1c, Qc] = ClusterizeGrayLevels(c, M, opts);
        
        S = Qc(Qc>1);
        if(isempty(S)==0)
            med=median(S);
            mad=median(abs(S-med));
            if(M==0)
                M = med+mad;
            else
                M = lambda*M + (1-lambda)*med;
            end
        end
        
        Mask(xstart:xstop, ystart:ystop) = (Mask(xstart:xstop, ystart:ystop) + Mask1c)>0;
        QuantImage(xstart:xstop, ystart:ystop) = max(QuantImage(xstart:xstop, ystart:ystop), Qc); 
            
    end
    ptmp = ((i-1)*size(A,2))/(size(A,1)*size(A,2))+0.01;
    waitbar(ptmp,hw,'Segmenting square patches...');
    
    if(verbose)
        figure(1);
        subplot(2,2,1);
        imshow(uint8(c));
        subplot(2,2,2);
        imshow(Mask1c,[]);
        subplot(2,2,3);
        imshow(QuantImage,[]);
        subplot(2,2,4);
        imshow(Mask,[]);
        drawnow;
    end
end

waitbar(1,hw,'Segmenting square patches...');
pause(0.1);
close(hw);

% we fill the holes in the foreground
Mask1 = (imfill(Mask, 'holes'))>0;

% we delete all connected components with size smaller than comp_size_th
% pixls in the foreground
[L,n] = bwlabel(Mask1);
Mask = zeros(size(Mask1));
for i=1:n
    c = (L==i);
    
    if(sum(sum(c))>comp_size_th)
        Mask = Mask + c;
    end
end
    
Results.Mask = Mask;
Results.QuantImage = QuantImage;
Results.Time = toc;




        
%% Side Functions

% NORMACELL
function [B] = NormCell(A)

A = double(A);

m = mean(A(:));
T = A>m;
Q = T.*double(A);
q = Q(Q(:)>0);
mn = mean(q);

f = max(0, 1-(A-mn)/(mn-m));
vmax = (255/max(f(:)));
B = floor(abs(255-vmax*f));

% CLUSTERIZEGRAYLEVELS
function [Mask, Q] = ClusterizeGrayLevels(A, thbin, opts)

th = opts.glc_th1; 

B = double(A); 

if(min(size(B,1), size(B,2))>8)
    C = adapthisteq(uint8(B));
else
    C = B;
end

if(opts.smooth_size > 1)
    C = medfilt2(uint8(C), [opts.smooth_size opts.smooth_size]);
end

h = histc(C(:), 0:255);
P = h/sum(h);

[SP] = ComputeSparseness(C);


T=zeros(size(B));
for i=1:size(B,1)
    for j=1:size(B,2)
        T(i,j)=SP(1+C(i,j));
    end
end

Pixels = double(unique(C(:)));

n = length(Pixels);

delta = 1; 
if(n > 1)
    M = zeros(n, n);
    for i=1:n
        for j=1:n
            M(i,j)=abs(P(Pixels(i)+1)*Pixels(i)-P(Pixels(j)+1)*Pixels(j))+delta*log(abs(Pixels(i)-Pixels(j))+1) + min(SP(i), SP(j))/(max(SP(i),SP(j))+1);
        end
    end
end

Clusters = zeros(1,n);
Assigned = zeros(1,n);
dist     = zeros(1,n);
nclusters = 0;
for i=1:n
    dist(:)=1000000;
    for j=1:nclusters
        k=1;
        while(Clusters(j,k)>0)
          dist(j) = min(dist(j), M(i,Clusters(j,k)));
          k = k+1;
        end
    end
    dminpos = find(dist==min(dist),1);
    dmin = min(dist);
    
    if(dmin < th)
        pos = find(Clusters(dminpos,:)==0,1);
        Clusters(dminpos, pos) = i;
        Assigned(i)=1;
    else
        nclusters = nclusters + 1;
        Clusters(nclusters, 1) = i;
        Assigned(i) = 1;
    end
end

Q = zeros(size(B));
for i=1:n
    c = (C==(Pixels(i)));
    [x,y] = find(Clusters==i);
    Q = Q + c*x;
end

S = Q(Q>1);
med=median(S);
mad=median(abs(S-med));
if(thbin==0)
    thbin = med+2*mad;
end

if(opts.glc_th2==0)
    entQ = entropy(uint8(Q));
    entA = entropy(uint8(A));

%       diff_ent = entQ - entA
%       if(diff_ent>1)
%         THfinal = log2(diff_ent)/4;
%       else
%         THfinal = 0.05;
%       end
%       THfinal
    THfinal = entQ/24;
%     THfinal = min(max(0,1-entA/entQ), 0.4); %entropy(uint8(Q))/24;
else
    THfinal = opts.glc_th2;
end

% Mask = Q>(thbin+0.5*mad);
Mask = Q>thbin;

for th = thbin:max(Q(:))
    if(sum(thbin:th)/sum(thbin:max(Q(:))) > THfinal)
        Mask = (Q>th);
        break;
    end
end


% COMPUTESPARSENESS
function [SP] = ComputeSparseness(A)

SP = zeros(1,256);

for i=0:255
    [x, y] = find(A==i);
    
    Xc = mean(x);
    Yc = mean(y);
    
    t  = sqrt(mean((x-Xc).^2 + (y-Yc).^2));
    if(isnan(t)==0)
        SP(i+1) = t;
    end
end



% QUASI_SIGMOIDAL
function val = quasi_sigmoidal(x, xmax)

a = 2 + sqrt(3);
b = 7 - 4*sqrt(3);

val = (1 - b.^(x./xmax))./(a*b.^(x./xmax) + 1);


function [ws, n, beta] = find_patch(A, sz, opts)

entp = 0;
sz2 = floor(sz/2);

Ptc = A;

for i=1:sz2:size(A,1)-sz
    for j=1:sz2:size(A,2)-sz
        p = A(i:i+sz-1, j:j+sz-1);
        if(entropy(uint8(p))> entp)
            Ptc = p;
        end
    end
end

[Mask1c, Q] = ClusterizeGrayLevels(Ptc, 0, opts);
entP = entropy(uint8(Ptc));
entQ = entropy(uint8(Q));

diff_ent = entQ - entP;
if(diff_ent>1)
    ws = log2(diff_ent)/4;
else
    ws = 0;
end

% ws = min(ws, 0.5);
n = round((1-1.5*ws)*512);
n = max(n,128);
beta = round(3/4*n);

