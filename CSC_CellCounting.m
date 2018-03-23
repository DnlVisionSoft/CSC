% [Results] = CSC_CellCounting(Mask, opts) 
% is a function of the CSC package that counts cells in a binary image, 
% where the foreground represents segmented cells. 
% This function implements a counting process that is based on distance 
% transform and curvature approximation. For mor details about this 
% algorithm, please refer to the following paper:
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
%   Mask    A binary image of size HxW, where foreground pixels are 1 and
%           background pixels are 0
%   
%   opts    a structure of options with the following optional fields:
%
%       - max_num_seeds (default 1000) - the max number of seeds allowed
%           that are produced by the distance transform.
%
%           It corresponds to the parameter "delta" in [1].
%       - clp_th (default 1) - seeds, whose distance is lower than 
%           this threshold are assigned to the same cluster by the
%           incremental clustering algorithm.
%           It corresponds to the parameter "delta" in [1].
%
%       - verbose (default 0) - it specifies if messages and partial results 
%           must be shown(1) or not(0) during the counting process.
%
%
%   Output
%   ------------------ 
%   Results     structure of the segmented vessels with the following
%               fields:
%
%                   - Labels (size HxW) a matrix where each cell is 
%                       represented with a different label.
%                   - Annotated Mask (size HxW) color image where cells in
%                       the foreground mask are annotated by a red ellipse.
%                   - AvgCellSize the average size of a cell
%                   - StdCellSize the standard varition of cell size
%                   - nEstimedCells the number of seeds that correspond to
%                        the estimed number of cells
%                   - nCells (scalar) the number of cells found in the
%                       input image.
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

function [Results] = CSC_CellCounting(Mask, opts)

% we measure the execution time of the function in seconds
tic;

% we expand the original mask to deal with cells on the image borders
s = size(Mask);
stp = 50;

Mask_tmp = zeros(size(Mask,1)+2*stp, size(Mask,2)+2*stp);
Mask_tmp(stp:size(Mask,1)+stp-1, stp:size(Mask,2)+stp-1) = Mask;

Mask = Mask_tmp;

% we set up all images produced by this function as intermediate and/or
% finl result
C = zeros(size(Mask));
Cmp = zeros(size(Mask));
CellsImage = zeros(size(Mask));
TempImage = zeros(size(Mask));

% we apply the distance transform and locate regional maxia
U=bwulterode(Mask);
% [U] = CSC_DetectSeeds(Mask, opts);

% we get the spatial coordinate of the regional maxima
[X,Y]=find(U>0);

if(isfield(opts, 'max_num_seeds') == 0)
    opts.max_num_seeds = 1000;
end

if(isfield(opts, 'verbose') == 0)
    opts.verbose = 0;
end

if(isfield(opts, 'clp_th') == 0)
    opts.clp_th = 5;
end


% if the number of initial seeds exceed the maximum number of 1000 an error
% message is shown and empty images are returned by the function
if(length(X)>opts.max_num_seeds)
    Results.Labels = CellsImage;
    Results.AnnotatedMask = TempImage;
    Results.nCells = 0;
    sprintf('Too many Seeds found: %d', length(X))
    return;
end

if(opts.verbose)
    [Xf, Yf] = find(U>0);
    figure(3);
    imshow(Mask, []);
    hold on;
    plot(Yf, Xf, 'r^', 'MarkerSize', 8, 'MarkerEdgeColor',  'k', 'MarkerFaceColor', 'r');  
    
    TempIm = frame2im(getframe(gcf)); 
    Results.DistTransSeedsIm = TempIm;
end
    
% we apply the incremental clustering to the original seeds to collapse
% seeds whose distance is lower than the threshold opts.clp_th
P=[X,Y];
[Seeds, Clusters, nclusters] = CSC_ClusterizePoints(P, opts.clp_th);

SingleSeeds = Seeds;

%we generate a black image with white pixels only in correspondence of the
%seeds
Q=zeros(size(Mask));
for i=1:size(Seeds,1)
    Q(Seeds(i,1), Seeds(i,2))=i;
end

if(opts.verbose)
    [Xf, Yf] = find(Q>0);
    figure(3);
    hold on;
    plot(Yf, Xf, 'go', 'MarkerSize', 5, 'MarkerEdgeColor',  'k', 'MarkerFaceColor', 'g');
    
    TempIm = frame2im(getframe(gcf)); 
    Results.ClusteredSeedsIm = TempIm;
end


Q1 = Q;
 

% we first process connected components including only one seed, that means
% this component surely corresponds to a single cell. Moreover, we generate
% a new edge image, where boundry points are labeled with the label of the
% assigned seed
crem = bwmorph(Mask, 'remove');
CellSize = 0;
t=[];
[L,n] = bwlabel(Mask>0);
for k=1:n
    labs = Q.*(L==k);
    nlabs = sum(labs(:)>0);

    if(nlabs==1)
        clab = unique(labs(labs>0));
        
        CellSize(k) = sum(sum(crem.*(L==k)));
        t = [t, clab];
        CellsImage = CellsImage + clab.*(L==k);
    end
    
    if(opts.verbose)
        figure(2);
        subplot(1,2,1);
        imshow(128*(Q>0)+(128*double(bwmorph(CellsImage>0, 'remove'))+32*double(Mask>0)), []);
        title(sprintf('Cells found: %d wrt estimed: %d', length(t), sum(Q(:)>0)));
        subplot(1,2,2);
        imshow(label2rgb(CellsImage));
        title('Cells body represented by different colours');
        drawnow;
    end
end 


% we compute the average and standard deviation over the size of single
% cells
AvgCellSize = mean(CellSize);
StdCellSize = std(CellSize);

Results.AvgCellSize = AvgCellSize;
Results.StdCellSize = StdCellSize;

size(Seeds)
nst = max(size(Seeds));
MaskClust = Mask.*(1-(CellsImage>0));
Q = Q.*(1-(CellsImage>0));
[U] = CSC_DetectSeeds(MaskClust, opts, 0);
[Xd,Yd]=find(U>0);
for i=1:length(Xd)
    d = sqrt((Xd(i)-Seeds(:,1)).^2 + (Yd(i)-Seeds(:,2)).^2);
    if(min(d)>(AvgCellSize+StdCellSize/2)/(2*pi))
        Seeds = [Seeds; [Xd(i), Yd(i)]];
        Q(Xd(i), Yd(i))=size(Seeds,1);
    end
end
size(Seeds)

if(opts.verbose)
    [Xf, Yf] = find(Q>=nst);
    figure(3);
    hold on;
    plot(Yf, Xf, 'bs', 'MarkerSize', 5, 'MarkerEdgeColor',  'k', 'MarkerFaceColor', 'b');
    
    TempIm = frame2im(getframe(gcf)); 
    Results.ClusteredSeedsIm = TempIm;
end


% Next, we process connected components including more than one seed, which
% correponds to cluster of cells.
[L,n] = bwlabel(MaskClust>0);
d = zeros(1,size(Seeds,1));
for k=1:n
    nlabs = sum(sum((Q.*(L==k))>0));
    if(nlabs>1)

        labs = Q.*(L==k);
        labs = setdiff(unique(labs(:)), 0);
        p = labs;
        t = [t, p'];
   
        % boundary points are assigned to the closest seed in the connected
        % component. Only seeds that have passed the test on the number of
        % assigned points are considered in this phase
%         if(length(p)>0)
            [X,Y] = find((L==k)>0);
            d(:)=0;
            for i=1:length(X)
                d = abs(X(i)-Seeds(p,1))+abs(Y(i)-Seeds(p,2));
                lab = find(d==min(d), 1 );
                CellsImage(X(i),Y(i)) = p(lab);
            end
%         end
    end
    
    if(opts.verbose)
        figure(2);
        subplot(1,2,1);
        imshow(128*((Q+Q1)>0)+(128*double(bwmorph(CellsImage>0, 'remove'))+32*double(Mask>0)), []);
        title(sprintf('Cells found: %d wrt estimed: %d', length(unique(CellsImage(:)))-1, sum(Q(:)>0)));
        subplot(1,2,2);
        imshow(label2rgb(CellsImage));
        title('Cells body represented by different colours');
        drawnow;
    end
end 

% we only consider seeds that passed the previous test regarding the number
% of assigned boundary points
Seeds = Seeds(t,:);

Results.nEstimedCells = size(Seeds,1);

% we generate again a black image with white points corresponding to
% detected seeds
Q=zeros(size(Mask));
for i=1:size(Seeds,1)
    Q(Seeds(i,1), Seeds(i,2))=1;
end

TempImage=zeros(size(Mask,1), size(Mask,2), 3);
for k=1:max(CellsImage(:))
    c = (CellsImage == k);
    if(sum(c(:))>0)
        TempImage(:,:,1) = TempImage(:,:,1)+255*double(bwmorph(c>0, 'remove'));
    end
end

TempImage(:,:,2) = ((Q+Q1)>0);
TempImage(:,:,3) = Mask.*(1-TempImage(:,:,1)).*(1-TempImage(:,:,2));


TempImage = TempImage(stp:s(1)+stp-1, stp:s(2)+stp-1, :);
CellsImage = CellsImage(stp:s(1)+stp-1, stp:s(2)+stp-1);

Results.Labels = CellsImage;
Results.AnnotatedMask = TempImage;
Results.nCells = length(unique(CellsImage(:)))-1;
Results.Time = toc;

%% %%%%%%%%%%%%%%%%%%%%%%% Side Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [C] = CSC_DetectSeeds(A, opts, CellSize)

C = zeros(size(A));
DF = zeros(size(A));

[L,n] = bwlabel(A);
for h=1:n
    c = (L==h);
    crem = bwmorph(c, 'remove');
    [x,y] = find(crem>0);
    [xt,yt] = sortpixels(x,y,8);
    
    if(CellSize<=0)
        t = round(4*log(sum(crem(:))));
    else
        t = round(CellSize);
    end
    
    [N] = calcolacurvatura(xt,yt,c>0,t);
    
    N = N.*(N>0);
    
    for j=1:length(xt)
        DF(xt(j), yt(j)) = N(j);
    end
end


SeedLab = C;
[L,n] = bwlabel(DF>0);
N = 500;
for i=1:n
    
    [X,Y] = find(L==i);
    
    if(length(X)>10)
        try
            [a, par] = fitellip(X,Y);
            
            [X, Y] = drawellip(par, N);
            x=round(mean(X));
            y=round(mean(Y));
            
            C(x,y) = 1;
            SeedLab(x,y) = i;
            
        catch
            sprintf('X,Y vectors to small (%d points) to compute the ellipse', length(X))
        end
    end
    
end


function [Seeds, Clusters,nclusters] = CSC_ClusterizePoints(P, th)

n = max(size(P));

if(n > 1)
    M = zeros(n, n);
    for i=1:n
        for j=1:n
            M(i,j)=abs(P(i,1)-P(j,1))+abs(P(i,2)-P(j,2));
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
    
    if(dmin <= th)
        pos = find(Clusters(dminpos,:)==0,1);
        Clusters(dminpos, pos) = i;
        Assigned(i)=1;
    else
        nclusters = nclusters + 1;
        Clusters(nclusters, 1) = i;
        Assigned(i) = 1;
    end
end

Seeds = zeros(nclusters,2);

for i=1:nclusters
    pos = find(Clusters(i,:)==0,1)-1;
    Seeds(i,1) = round(mean(P(Clusters(i,1:pos),1)));
    Seeds(i,2) = round(mean(P(Clusters(i,1:pos),2)));
end


% drawellip
function [X, Y] = drawellip(a, N)

   %  convert to standard form: ((x-cx)/r1)^2 + ((y-cy)/r2)^2 = 1
   % rotated by theta
   v = solveellipse(a);

   % draw ellipse with N points   
   dx = 2*pi/N;
   theta = v(5);
   R = [ [ cos(theta) sin(theta)]', [-sin(theta) cos(theta)]'];
   for i = 1:N
        ang = i*dx;
        x = v(1)*cos(ang);
        y = v(2)*sin(ang);
        d1 = R*[x y]';
        X(i) = d1(1) + v(3);
        Y(i) = d1(2) + v(4);
   end
   
   
% fitellip 
function [a, par] = fitellip(X,Y)

   % normalize data
   mx = mean(X);
   my = mean(Y);
   sx = (max(X)-min(X))/2;
   sy = (max(Y)-min(Y))/2;
   x = (X-mx)/sx;
   y = (Y-my)/sy;
   
   % Build design matrix
   D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ];

   % Build scatter matrix
   S = D'*D;

   % Build 6x6 constraint matrix
   C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;
   
   % Solve eigensystem
   [gevec, geval] = eig(S,C);

   % Find the negative eigenvalue
   [NegR, NegC] = find(geval < 0 & ~isinf(geval));

   % Extract eigenvector corresponding to positive eigenvalue
   A = gevec(:,NegC);
   
   if(isempty(A))
       par = [];
       return;
   end

   % unnormalize
   par = [
        A(1)*sy*sy,   ...
        A(2)*sx*sy,   ...
        A(3)*sx*sx,   ...
        -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
        -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
        A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
                - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
                + A(6)*sx*sx*sy*sy   ...
       ]';


   % Convert to geometric radii, and centers


thetarad = 0.5*atan2(par(2),par(1) - par(3));
cost = cos(thetarad);
sint = sin(thetarad);
sin_squared = sint.*sint;
cos_squared = cost.*cost;
cos_sin = sint .* cost;

Ao = par(6);
Au =   par(4) .* cost + par(5) .* sint;
Av = - par(4) .* sint + par(5) .* cost;
Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;

% ROTATED = [Ao Au Av Auu Avv]

tuCentre = - Au./(2.*Auu);
tvCentre = - Av./(2.*Avv);
wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;

uCentre = tuCentre .* cost - tvCentre .* sint;
vCentre = tuCentre .* sint + tvCentre .* cost;

Ru = -wCentre./Auu;
Rv = -wCentre./Avv;

Ru = sqrt(abs(Ru)).*sign(Ru);
Rv = sqrt(abs(Rv)).*sign(Rv);

a = [uCentre, vCentre, Ru, Rv, thetarad];



% solveellipse
function v = solveellipse(a)

        % get ellipse orientation
        theta = atan2(a(2),a(1)-a(3))/2;

        % get scaled major/minor axes
        ct = cos(theta);
        st = sin(theta);
        ap = a(1)*ct*ct + a(2)*ct*st + a(3)*st*st;
        cp = a(1)*st*st - a(2)*ct*st + a(3)*ct*ct;

        % get translations
        T = [[a(1) a(2)/2]' [a(2)/2 a(3)]'];
        t = -inv(2*T)*[a(4) a(5)]';
        cx = t(1);
        cy = t(2);

        % get scale factor
        val = t'*T*t;
        scale = 1 / (val- a(6));

        % get major/minor axis radii
        r1 = 1/sqrt(scale*ap);
        r2 = 1/sqrt(scale*cp);
        v = [r1 r2 cx cy theta]';

        function [Xout, Yout] = sortpixels(X,Y,lftn)

Xout = 0;
Yout = 0;

L = zeros(max(X), max(Y));

for i=1:length(X)
    L(X(i), Y(i)) = 1;
end

% B = bwtraceboundary(L, [X(1) Y(1)], 'S');
P = bwboundaries(L, 8, 'noholes');

B = P{1};

Xout = B(:,1)';
Yout = B(:,2)';

% figure(1001);
% clf
% hold on;
% plot(Xout, Yout, 'b');

if(lftn>2)
    Xout = liftvector(Xout, lftn);
    Yout = liftvector(Yout, lftn);
end

% plot(Xout, Yout, 'g');

function W = liftvector(V, k)

n = length(V);
p = floor(k/2);

Vext = [V(n-p:n),V, V(1:p)];
W = V;

for i=1:n
    v = Vext(i:i+p);
    m = mean(v);
    W(i)=round(m);
end

        
function [curva] = calcolacurvatura(X,Y, Ima, stp)

n = length(X);

curva = 0;  

for k=1:n
    
    p=rem(k+stp, n)+1;

    xm1 = round((X(k)+X(p))/2);
    ym1 = round((Y(k)+Y(p))/2);
    
    p2 = rem(round(k+stp/2), n)+1; 
    xm2 = X(p2);
    ym2 = Y(p2);
    
    if(Ima(xm1, ym1)==0)
        segno = -1;
    else
        segno = +1;
    end
    
    curva(p2) = segno * sqrt((xm1-xm2)^2 + (ym1-ym2)^2);
    
end


if(isempty(max(curva(curva>0)))==0)
    curva = curva/max(curva(curva>0));
end
