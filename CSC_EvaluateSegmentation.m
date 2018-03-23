% [Measures] = CSC_EvaluateSegmentation(A, GT, indices)
% is a function of the CSC package that evaluate the performance of a
% segmentation algorithm by comparing the input image A with a ground truth
% GT.
% For mor details about this algorithm, please refer
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
%   A           A binary image or a label image of size HxW.
%   
%   GT          A binary image or a label image of size HxW.
%
%   indices     A vector of integers (0/1). If indices(i) is 1 the
%               corresponding measure is computed. The sorted list of
%               measures is the following:
%
%               1) Fmeasure
%               2) Random Index
%               3) Jaccard Distance
%               4) Jaccard Similarity Coefficient
%               5) Housdorff MEtric
%               6) Normalized S Distance
%               7) Number of Cells Split
%               8) Number of Cells Merged
%               9) Number of Cells Added
%               10) Number of Cells Missing
%               11) Recognition Error Rate
%               12) Global Consistency of Edges
%               13) Variation of Information
%               14) Boundary Displacement Error
%
%
%   Output
%   ------------------ 
%   Measures    A structure with 18 fields, where the ith measure
%               has value Nan if indices(i)=0 or an integer/real number if 
%               indices(i)=1. Notice that Accuracy, Precision, Sensitivity 
%               and Specificity are always computed and they do not have a 
%               corresponding position in the vector indices.
%
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

function [Measures] = CSC_EvaluateSegmentation(A, GT, indices)

Measures = struct;

A = double(A);
GT = double(GT);

% statistiche sulla prima ground trouth
TP1 = (A>0).*(GT>0);
TN1 = (A==0).*(GT==0);
FP1 = (A>0).*(GT==0);
FN1 = (A==0).*(GT>0);


TP1 = sum(sum(TP1)); 
TN1 = sum(sum(TN1));
FP1 = sum(sum(FP1)); 
FN1 = sum(sum(FN1)); 

% Measures.TP1 = TP1;
% Measures.TN1 = TN1;
% Measures.FP1 = FP1;
% Measures.FN1 = FN1;

Measures.Accuracy = (TP1+TN1)/(TP1+FN1+TN1+FP1);
Measures.Precision = TP1/(TP1+FP1);
Measures.Sensitivity = TP1/(TP1+FN1);
Measures.Specificity = TN1/(FP1+TN1);

% 1) Fmeasure
if(indices(1))
    Measures.Fmeasure = 2*Measures.Precision*Measures.Sensitivity / (Measures.Precision+Measures.Sensitivity);
else
    Measures.Fmeasure = NaN;
end

% 2) Rand Index
if(indices(2))
    Measures.RI = RandIndex(A, GT);
else
    Measures.RI = NaN;
end

% 3) Jaccard Distance
if(indices(3))
%     q = 1-jaccard(A,GT);
%     q = q(isnan(q)==0);
%     Measures.JD = mean(q);
    Measures.JD = 1-sum(sum(A.*GT))/sum(sum((A+GT)>0));
else
    Measures.JD = NaN;
end

% 4) Jaccard Similarity Coefficient
if(indices(4))
    Measures.JI = JaccardIndex(A,GT);
else
    Measures.JI = NaN;
end

[HM, NSD, Split, Merged, Added, Missing]= ComputeMetrics(A,GT);

% 5) Housdorff MEtric
if(indices(5))
    Measures.HM = HM;
else
    Measures.HM = NaN;
end

% 6) Normalized S Distance
if(indices(6))
    Measures.NSD = NSD;
else
    Measures.NSD = NaN;
end

% 7) Number of Cells Split
if(indices(7))
    Measures.Split = Split;
else
    Measures.Split = NaN;
end

% 8) Number of Cells Merged
if(indices(8))
    Measures.Merged = Merged;
else
    Measures.Merged = NaN;
end

% 9) Number of Cells Added
if(indices(9))
    Measures.Added = Added;
else
    Measures.Added = NaN;
end

% 10) Number of Cells Missing
if(indices(10))
    Measures.Missing = Missing;
else
    Measures.Missing = NaN;
end

% 11) Recognition Error Rate
if(indices(11))
    Measures.RER = RecognitionErrorRate(A, GT);
else
    Measures.RER = NaN;
end

% 12) Global Consistency of Edges
if(indices(12))
    Measures.GCE = GlobalConsistencyEdges(A,GT);
else
    Measures.GCE = NaN;
end

% 13) Variation of Information
if(indices(13))
    Measures.VOI = VariationOfInformation(A, GT);
else
    Measures.VOI = NaN;
end

% 14) Boundary Displacement Error
if(indices(14))
    Measures.BDE = BoundaryDisplacementError(A,GT);
else
    Measures.BDE = NaN;
end

%% SIDE FUNCTIONS
function [JI, RI] = JaccardIndex(R,S)

a=0;
b=0;
c=0;
RLab = unique(R(:));
SLab = unique(S(:));
n = length(R(:));

m = max(length(RLab), length(SLab));
for i=1:m
    if(i<=length(RLab))
        rl = find(R(:)==RLab(i));
    else
        rl=[];
    end
    
    if(i<=length(SLab))
        sl = find(S(:)==SLab(i)); 
    else
        sl = [];
    end
    
    Rr = R(sl);
    rok = (length(sl)>0);
    Sr = S(rl);
    sok = (length(rl)>0);
    rlen = length(sl);
    slen = length(rl);

    t = max(rlen, slen);

    for j=1:t
        for h=j+1:t
            if(sok && j<=slen && h<=slen)
                if(Sr(j)==Sr(h))
                    a=a+1;
                else
                    c=c+1;
                end
            end
            
             if(rok && j<=rlen && h<=rlen)
                if(Rr(j)~=Rr(h))
                    b=b+1;
                end
            end
            
        end
    end
end

d = n*(n-1)/2 - (a+b+c);
JI = (a+d)/(b+c+d);
RI = (a+d)/(a+b+c+d);



function [RI] = RandIndex(A, GT)
[W,H]=size(A);

if min(min(A)) < 1
    A = A - min(min(A)) + 1;
end
if min(min(GT)) < 1
    GT = GT - min(min(GT)) + 1;
end

sc1=max(max(A));
sc2=max(max(GT));

n=zeros(sc1,sc2);

for i=1:W
    for j=1:H
        u=A(i,j);
        v=GT(i,j);
        n(u,v)=n(u,v)+1;
    end
end

N = sum(sum(n));
n_u=sum(n,2);
n_v=sum(n,1);
N_choose_2=N*(N-1)/2;
RI = 1 - ( sum(n_u .* n_u)/2 + sum(n_v .* n_v)/2 - sum(sum(n.*n)) )/N_choose_2;


function [GCE] = GlobalConsistencyEdges(A, GT)
[W,H]=size(A);

if min(min(A)) < 1
    A = A - min(min(A)) + 1;
end
if min(min(GT)) < 1
    GT = GT - min(min(GT)) + 1;
end

sc1=max(max(A));
sc2=max(max(GT));

n=zeros(sc1,sc2);

for i=1:W
    for j=1:H
        u=A(i,j);
        v=GT(i,j);
        n(u,v)=n(u,v)+1;
    end
end

N = sum(sum(n));
m1 = sum(n,2);
m2 = sum(n,1);

E1 = 1 - sum( sum(n.*n,2) ./ (m1 + (m1 == 0)) ) / N;
E2 = 1 - sum( sum(n.*n,1) ./ (m2 + (m2 == 0)) ) / N;
GCE = min( E1, E2 );


function [VOI] = VariationOfInformation(A, GT)
[W,H]=size(A);

if min(min(A)) < 1
    A = A - min(min(A)) + 1;
end
if min(min(GT)) < 1
    GT = GT - min(min(GT)) + 1;
end

sc1=max(max(A));
sc2=max(max(GT));

n=zeros(sc1,sc2);

for i=1:W
    for j=1:H
        u=A(i,j);
        v=GT(i,j);
        n(u,v)=n(u,v)+1;
    end
end

N = sum(sum(n));
jnt = n / N; 
m2 = sum(jnt,1); 
m1 = sum(jnt,2); 
H1 = - sum( m1 .* log2(m1 + (m1 == 0) ) ); 
H2 = - sum( m2 .* log2(m2 + (m2 == 0) ) ); 
MI = sum(sum( jnt .* Log2QT( jnt, m1*m2 )  )); 
VOI = H1 + H2 - 2 * MI; 


function [log2qt] = Log2QT( A, GT )
log2qt = log2( (A + ((A==0).*GT) + (GT==0)) ./ (GT + (GT==0)) );



function [HM, NSD, Split, Merged, Added, Missing]= ComputeMetrics(R,S)

Rlab = unique([0 R(:)']);
Slab = unique([0 S(:)']);
StRlab = zeros(1, length(Slab));

n = length(Rlab);
m = length(Slab);
MRS = zeros(n,m);
MSR = zeros(m,n);

for k=2:n
    
    g = (R==Rlab(k));
    if(sum(g(:))>0)
        
        c = (S.*g);
        
        h = histc(c(g>0), Slab);
        mlab = min(find(h==max(h)));
        
        MRS(k,mlab)=1;
    end
end

StRlab(1)=1;
for k=2:m
    
    g = (S==Slab(k));
    c = R.*g;
    
    h = histc(c(g>0), Rlab);
    nlab = min(find(h==max(h)));
        
    MSR(k,nlab)=1;

    StRlab(k)=nlab;   

end

smrs = sum(MRS);
smsr = sum(MSR);
Split = sum(1.*(smsr>1));
Merged = sum(1.*(smrs>1));
Added = sum(MSR(:,1));
Missing = sum(MRS(:,1));

SplitLab = find(smsr>1);
MergedLab = find(smrs>1);
AddedLab = find(MSR(:,1)>0)';
MissingLab = find(MRS(:,1)>0)';

Svalidlab = ones(1, length(Slab));
Svalidlab(1)=0;

Rvalidlab = ones(1, length(Rlab));
Rvalidlab(1)=0;

cnt = 0;
HM = zeros(1, length(Slab));
NSD = zeros(1, length(Slab));
for k=2:length(Slab)
    
    num = 0;
    den = 0;
    HMtmp = 0;
    

    if(Svalidlab(k)==0 || Rvalidlab(StRlab(k))==0)
        continue;
    end
    
    Sb = (S==Slab(k));
    Rb = (R==Rlab(StRlab(k)));
    
    RB =  bwdist(1-Rb) + bwdist(Rb);

    [X,Y] = find((Rb+Sb)>0);
    
    cnt = cnt + 1;
    
    for i=1:length(X)
        if(Rb(X(i),Y(i))~=Sb(X(i),Y(i)))
            num = num+RB(X(i),Y(i));
            HM(cnt) = max(HM(cnt), RB(X(i),Y(i)));
        end
        den = den + RB(X(i),Y(i));
    end
    
    NSD(cnt) = num/den;

end

HM = sum(HM(1:cnt))/cnt;
NSD = sum(NSD(1:cnt))/cnt;

function [BDE] = BoundaryDisplacementError(A, GT)

if isempty(find(A~=A(1)))
    b1 = zeros(size(A));
    b1(1,:) = 1;
    b1(:,1) = 1;
    b1(end,:) = 1;
    b1(:,end) = 1;
else
    [cx,cy] = gradient(A);
    [bPx{1},bPy{1}] = find((abs(cx)+abs(cy))~=0);
    
    b1 = abs(cx) + abs(cy) > 0;
end

if isempty(find(GT~=GT(1)))
    b2 = zeros(size(GT));
    b2(1,:) = 1;
    b2(:,1) = 1;
    b2(end,:) = 1;
    b2(:,end) = 1;    
else    
    [cx,cy] = gradient(GT);
    [bPx{2},bPy{2}] = find((abs(cx)+abs(cy))~=0);
    
    b2 = abs(cx) + abs(cy) > 0;
end

D1 = bwdist(b1);
D2 = bwdist(b2);

d12 = sum(sum(b1 .* D2 ));
d21 = sum(sum(b2 .* D1 ));

avge12 = d12 / sum(sum(b1));
avge21 = d21 / sum(sum(b2));

BDE = (avge12 + avge21)/2;

function [rer] = RecognitionErrorRate(A, GT)

Tcell = (A-GT)>0;
Tback = (A-GT)<0;

n = prod(size(GT));

rer = sum(sum((Tcell+Tback)))/n * 100;

