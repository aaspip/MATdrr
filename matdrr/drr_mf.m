function [ D1 ] = das_mf(D,nfw,ifb,axis)
%DRR_MF: median filter along first or second axis for 2D profile
%  IN   D:   	intput data 
%       nfw:    window size 
%       ifb:    if use padded boundary (if not, zero will be padded)
%       axis:    temporal sampling interval
%      
%  OUT   D1:  	output data
% 
% MIT License
% 
% Copyright (C) 2014 Yangkang Chen
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
% References
% Huang et al., 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.
% Chen et al., 2020, Deblending of simultaneous-source data using a structure-oriented space-varying median filter, Geophysical Journal International, 222, 1805–1823.
% Gan et al., 2016, Separation of simultaneous sources using a structural-oriented median filter in the flattened dimension, Computers & Geosciences, 86, 46-54.
% Chen, Y., 2015, Deblending using a space-varying median filter, Exploration Geophysics, 46, 332-341.

if nargin==0
 error('Input data must be provided!');
end

if nargin==1
 nfw=7;
 ifb=1;
 axis=2;
end;

if nargin==2
 ifb=1;
 axis=2;    
end

% nfw should be odd
if mod(nfw,2)==0
    nfw=nfw+1;
end

if axis==2
   D=D.'; 
end

[n1,n2]=size(D);
nfw2=(nfw-1)/2;

if ifb==1
    D=[flipud(D(1:nfw2,:));D;flipud(D(n1-nfw2+1:n1,:))];
else
    D=[zeros(nfw2,n2);D;zeros(nfw2,n2)];    
end

% output data
D1=zeros(n1,n2);
for i2=1:n2
   for i1=1:n1
      D1(i1,i2)=median(D(i1:i1+nfw-1,i2)); 
   end 
end
if axis==2
    D1=D1.';
end
return
