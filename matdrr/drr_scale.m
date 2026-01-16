function [ D1 ] = drr_scale(D,N,dscale)
%drr_scale: Scale the data up to the Nth dimension = sfscale axis=N
% IN   D:   	intput data
%      N:      number of dimension for scaling
%              default: N=2
% 		dscale:  Scale by this factor
%      (does not include the rscale and pclip functions (not convenient actually))
%
% OUT   D1:  	output data
%
% Copyright (C) 2015 The University of Texas at Austin
% Copyright (C) 2015 Yangkang Chen
% Modified by Yangkang Chen on Jan, 2020
% 
% MIT License
% 
% Copyright (C) 2015 Yangkang Chen
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

if nargin==0
    error('Input data must be provided!');
end

if nargin==1
    N=2;
    dscale=1.0;
end

if nargin==2
    dscale=1.0;
end


[n1,n2,n3]=size(D);
D1=D;
switch N
    case 1
        for i3=1:n3
            for i2=1:n2
                D1(:,i2,i3)=D1(:,i2,i3)/max(abs(D1(:,i2,i3)));
            end
        end
    case 2
        for i3=1:n3
            D1(:,:,i3)=D1(:,:,i3)/max(max(abs(D1(:,:,i3))));
        end
    case 3
        D1=D1/max(max(max(abs(D1))));
        
    case 0
        
        D1=D1*dscale;
        
    otherwise
        error('Invalid argument value N.')
end


return
