function [] = drr_imagesc(data,pclip,mode,x,z)
%drr_imagesc: fast plot data using pclip
% 
% by Yangkang Chen
% Jan, 2016
% 
% Input:
%  data: input data
%  pclip: clip value (percential or exact)
%  mode=1: pclip; mode=2: clip
% 
% MIT License
% 
% Copyright (C) 2016 Yangkang Chen
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

if nargin==1
    pclip=99;
    mode=1;%using pclip;
end

if nargin==2
    mode=1;
end

% mi=min(min(abs(data)));
% ma=max(max(abs(data)));

if mode==1
    t=prctile(abs(data(:)),pclip);
    % figure;
    if nargin==5
        imagesc(x,z,data);caxis([-t,t]);colormap(cseis);
    else
        imagesc(data);caxis([-t,t]);colormap(cseis);
    end
    
else
    if nargin==5
        imagesc(x,z,data);caxis([-pclip,pclip]);colormap(cseis);
    else
        imagesc(data);caxis([-pclip,pclip]);colormap(cseis);
    end
    
end

end

function [map]=cseis()

map = [[0.5*ones(1,40),linspace(0.5,1,88),linspace(1,0,88),zeros(1,40)]',[0.25*ones(1,40),linspace(0.25,1,88),linspace(1,0,88),zeros(1,40)]',[zeros(1,40),linspace(0.,1,88),linspace(1,0,88),zeros(1,40)]'];

end



