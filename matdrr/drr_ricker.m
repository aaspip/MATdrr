function [w,tw] = drr_ricker(f,dt,tlength)
% drr_ricker: Ricker wavelet of central frequency f0.
%
% INPUT:
% f : central freq. in Hz (f <<1/(2dt) )
% dt: sampling interval in sec
% tlength : the duration of wavelet in sec
%
% OUTPUT: 
% w:  the Ricker wavelet
% tw: time axis
%
% Example
%
%  [w,tw] = drr_ricker(10,0.004,0.2);
%   plot(tw,w);
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

if nargin==3
    nw=floor(tlength/dt)+1;
else
    nw=2.2/f/dt;
    nw=2*floor(nw/2)+1;
end
nc=floor(nw/2);
w = zeros(nw,1);

k=[1:1:nw]';

alpha = (nc-k+1).*f*dt*pi;
beta=alpha.^2;
w = (1.-beta.*2).*exp(-beta);

if nargout>1;
    tw = -(nc+1-[1:1:nw])*dt;
end

