function [D,f,k]=drr_fk2d(d,dt,dx)
% drr_fk2d: FK spectrum (2D) for spectrum comparison
% 
% INPUT  
% d:   	intput data
% dt:     time sampling
% dx:     space sampling
%          
% OUTPUT  
% D:  	output spectrum
%
% DEMO
% demos/test_matdrr_drr2d_win.m
%
% 
% MIT License
% 
% Copyright (C) 2020 Yangkang Chen
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
   dt=0.004;
   dx=1;
end

[n1,n2]=size(d);
nf = 2^nextpow2(n1);nf2=nf/2;
nk = 2^nextpow2(n2);

D=fftshift(fft(fft(d,nf,1),nk,2));
D=D(nf/2+1:end,:);

f_nq=1/dt/2;
k_nq=1/dx/2;
f=linspace(0,f_nq,nf2);
k=linspace(-k_nq,k_nq,nk);


return


