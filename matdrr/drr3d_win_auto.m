function [ D1 ] = drr3d_win_auto(D,flow,fhigh,dt,N,K,verb,n1win,n2win,n3win,r1,r2,r3,mode)
% DRR3D_WIN: DRR3D in windows
%
% IN   D:   	intput 3D data
%      flow:   processing frequency range (lower)
%      fhigh:  processing frequency range (higher)
%      dt:     temporal sampling interval
%      N:      number of singular value to be preserved
%      K:      damping factor
%      verb:   verbosity flag (default: 0)
%      n1win:  window size
%      n2win:  window size
%      n3win:  window size
%      r1:  	overlapping ratio (default, 0.5)
%      r2:  	overlapping ratio (default, 0.5)
%      r3:  	overlapping ratio (default, 0.5)
%
% OUT  D1:  	output data
%
% Copyright (C) 2013 The University of Texas at Austin
% Copyright (C) 2013 Yangkang Chen
% Modified 2015 by Yangkang Chen
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
% REFERENCES
% Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
% Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
% Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
% Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
% Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
% Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
% Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
% Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
% Chen et al., 2023, DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology, Computers & Geosciences, 180, 105440.
%
% DEMO: demos/test_matdrr_drr3d_win.m

if nargin==0
    error('Input data must be provided!');
end

[nt,nx,ny]=size(D);

if nargin==1
    flow=1;
    fhigh=124;
    dt=0.004;
    N=1;
    K=2;
    verb=0;
    n1win=nt;
    n2win=nx;
    n3win=ny;
    r1=0.5;
    r2=0.5;
    r3=0.5;
    mode=2;
end

param.dt=dt;
param.flow=flow;
param.fhigh=fhigh;
param.N=N;
param.K=K; %The N and K follows the definitions in Bai et al. (2020)
param.verb=verb;
param.mode=mode;

D1=drr_win3d(@localdrr3d_auto, param, D, n1win, n2win, n3win, r1, r2, r3);


return



