function [ D1 ] = drr3drecon_win(D,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a,n1win,n2win,n3win,r1,r2,r3)
%  DRR3D_WIN: DRR3D in windows
%
%  IN   D:   	intput 3D data
%       MASK:   sampling mask (consistent with the POCS based approaches)
%       flow:   processing frequency range (lower)
%       fhigh:  processing frequency range (higher)
%       dt:     temporal sampling interval
%       N:      number of singular value to be preserved
%       K:      damping factor
%       Niter:  number of maximum iteration
%       eps:    tolerence (||S(n)-S(n-1)||_F<eps)
%       verb:   verbosity flag (default: 0)
%       mode:   mode=1: denoising and reconstruction
%               mode=0: reconstruction only
%       a:      weight vector
%       n1win:  window size
%       n2win:  window size
%       n3win:  window size
%       r1:  	overlapping ratio (default, 0.5)
%       r2:  	overlapping ratio (default, 0.5)
%       r3:  	overlapping ratio (default, 0.5)
%
%  OUT  D1:  	output data
%
%  Copyright (C) 2013 The University of Texas at Austin
%  Copyright (C) 2013 Yangkang Chen
%  Modified 2015 by Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%  REFERENCES
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  Chen et al., 2023, DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology, Computers & Geosciences, 180, 105440.

% DEMO: demos/test_matdrr_drr3d_win.m

if nargin==0
    error('Input data must be provided!');
end

[nt,nx,ny]=size(D);

if nargin==2
    flow=1;
    fhigh=124;
    dt=0.004;
    N=1;
    K=2;
    Niter=10;
    a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing;
    eps=0.00001;
    mode=1;
    verb=0;
    n1win=nt;
    n2win=nx;
    n3win=ny;
    r1=0.5;
    r2=0.5;
    r3=0.5;
end

param.dt=dt;
param.flow=flow;
param.fhigh=fhigh;
param.N=N;
param.K=K; %The N and K follows the definitions in Bai et al. (2020)
param.verb=verb;
param.mode=mode;
param.niter=Niter;
param.a=a;
param.eps=eps;

D1=drr_win3dmask(@localdrr3drecon, mask, param, D, n1win, n2win, n3win, r1, r2, r3);


return



