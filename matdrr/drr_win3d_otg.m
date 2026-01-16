function [ dout ] =drr_win3d_otg(oper, param, din, n1win, r1)
% DRR_WIN3D_OTG: Processing in 3D windows for OTG case (only time windows are allowed)
%
% din:          input data (2D)
% oper:         operator
% param:        parameters of operator
% n1win:        first window length
% r1:           first overlapping ratio
%
% dout:     output data
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
% see also
% win2d.m,win3d_mask.m
%
% Example:
% test/test_win2d_fxmssa.m
% test/test_win3d_fxymssa.m
%
% REFERENCES
% Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
% Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
% Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
% Chen et al., 2023, DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology, Computers & Geosciences, 180, 105440.



[n1,n2]=size(din);

if nargin==3
    n1win=n1;
    r1=0.5;
end

if nargin==4
    r1=0.5;
end

nov1=(1-r1)*n1win;  % non-overlapping size 1
ov1=r1*n1win;       % overlapping size 1


n1pad=n1win;        %padding size 1
nw1=1;
while n1pad<n1
    n1pad=n1pad+nov1;nw1=nw1+1;
end

D1=zeros(n1pad, n2);D1(1:n1,1:n2)=din;
D2=zeros(n1pad,param.nx,param.ny);

for iw1=0:nw1-1
    s1=iw1*nov1;%s2=iw2*nov2;s3=iw3*nov3;
    dtmp=D1(s1+1:s1+n1win,:);

    %uncomment this line for checking the correctness (checked 100% correct)
    dtmp = feval(oper,dtmp,param);
    %only valid for space-independent param
    %for reconstruction, with mask, param needs to be changed

    dtmp=win_weight3d(dtmp,iw1,0,0,nw1,1,1,n1win,param.nx,param.ny,ov1,0,0);

    D2(s1+1:s1+n1win,:,:)=D2(s1+1:s1+n1win,:,:)+dtmp;
end

dout=D2(1:n1,:,:);

return


function [dout ]= win_weight3d(din,iw1,iw2,iw3,nw1,nw2,nw3,n1win,n2win,n3win,ov1,ov2,ov3)
%      follow exactly usr/cm/win.c
%      float din /*input data*/,
% 		int iw1 /*starting window 1 in dst*/,
% 		int iw2 /*starting window 2 in dst*/,
% 		int iw3 /*starting window 3 in dst*/,
% 		int nw1 /*no of windows 1 in src*/,
% 		int nw2 /*no of windows 2 in src*/,
% 		int nw3 /*no of windows 3 in src*/,
% 		int n1win /*window length 1 in src*/,
% 		int n2win /*window legnth 2 in src*/,
% 		int n3win /*window legnth 3 in src*/,
% 		int ov1 /*copy length in axis1*/,
% 		int ov2 /*copy length in axis2*/,
% 		int ov3 /*copy length in axis3*/)


if iw3~=0
    for i1=0:n1win-1
        for i2=0:n2win-1
            for i3=0:ov3-1
                din(i1+1,i2+1,i3+1)=din(i1+1,i2+1,i3+1)*(i3+1)/(ov3+1);
            end
        end
    end

end

if iw3~=nw3-1
    for i1=0:n1win-1
        for i2=0:n2win-1
            for i3=0:ov3-1
                din(i1+1,i2+1,n3win-ov3+i3+1) = din(i1+1,i2+1,n3win-ov3+i3+1)*(ov3-i3)/(ov3+1);
            end
        end
    end
end



if iw2~=0
    for i3=0:n3win-1
        for i1=0:n1win-1
            for i2=0:ov2-1
                din(i1+1,i2+1,i3+1)=din(i1+1,i2+1,i3+1)*(i2+1)/(ov2+1);
            end
        end
    end

end

if iw2~=nw2-1
    for i3=0:n3win-1
        for i1=0:n1win-1
            for i2=0:ov2-1
                din(i1+1,n2win-ov2+i2+1,i3+1) = din(i1+1,n2win-ov2+i2+1,i3+1)*(ov2-i2)/(ov2+1);
            end
        end
    end
end

if iw1~=0
    for i3=0:n3win-1
        for i2=0:n2win-1
            for i1=0:ov1-1
                din(i1+1,i2+1,i3+1)=din(i1+1,i2+1,i3+1)*(i1+1)/(ov1+1);
            end
        end
    end

end

if iw1~=nw1-1
    for i3=0:n3win-1
        for i2=0:n2win-1
            for i1=0:ov1-1
                din(n1win-ov1+i1+1,i2+1,i3+1)=din(n1win-ov1+i1+1,i2+1,i3+1)*(ov1-i1)/(ov1+1);
            end
        end
    end
end

dout=din;
return












