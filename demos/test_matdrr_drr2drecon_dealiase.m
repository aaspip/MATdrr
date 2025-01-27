% Demonstration script for 
% 3D dealiased reconstruction via the damped rank-reduction method
%
%  Copyright (C) 2023 The University of Texas at Austin
%  Copyright (C) 2023 Yangkang Chen
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
%  References:   
%
%  [0] Huang, W., D. Feng, and Y. Chen, 2020, De‐aliased and de‐noise Cadzow filtering for seismic data reconstruction, Geophysical Prospecting, 68, 443-571.
%  [1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  [2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  [4] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  [5] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  [6] Chen et al., 2023, DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology, Computers & Geosciences, 180, 105440.

clc;clear;close all;
addpath(genpath('../matdrr'));


%% generate 3D synthetic data
a1=zeros(300,40);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
% b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
  t1(i)=round(140);
  t3(i)=round(-3*i+180);
  t4(i)=round(3*i+10);
  a1(t1(i):t1(i)+k-1,i)=b1; 
  a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b1;
end

d=a1+a3+a4;

%% add band-limitted noise and decimate data
randn('state',201315);
[nt,nx,ny]=size(d);
noise=randn(nt,nx,ny);
noise=drr_bandpass(noise,0.004,0,60);
dn=0.1*noise+d;

%% Doing the aliased reconstruction (2D densification)
d1=dn(:,1:2:end);
[n1,n2,n3]=size(d1);
d3=drr3drecon_dealiase(d1,0,50,0.004,3,4,20,2,2,1);%takes about ? minutes
d3=d3(:,:,1);


%% professional plot
dt=0.004;dx=1;
figure('units','normalized','Position',[0.2 0.4 0.5, 0.8]);
subplot(2,2,1)
t=[0:n1-1]*dt;
x=[1:n2]*dx;
drr_imagesc(d1,1,2,x,t);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
title('Under-sampled','Fontsize',16,'fontweight','bold');
text(-2.5,-0.1,'a)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','center');

subplot(2,2,2)
t=[0:n1-1]*dt;
x=[1:n2]*dx;
drr_imagesc(d3,1,2,[1:n2*2]*dx,t);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
title('Densified','Fontsize',16,'fontweight','bold');
text(-5,-0.1,'b)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','center');

subplot(2,1,2);
dcomp=[d,dn,d3,dn-d3];
drr_imagesc(dcomp,1,2,[1:n2*4]*dx,t);
% imagesc(1:n2,[0:n1-1]*dt,dcomp);colormap(seis);
ylabel('Time (s)','Fontsize',16,'fontweight','bold');
xlabel('Trace','Fontsize',16,'fontweight','bold');
title(sprintf('2D densification performance (Clean|Noisy|Recon|Error; SNR=%g dB)',yc_snr(d,d3)),'Fontsize',16,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
text(-5,-0.1,'c)','color','k','Fontsize',24,'fontweight','bold','HorizontalAlignment','center');

print(gcf,'-dpng','-r300','test_matdrr_drr2drecon_dealiase.png');
print(gcf,'-depsc','-r200','test_matdrr_drr2drecon_dealiase.eps');
