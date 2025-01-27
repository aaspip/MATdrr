% Demonstration script for 
% State-of-the-art DAS denoising via the localized damped rank-reduction method (LDRR)
% This script takes about 10 minutes
% 
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
%  Chen et al., 2022, 3D seismic diffraction separation and imaging using the local rank-reduction method, IEEE Transactions on Geoscience and Remote Sensing, 60, 4507110.
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497–V506.
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  Chen et al., 2023, DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology, Computers & Geosciences, 180, 105440.
% 
% Dependency MATdrr
% svn co https://github.com/aaspip/MATdrr/trunk ./MATdrr
% or git clone https://github.com/aaspip/MATdrr ./
%
% Related pacakge:
% https://github.com/chenyk1990/dasmrrcoh
% https://github.com/chenyk1990/dasmrrcoh-dataonly

clc;clear;close all;
addpath(genpath('../matdrr'));

if ~isdir('fig')
    mkdir('fig');
end

if ~isdir('processed')
    mkdir('processed');
end

%% Download data first
%% The whole dataset in the folder "raw" can be downloaded from https://github.com/chenyk1990/dasmrrcoh-dataonly
%% A complete processing workflow can be found at https://github.com/chenyk1990/dasmrrcoh-dataonly or https://github.com/chenyk1990/dasmrrcoh
%
% https://github.com/chenyk1990/dasmrrcoh-dataonly/tree/main/raw/2017-06-23T22:03:22.380000Z_mag2.17.mat

name='2017-06-23T22:03:22.380000Z_mag2.17.mat'
ieq=3;
load(name);
eq=data;
d_bp=drr_bandpass(eq',1/250,0,20)';
d_bpmf=drr_mf(d_bp,5,1,1);
%% LDRR
n1win=1024;n2win=800;n3win=1;
n1win=512;n2win=200;n3win=1;
r1=0.5;r2=0.5;r3=0.5;
d_bpmfmrr=drr3d_win(d_bpmf',0,50,1/250,2,4,0,n1win,n2win,n3win,r1,r2,r3)';
save(sprintf('processed/eq%d.mat',ieq),'d_bp','d_bpmf','d_bpmfmrr');

[n1,n2]=size(data);
t=[0:n2-1]*(1/250);
x=1:n1;

figure('units','normalized','Position',[0.2 0.4 0.7, 0.6],'color','w');
ax1=subplot(2,2,1);
drr_imagesc(eq,95,1,t,x);colormap(ax1,seis);
name(end-3:end)=[];
title(name,'Interpreter', 'none','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'a)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

ax2=subplot(2,2,2);
drr_imagesc(d_bp,95,1,t,x);colormap(ax2,seis);
title('BP','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
% ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'b)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

ax3=subplot(2,2,3);
drr_imagesc(d_bpmf,95,1,t,x);colormap(ax3,seis);title('BP+MF');
title('BP+MF','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'c)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
xlabel('Time (s)','Fontsize',14,'fontweight','bold');

ax4=subplot(2,2,4);
drr_imagesc(d_bpmfmrr,95,1,t,x);colormap(ax4,seis);title('BP+MF+MRR');
title('BP+MF+MRR','Fontsize',14,'fontweight','bold');
% xlabel('Time (s)','Fontsize',14,'fontweight','bold');
% ylabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-5,-100,'d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
xlabel('Time (s)','Fontsize',14,'fontweight','bold');
print(gcf,'-depsc','-r300','test_matdrr_drr2d_win_dassafod.eps');
print(gcf,'-dpng','-r300','test_matdrr_drr2d_win_dassafod.png');
