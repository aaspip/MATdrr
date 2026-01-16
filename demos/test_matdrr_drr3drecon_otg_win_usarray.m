% Demo script for OTG teleseismic denoising and reconstruction
% as introduced in Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
% Another version of this example is at 
% https://github.com/chenyk1990/reproducible_research/tree/master/nonMada/usarray
% 
% This takes about 10 minutes
% 
% Written by Yangkang Chen
% Feb, 2018
% Modified on Dec, 2020
% Further polished on July, 2022, Feb, 2023
%
% MIT License
% 
% Copyright (C) 2018 Yangkang Chen
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
% Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
% Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497?V506.
% Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
% Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
% Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
% Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
% Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
% Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
% Chen et al., 2023, DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology, Computers & Geosciences, 180, 105440.

clc;clear;close all;
addpath(genpath('../matdrr'));

%% Please download data from https://github.com/aaspip/data/blob/main/usarray_200901181411_wfm.mat
load('usarray_200901181411_wfm.mat');
% d, dists(shot receiver distance/offset in degree), stla, stlo, t 

% figure;drr_imagesc([d(:,:)]);

%% rm bad trace
inds=[18,41,70];
d(:,inds)=[];d=drr_scale(d);
stlo(inds)=[];
stla(inds)=[];
d0=d(:,105:433);
stlo0=stlo(105:433);
stla0=stla(105:433);



d0=d0(1001:3000,:);


%% 3D processing/reconstruction
mla=[33,49];
mlo=[-116,-102];
%binning
[d3d,x1,y1,mask]=drr_bin3d(d0,stlo0,stla0,16,28,mlo(1),mla(1),mlo(2),mla(2));
[stlo1,stla1]=meshgrid(x1,y1);

%% global processing
flow=0;fhigh=0.5;dt=1;N=8;Niter=10;mode=1;verb=1;eps=0.00001;K=4;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=drr3drecon(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a);


% figure;drr_imagesc([d3d(:,:),d1(:,:)]);
d2=drr3drecon_otg(d0,stlo0,stla0,16,28,mlo(1),mla(1),mlo(2),mla(2),flow,fhigh,dt,N,K,Niter,eps,verb,mode);


%% fixed-rank windowed reconstruction
[n1,n2,n3]=size(d3d);
n1win=200;r1=0.5;N=4;K=4;
% d11=drr3drecon_win(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a,n1win,n2,n3,r1,0,0);
d22=drr3drecon_otg_win(d0,stlo0,stla0,16,28,mlo(1),mla(1),mlo(2),mla(2),flow,fhigh,dt,N,K,Niter,eps,verb,mode,n1win,r1);


% figure;drr_imagesc([d3d(:,:),d1(:,:),d2(:,:)]);


%% plot the grid
figure;
plot(stlo0,stla0,'o','linewidth',2);hold on;
plot(stlo1(:),stla1(:),'*','linewidth',2);
ylim([32,50]);xlim([-117,-85]);
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Longitude (^o)','Fontsize',20);
title('Irregular V.S. regular grids','Fontsize',20);
plot([stlo1(1,5),stlo1(1,5)],[32,50],'linewidth',2,'color','r');
plot([stlo1(1,12),stlo1(1,12)],[32,50],'linewidth',2,'color','g');
legend('Irregular','Regular','Slice 1','Slice 2','Fontsize',20,'location','best');
print(gcf,'-dpng','-r300','test_matdrr_drr3drecon_otg_usarray_grid.png');
print(gcf,'-depsc','-r300','test_matdrr_drr3drecon_otg_usarray_grid.eps');



%% plot the results
x1=750;
x2=1200;
y1=45;
y2=49.5;

ilon=12;
dtest=squeeze(d3d(:,ilon,:));
figure('units','normalized','Position',[0.2 0.4 1.0, 1.2],'color','w');
subplot(2,4,1);
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Raw','Fontsize',20);
drr_framebox(x1,x2,y1,y2,'r',2);
text(0,52,'a)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,2);
dtest=squeeze(d1(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('DRR','Fontsize',20);
drr_framebox(x1,x2,y1,y2,'r',2);
text(0,52,'b)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,3);
dtest=squeeze(d2(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('OTG','Fontsize',20);
drr_framebox(x1,x2,y1,y2,'r',2);
text(0,52,'c)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,4);
dtest=squeeze(d22(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([31.5,50.2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Windowed OTG','Fontsize',20);
drr_framebox(x1,x2,y1,y2,'r',2);
text(0,52,'d)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

%% zoomed comparison
subplot(2,4,5);
dtest=squeeze(d3d(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([y1,y2]);xlim([x1,x2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Raw','Fontsize',20);
text(660,50,'e)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

% print(gcf,'-depsc','-r400','us_lon1_z.eps');  

subplot(2,4,6);
dtest=squeeze(d1(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([y1,y2]);xlim([x1,x2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('DRR','Fontsize',20);
text(660,50,'f)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,7);
dtest=squeeze(d2(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([y1,y2]);xlim([x1,x2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('OTG','Fontsize',20);
text(660,50,'g)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,8);
dtest=squeeze(d22(:,ilon,:));
drr_wigbh(dtest,stla1,t(1001:3000),20);
ylim([y1,y2]);xlim([x1,x2]);
ylabel('Latitude (^o)','Fontsize',20);
xlabel('Time (s)','Fontsize',20);
title('Windowed OTG','Fontsize',20);
text(660,50,'g)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');


%% following are all annotations
% Create textbox
print(gcf,'-dpng','-r300','test_matdrr_drr3drecon_otg_win_usarray.png');
print(gcf,'-depsc','-r300','test_matdrr_drr3drecon_otg_win_usarray.eps');


%% plot 3D 
% figure('units','normalized','Position',[0.2 0.4 0.5, 0.8],'color','w');
% z=1:2000;x=1:16;y=1:28;
% subplot(3,2,1);drr_plot3d(d3d,[1000,8,14],z,x,y);caxis([-0.1,0.1]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Data','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'a)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
% subplot(3,2,3);drr_plot3d(d1,[1000,8,14],z,x,y);caxis([-0.1,0.1]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Ground-truth diffraction','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'b)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
% % subplot(3,2,4);drr_plot3d(data-diffr,[100,20,20],z,x,y);caxis([-0.1,0.1]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('Ground-truth reflection','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'c)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
% subplot(3,2,5);drr_plot3d(d2,[1000,8,14],z,x,y);caxis([-0.1,0.1]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('LDRR diffraction','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'d)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');
% % subplot(3,2,6);drr_plot3d(data-diffr,[100,20,20],z,x,y);caxis([-0.1,0.1]);xlabel('X (km)','Fontsize',15);ylabel('Y (km)','Fontsize',15);zlabel('Time (s)','Fontsize',15);title('LDRR reflection','Fontsize',15,'fontweight','normal');set(gca,'Linewidth',2,'Fontsize',15);text(-0.5,-0.2, -0.3,'e)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

% 
% figure('units','normalized','Position',[0.2 0.4 0.5, 0.8],'color','w');
% subplot(1,3,1);drr_imagesc(squeeze(d3d(1000,:,:)));
% 
% subplot(1,3,2);drr_imagesc(squeeze(d1(1000,:,:)));
% subplot(1,3,3);drr_imagesc(squeeze(d2(1000,:,:)));



