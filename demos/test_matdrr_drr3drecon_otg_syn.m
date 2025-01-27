%  Demo script for OTG teleseismic denoising and reconstruction (synthetic example)
% 
%  Written by Yangkang Chen
%  April, 2024
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
% REFERENCES
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497?V506.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  Chen et al., 2023, DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology, Computers & Geosciences, 180, 105440.

clc;clear;close all;
addpath(genpath('../matdrr'));

%% rm bad trace
%% generate 3D synthetic data
a1=zeros(300,20);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
  t1(i)=round(140);
  t3(i)=round(-6*i+180);
  t4(i)=round(6*i+10);
  a1(t1(i):t1(i)+k-1,i)=b1; 
  a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b1;
end

temp=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
for j=1:20
    a4=zeros(300,20);
    for i=1:m
  t4(i)=round(6*i+10+3*j); 
  a4(t4(i):t4(i)+k-1,i)=b1;
  
  t1(i)=round(140-2*j);
  a1(t1(i):t1(i)+k-1,i)=b1;
    end
    shot(:,:,j)=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
end
plane3d=shot;
d=plane3d/max(max(max(plane3d)));

%% without noise
dn=d;

%% decimate
[nt,nx,ny]=size(d);
% ratio=0.5;
% mask=drr_genmask(reshape(d,nt,nx*ny),ratio,'c',201415);
% mask=reshape(mask,nt,nx,ny);
% d0=dn.*mask;
% 

% drr_plot3d(d0);

% stlo0=stlo(105:433);
% stla0=stla(105:433);
% d0=d0(1001:3000,:);

xx=1:nx;
yy=1:ny;
[stlo1,stla1]=meshgrid(xx,yy);

randn('state',2021);
tmpx=randn(nx*ny,1);
randn('state',2022);
tmpy=randn(nx*ny,1);
for ix=1:nx
    for iy=1:ny
        ii=(ix-1)*ny+iy;
        x(ii)=xx(ix)+tmpx(ii);
        y(ii)=yy(iy)+tmpy(ii);
    end
end

%% plot the grid
figure;
plot(x,y,'o','linewidth',2);hold on;
plot(stlo1(:),stla1(:),'*','linewidth',2);
ylim([-1,21]);xlim([-1,21]);
set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
ylabel('Y (trace)','Fontsize',16);
xlabel('X (trace)','Fontsize',16);
title('Irregular V.S. regular grids','Fontsize',16);
plot([stlo1(1,5),stlo1(1,5)],[-1,21],'linewidth',2,'color','r');
plot([stlo1(1,12),stlo1(1,12)],[-1,21],'linewidth',2,'color','g');
legend('Irregular','Regular','Slice 1','Slice 2','Fontsize',16,'location','best');
print(gcf,'-dpng','-r300','test_matdrr_drr3drecon_otg_syn_grid.png');
print(gcf,'-depsc','-r300','test_matdrr_drr3drecon_otg_syn_grid.eps');


%% generate the observed data (2D array)
par.x=x;
par.y=y;
par.nx=nx;
par.ny=ny;
par.ox=1;
par.oy=1;
par.mx=nx;
par.my=ny;

stlo0=x;
stla0=y;
d0=zeros(nt,length(x));

for ii=1:nt
d0(ii,:) = drr_inter_op(squeeze(dn(ii,:,:)),par,0);
end

figure;drr_imagesc(d0);

%% 3D processing/reconstruction
mla=[1,20];
mlo=[1,20];
%binning
[d3d,x1,y1,mask]=drr_bin3d(d0,stlo0,stla0,20,20,mlo(1),mla(1),mlo(2),mla(2));
[stlo1,stla1]=meshgrid(x1,y1);

%% global processing
flow=0;fhigh=0.5;dt=1;N=4;Niter=10;mode=1;verb=1;eps=0.00001;K=4;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=drr3drecon(d3d,mask,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a);


% figure;drr_imagesc([d3d(:,:),d1(:,:)]);

d2=drr3drecon_otg(d0,stlo0,stla0,20,20,mlo(1),mla(1),mlo(2),mla(2),flow,fhigh,dt,N,K,Niter,eps,verb,mode);

% figure;drr_imagesc([d3d(:,:),d1(:,:),d2(:,:)]);




%% plot the results
t=[0:300-1]*0.004;
x1=t(100);
x2=t(200);
y1=15;
y2=20;

ilon=12;

dtest=squeeze(d(:,ilon,:));
figure('units','normalized','Position',[0.2 0.4 1.0, 1.2],'color','w');
subplot(2,4,1);
drr_wigbh(dtest,stla1,t,5);
ylim([0,22]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('Ground truth','Fontsize',16);
drr_framebox(x1,x2,y1,y2,'r',2);
text(-0.2,23,'a)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,2);
dtest=squeeze(d3d(:,ilon,:));
drr_wigbh(dtest,stla1,t,5);
ylim([0,22]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('Observed','Fontsize',16);
drr_framebox(x1,x2,y1,y2,'r',2);
text(-0.2,23,'b)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,3);
dtest=squeeze(d1(:,ilon,:));
drr_wigbh(dtest,stla1,t,5);
ylim([0,22]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('DRR','Fontsize',16);
drr_framebox(x1,x2,y1,y2,'r',2);
text(-0.2,23,'c)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,4);
dtest=squeeze(d2(:,ilon,:));
drr_wigbh(dtest,stla1,t,5);
ylim([0,22]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('OTG','Fontsize',16);
drr_framebox(x1,x2,y1,y2,'r',2);
text(-0.2,23,'d)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

%% zoomed comparison
subplot(2,4,5);
dtest=squeeze(d(:,ilon,:));
drr_wigbh(dtest,stla1,t,1);

ylim([14,21]);xlim([x1,x2]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('Ground truth','Fontsize',16);
text(0.33,21.4,'e)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');


subplot(2,4,6);
dtest=squeeze(d3d(:,ilon,:));
drr_wigbh(dtest,stla1,t,1);
ylim([14,21]);xlim([x1,x2]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('Observed','Fontsize',16);
text(0.33,21.4,'f)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

% print(gcf,'-depsc','-r400','us_lon1_z.eps');  

subplot(2,4,7);
dtest=squeeze(d1(:,ilon,:));
drr_wigbh(dtest,stla1,t,1);
ylim([14,21]);xlim([x1,x2]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('DRR','Fontsize',16);
text(0.33,21.4,'g)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,8);
dtest=squeeze(d2(:,ilon,:));
drr_wigbh(dtest,stla1,t,1);
ylim([14,21]);xlim([x1,x2]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('OTG','Fontsize',16);
text(0.33,21.4,'h)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');


%% following are all annotations
% Create textbox
print(gcf,'-dpng','-r300','test_matdrr_drr3drecon_otg_syn.png');
print(gcf,'-depsc','-r300','test_matdrr_drr3drecon_otg_syn.eps');


%% another slice
ilon=5;

dtest=squeeze(d(:,ilon,:));
figure('units','normalized','Position',[0.2 0.4 1.0, 1.2],'color','w');
subplot(2,4,1);
drr_wigbh(dtest,stla1,t,5);
ylim([0,22]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('Ground truth','Fontsize',16);
drr_framebox(x1,x2,y1,y2,'r',2);
text(-0.2,23,'a)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,2);
dtest=squeeze(d3d(:,ilon,:));
drr_wigbh(dtest,stla1,t,5);
ylim([0,22]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('Observed','Fontsize',16);
drr_framebox(x1,x2,y1,y2,'r',2);
text(-0.2,23,'b)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,3);
dtest=squeeze(d1(:,ilon,:));
drr_wigbh(dtest,stla1,t,5);
ylim([0,22]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('DRR','Fontsize',16);
drr_framebox(x1,x2,y1,y2,'r',2);
text(-0.2,23,'c)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,4);
dtest=squeeze(d2(:,ilon,:));
drr_wigbh(dtest,stla1,t,5);
ylim([0,22]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('OTG','Fontsize',16);
drr_framebox(x1,x2,y1,y2,'r',2);
text(-0.2,23,'d)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

%% zoomed comparison
subplot(2,4,5);
dtest=squeeze(d(:,ilon,:));
drr_wigbh(dtest,stla1,t,1);

ylim([14,21]);xlim([x1,x2]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('Ground truth','Fontsize',16);
text(0.33,21.4,'e)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');


subplot(2,4,6);
dtest=squeeze(d3d(:,ilon,:));
drr_wigbh(dtest,stla1,t,1);
ylim([14,21]);xlim([x1,x2]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('Observed','Fontsize',16);
text(0.33,21.4,'f)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

% print(gcf,'-depsc','-r400','us_lon1_z.eps');  

subplot(2,4,7);
dtest=squeeze(d1(:,ilon,:));
drr_wigbh(dtest,stla1,t,1);
ylim([14,21]);xlim([x1,x2]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('DRR','Fontsize',16);
text(0.33,21.4,'g)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,8);
dtest=squeeze(d2(:,ilon,:));
drr_wigbh(dtest,stla1,t,1);
ylim([14,21]);xlim([x1,x2]);
ylabel('Y (trace)','Fontsize',16);
xlabel('Time (s)','Fontsize',16);
title('OTG','Fontsize',16);
text(0.33,21.4,'h)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','left');


%% following are all annotations
% Create textbox
print(gcf,'-dpng','-r300','test_matdrr_drr3drecon_otg_syn2.png');
print(gcf,'-depsc','-r300','test_matdrr_drr3drecon_otg_syn2.eps');

