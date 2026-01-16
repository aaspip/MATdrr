% Demonstration script for 
% State-of-the-art deblending via the localized damped rank-reduction method (LDRR)
% 
% MIT License
% 
% Copyright (C) 2023 Yangkang Chen
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
% References:   
% 
% Chen et al., 2022, 3D seismic diffraction separation and imaging using the local rank-reduction method, IEEE Transactions on Geoscience and Remote Sensing, 60, 4507110.
% Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, 85, V497–V506.
% Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
% Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
% Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
% Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
% Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
% Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
% Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
% Chen et al., 2023, DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology, Computers & Geosciences, 180, 105440.


clc;clear;close all;

%% Please change the directory path
% requiring the DRR package
% https://github.com/chenyk1990/MATdrr
addpath(genpath('~/MATdrr'));
addpath(genpath('~/deblending-time/subroutines'));

%% please download data from https://drive.google.com/file/d/1ge0Mn_SB4LUsVgOBvATh0iISwGQahKh4/view?usp=sharing
load yc_fieldsr.mat
%% in this dataset
%there are two sources data3d(:,:,1:60) and data3d(:,:,61:120)
%each source contains 120 shots
%there are 60 receivers
%

d1=data3d(:,:,1);
d2=data3d(:,:,60);
figure;
subplot(1,2,1);dbt_imagesc(d1);
subplot(1,2,2);dbt_imagesc(d2);

h1=1:120;
h2=1:120;

dt=0.004;
t=[0:1500-1]*dt;
nt=1500;
nx=120;
%% apply dbt_dithering
randn('state',202122);
shift1=floor(0.1*randn(1,size(d1,2))/dt);   % shift of data1 to data2
shift2=-shift1;                             % shift of data2 to data1

d1shift=dbt_dither(d1,shift1);
d2shift=dbt_dither(d2,shift2);

figure;
subplot(1,2,1);imagesc(h1,t,d1shift);
subplot(1,2,2);imagesc(h2,t,d2shift);

%% blend
d1b=d1+d2shift;
d2b=d2+d1shift;

%% Correct time
del=shift2;

%% mask1
dd=[fliplr(d1),d1];
d1m=dbt_mutter(dd,120,42,505);
% figure;dbt_imagesc([dd,d1m,dd-d1m]);
mask1=ones(size(dd));
mask1=dbt_mutter(mask1,120,42,505);
mask1=mask1(:,121:end);
%mask2
mask2=ones(size(d2));
mask2=dbt_mutter(mask2,60,40,310);

% figure;subplot(1,2,1);imagesc(mask1);subplot(1,2,2);imagesc(mask2);


D1=zeros(nt,nx);
D2=zeros(nt,nx);
for iter=1:10
    fprintf('\n Iter %d \n',iter);
    D1T=dbt_dither(D1,-del);
    D2T=dbt_dither(D2,del);
    D1u = D1 + 0.5*(d1b-(D1+D2T));    % updated model
    D2u = D2 + 0.5*(d2b-(D1T+D2));    % updated model
    D1u=D1u.*mask1;
    D2u=D2u.*mask2;
    D1=drr3d_win(D1u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);
    D2=drr3d_win(D2u,0,80,0.004,2,4,0,100,20,1,0.5,0.5,0.5);
    D1=D1.*mask1;
    D2=D2.*mask2;
    
    D1=D1+(d1b-D1).*(1-mask2).*mask1;
    D2=D2+(d2b-D2).*(1-mask1).*mask2;
    
    snr2(iter)=10*log10(sum(sum(d1.*d1))/sum(sum((d1-D1).*(d1-D1))));
    snr22(iter)=10*log10(sum(sum(d2.*d2))/sum(sum((d2-D2).*(d2-D2))));
    fprintf('iter1=%d,SNR=%g,SNR2=%g\n',iter,snr2(iter),snr22(iter));
end

D33=D1+(d1b-D1).*(1-mask2).*mask1;
D44=D2+(d2b-D2).*(1-mask1).*mask2;



%% Figure 6
[n1,n2]=size(d1);
ngap=10;
comp1=[d1,zeros(n1,ngap),d1b,zeros(n1,ngap),D33,zeros(n1,ngap),d1b-D33,zeros(n1,ngap),d1-D33]; 
comp2=[d2,zeros(n1,ngap),d2b,zeros(n1,ngap),D44,zeros(n1,ngap),d2b-D44,zeros(n1,ngap),d2-D44]; 
t=[0:1500-1]*0.004;
x=1:size(comp1,2);
xts1=[30,60,90];
xts2=xts1+ngap+nx;
xts3=xts1+ngap*2+nx*2;
xts4=xts1+ngap*3+nx*3;
xts5=xts1+ngap*4+nx*4;
xts=[xts1,xts2,xts3,xts4,xts5];
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(4,1,1:2);dbt_imagesc(comp1(1:1000,:),98,1,x,t(1:1000));
text(-50,-0.2,'a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
xticks(xts);
set(gca,'xticklabel',{'30','60','90'});
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(nx/2,0.16,'Unblended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx+ngap+nx/2,0.16,'Blended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*2+nx/2,0.16,'Deblended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*3+nx/2,0.16,'Blending noise','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*4+nx/2,0.16,'Deblending error','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

text(nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx+ngap+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*2+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*3+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*4+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

subplot(4,1,3:4);dbt_imagesc(comp2(1:1000,:),98,1,x,t(1:1000));
text(-50,-0.2,'b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','left');
xticks(xts);
set(gca,'xticklabel',{'30','60','90'});
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(nx/2,0.16,'Unblended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx+ngap+nx/2,0.16,'Blended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*2+nx/2,0.16,'Deblended','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*3+nx/2,0.16,'Blending noise','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*4+nx/2,0.16,'Deblending error','color','b','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(nx+ngap+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*2+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*3+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text((nx+ngap)*4+nx/2,4.35,'Shot NO #','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

print(gcf,'-depsc','-r300','test_matdrr_drr3d_win_deblending.eps');
print(gcf,'-dpng','-r300','test_matdrr_drr3d_win_deblending.png');
