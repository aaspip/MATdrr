function [ D1 ] = drr3drecon_otg(D,x,y,nx,ny,ox,oy,mx,my,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a)
% DRR3DRECON_OTG: OTG DRR for 3D seismic reconstruction 
%
% IN   D:   	 intput 2D data
%      x:     input x coordinates
%      y:     input y coordinates
%      nx:    input number of binned x points
%      ny:    input number of binned y points
%      ox:    min of x
%      oy:    min of y
%      mx:    max of x
%      my:    max of y
%      flow:   processing frequency range (lower)
%      fhigh:  processing frequency range (higher)
%      dt:     temporal sampling interval
%      N:      number of singular value to be preserved
%      K:     damping factor (default: 4)
%      Niter:  number of maximum iteration
%      eps:    tolerence (||S(n)-S(n-1)||_F<eps)
%      verb:   verbosity flag (default: 0)
%      mode:   mode=1: denoising and reconstruction
%              mode=0: reconstruction only
%      a:      weight vector
%
% OUT  D1:  	output data
% 
% MIT License
% 
% Copyright (C) 2021 Yangkang Chen
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
% [1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
% [2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
% [3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
% [4] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
% [5] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
% [6] Chen et al., 2023, DRR: an open-source multi-platform package for the damped rank-reduction method and its applications in seismology, Computers & Geosciences, 180, 105440.

if nargin==0
    error('Input data must be provided!');
end

if nargin==1
    error('Samping mask should be given');
end

if nargin==2
    flow=1;
    fhigh=124;
    dt=0.004;
    N=1;
    K=4;
    Niter=30;
    eps=0.00001;
    verb=0;
    mode=1;
end;

if mode==0;
    a=ones(1,Niter);
end


nt=size(D,1);
D1=zeros(nt,nx,ny);

nf=2^nextpow2(nt);

% Transform into F-X domain
DATA_FX=fft(D,nf,1);
DATA_FX0=zeros(nf,nx,ny);

% First and last nts of the DFT.
ilow  = floor(flow*dt*nf)+1;

if ilow<1;
    ilow=1;
end;

ihigh = floor(fhigh*dt*nf)+1;

if ihigh > floor(nf/2)+1;
    ihigh=floor(nf/2)+1;
end

lx=floor(nx/2)+1;
lxx=nx-lx+1;
ly=floor(ny/2)+1;
lyy=ny-ly+1;
M=zeros(lx*ly,lxx*lyy);

% construct par
par.x=x;
par.y=y;
par.nx=nx;
par.ny=ny;
par.ox=ox;
par.oy=oy;
par.mx=mx;
par.my=my;
s=0.5*ones(Niter,1);
% main loop
for k=ilow:ihigh
    S_obs=squeeze(DATA_FX(k,:)).'; %1D vector  
    Sn_1=zeros(nx,ny);
%    Sn_1=drr_inter_op(DATA_FX(k,:),par,-1);
    for iter=1:Niter
        
%        size(Sn_1)
%        size(S_obs)
%        size(inter_op(Sn_1,par,1))
        Sn=Sn_1-s(iter)*drr_inter_op(drr_inter_op(Sn_1,par,-1)-S_obs,par,1);
        
        
        M=P_H(Sn,lx,ly);
        M=P_RD(M,N,K);
        Sn=P_A(M,nx,ny,lx,ly);
        
%        Sn=a(iter)*S_obs+(1-a(iter))*mask.*Sn+(1-mask).*Sn;
        if norm(Sn-Sn_1,'fro')<eps
            break;
        end
        Sn_1=Sn;
    end
    
    for j=1:ny
        DATA_FX0(k,:,j) = DATA_FX0(k,:,j)+reshape(Sn(:,j),1,nx);
    end
    
    if(mod(k,5)==0 && verb==1)
        fprintf( 'F %d is done!\n\n',k);
    end
    
end

% Honor symmetries
for k=nf/2+2:nf
    DATA_FX0(k,:,:) = conj(DATA_FX0(nf-k+2,:,:));
end

% Back to TX (the output)
D1=real(ifft(DATA_FX0,[],1));
D1=D1(1:nt,:,:);

return

function [dout]=P_H(din,lx,ly)
% forming block Hankel matrix
% size(din)
[nx,ny]=size(din);
lxx=nx-lx+1;
lyy=ny-ly+1;

for j=1:ny
    r=hankel(din(1:lx,j),[din(lx:nx,j)]);
    if j<ly
        for id=1:j
            dout(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx) = r;
        end
    else
        for id=1:(ny-j+1)
            dout((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx)=r;
        end
    end
end
return

function [dout]=P_RD(din,N,K)
% Rank reduction on the block Hankel matrix


%     [U,D,V]=svds(din,N); % a little bit slower for small matrix
%     dout=U*D*V';
% %
[U,D,V]=svd(din);
for j=1:N
    D(j,j)=D(j,j)*(1-D(N+1,N+1)^K/(D(j,j)^K+0.000000000000001));
end

dout=U(:,1:N)*D(1:N,1:N)*(V(:,1:N)');

return

function [dout]=P_A(din,nx,ny,lx,ly)
% Averaging the block Hankel matrix to output the result
lxx=nx-lx+1;
lyy=ny-ly+1;
dout=zeros(nx,ny);

for j=1:ny
    if j<ly
        for id=1:j
            dout(:,j) =dout(:,j)+ ave_antid(din(1+(j-1)*lx-(id-1)*lx:j*lx-(id-1)*lx,1+(id-1)*lxx:lxx+(id-1)*lxx))/j;
        end
    else
        for id=1:(ny-j+1)
            dout(:,j) =dout(:,j)+ ave_antid(din((ly-1)*lx+1-(id-1)*lx:ly*lx-(id-1)*lx,(j-ly)*lxx+1+(id-1)*lxx:(j-ly+1)*lxx+(id-1)*lxx))/(ny-j+1);
        end
    end
end
return


function [dout] =ave_antid(din);
% averaging along antidiagonals
[n1,n2]=size(din);
nout=n1+n2-1;
dout=zeros(nout,1);
for i=1:nout
    if i<n1
        for id=1:i
            dout(i)=dout(i) + din(i-(id-1),1+(id-1))/i;
        end
    else
        for id=1:nout+1-i
            dout(i)=dout(i) + din(n1-(id-1),1+(i-n1)+(id-1))/(nout+1-i);
        end
    end
end
return

