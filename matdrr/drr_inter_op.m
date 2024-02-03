function [U] = drr_inter_op(D,par,adj)
%%drr_inter_op: Forward bi-linear interpolation operator
% 
% By Yangkang Chen
% Jan, 2022
% Revised in April, 2024
% 
% INPUT
% 
% if adj = true, D: 1D data, U: 2D data
% if adj = false, D: 2D data, U: 1D data
% 
% par.x=x;
% par.y=y;
% par.nx=nx;
% par.ny=ny;
% par.ox=ox;
% par.oy=oy;
% par.mx=mx;
% par.my=my;

rx=par.x;
ry=par.y;

dx=(par.mx-par.ox)/(par.nx-1);
dy=(par.my-par.oy)/(par.ny-1);
xx=par.ox+[0:par.nx-1]*dx;
yy=par.oy+[0:par.ny-1]*dy;


Nu=length(rx);
Nx=length(xx);
Ny=length(yy);
if adj==1
    U=zeros(Nx,Ny);
else
    U=zeros(Nu,1);
end
for k=1:Nu
    ia=floor((rx(k)-xx(1))/dx)+1;ib=ia+1;
    ja=floor((ry(k)-yy(1))/dy)+1;jb=ja+1;
    if ib>1 & ib<=Nx & jb>1 & jb<=Ny
        t=(rx(k)-xx(ia))/dx;
        u=(ry(k)-yy(ja))/dy;
        if adj==1
            U(ia,ja)=U(ia,ja)+(1-t)*(1-u)*D(k);
            U(ib,ja)=U(ib,ja)+t*(1-u)*D(k);
            U(ib,jb)=U(ib,jb)+t*u*D(k);
            U(ia,jb)=U(ia,jb)+(1-t)*u*D(k);
        else
%            fprintf('size(D)=(%d,%d)\n',size(D,1),size(D,2));
            U(k)=(1-t)*(1-u)*D(ia,ja)+t*(1-u)*D(ib,ja)+...
                t*u*D(ib,jb)+(1-t)*u*D(ia,jb);
        end
   
    end
end
