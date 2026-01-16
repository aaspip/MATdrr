function  drr_wigbh(din,x,z,scale)
%WIGBH wigb horizontally
%din: size nz*nx
%
% Nature of clip: the plotting amplitude correspoding to the clip value is 1 
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
% Example:
% d=levents;figure;wigbh(d,1:50,0.004*[0:500],2);


[nz,nx]=size(din);
vmax=max(din(:));
if nargin==1
   x=1:nx;
   z=1:nz; 
   clip=vmax;
end
if nargin==3
   clip=vmax;
end

din=din*scale;


for i = 1 : nx
    plot(z, din(:,i)+ x(i),'k','LineWidth',1.5); hold on
end
% set(gca,'ytick',(25:5:95)*2);
% set(gca,'yticklabel',[25:5:95]);
set(gca,'FontSize',20);
%set(gca,'YDir','reverse');
set(gca,'linewidth',1.5);
axis tight
% xlabel('time/s');
% xlim([20 60]);
% ylabel('gcarc');
set(gcf,'color','white');

end

