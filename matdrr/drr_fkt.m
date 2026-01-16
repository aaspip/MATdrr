function [ D1 ] = drr_fkt(D,sorh,t)
%DRR_FKT: FK domain thrsholding (any dimension)
% IN   D:   	intput data 
%      sorh:   soft or hard thresholding (s,h,ps,ph)
%              s: soft
%              h: hard
%              ps: percentile soft
%              ph: percentile hard
%      t:      threshold value
%     
% OUT   D1:  	output data
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

if nargin==1 % for quick denoising test
   sorh='ps';
   t=5;
end
D1=ifftn(drr_pthresh(fftn(D),sorh,t));
% [n1,n2,n3]=size(D);
% n1h=floor(n1/2)+1;
% 
% Dfft=reshape(fftn(D),n1,n2*n3);
% Dtmp=reshape(Dfft(1:n1h,:),n1h,n2,n3); % thresholding only on the half space
% %Dfft(1:n1h,:)=reshape(pthresh(Dtmp,sorh,t),n1h,n2*n3);
% 
% Dfft(1:n1h,:)=reshape(Dtmp,n1h,n2*n3);
% 
% %% honor conjugate symmetry property
% if mod(n1,2)==0 % even
% Dfft(n1h+1:end,:)=conj(flipud(Dfft(2:n1h-1,:)));
% else            % odd
% Dfft(n1h+1:end,:)=conj(flipud(Dfft(2:n1h,:)));    
% end
% Dfft=reshape(Dfft,n1,n2,n3);
% %D1=ifftn(Dfft);
% D1=Dfft;
return
