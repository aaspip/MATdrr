function [y,thr] = drr_pthresh(x,sorh,t) 
%DRR_PTHRESH Perform soft or hard thresholding or percentile
%       soft or hard thresholding.  
%  Y = WTHRESH(X,SORH,T) returns soft (if SORH = 's') 
%  or hard (if SORH = 'h') T-thresholding  of the input  
%  vector or matrix X. T is the threshold value. 
% 
%  Y = WTHRESH(X,'s',T) returns Y = SIGN(X).(|X|-T)+, soft  
%  thresholding is shrinkage. 
% 
%  Y = WTHRESH(X,'h',T) returns Y = X.1_(|X|>T), hard 
%  thresholding is cruder. 
% 
%  See also WDEN, WDENCMP, WPDENCMP. 
 
%  M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96. 
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

 
switch sorh 
  case 's' 
    tmp = (abs(x)-t); 
    tmp = (tmp+abs(tmp))/2; 
    y   = sign(x).*tmp; 
  
  case 'h' 
    y   = x.*(abs(x)>t);
    
  case 'ps'
   tmp=reshape(abs(x),1,prod(size(x)));
   t=prctile(tmp,100-t);      
   if nargout==2
      thr=t; 
   end
    tmp = (abs(x)-t); 
    tmp = (tmp+abs(tmp))/2; 
    y   = sign(x).*tmp; 
  case 'ph'    
   tmp=reshape(abs(x),1,prod(size(x)));
   t=prctile(tmp,100-t);       
   if nargout==2
      thr=t; 
   end
    y   = x.*(abs(x)>t);   
  otherwise 
    error('Invalid argument value.') 
end
