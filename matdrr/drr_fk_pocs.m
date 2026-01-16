function [ D1 ] = drr_fk_pocs(D,D0,mask,thr,niter,type)
%DRR FK POCS: FK POCS interpolation (POCS + IST + WPOCS)
% IN   D:   	intput observed data 
%      D0:   	initial data
%              mask: masking operator
%              thr: percentile threshold
%              niter: number of iterations
%              type: POCS or IST or WPOCS
%              'pocs': POCS  (Abma, 2006)
%              'ist': IST (Chen et al., 2015, Seismic data interpolation using nonlinear shaping regularization, Journal of Seismic Exploration, 24, 327-342)
%              'wpocs': Weighted POCS (Gao et al., 2012, Convergence improvement and noise attenuation
%              considerations for beyond alias projection onto convex sets reconstruction)
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

D1=D0;
 switch type 
     case 'pocs'
         for iter=1:niter
            D1=D.*mask+(1-mask).*drr_fkt(D1,'ps',thr); 
         end
         
     case 'ist'
         for iter=1:niter
           Dtmp=D1+D-mask.*D1;  
           D1=drr_fkt(Dtmp,'ps',thr);  
         end

     case 'wpocs'
         a=(niter-(1:niter))/(niter-1);%if this option, WPOCS = IST, you can try, why?
         % a=(niter-(1:niter))/(niter-1);a=a.*a;%
         for iter=1:niter
            Dn=drr_fkt(D1,'ps',thr);
            D1=a(iter)*D.*mask+(1-a(iter))*mask.*Dn+(1-mask).*Dn; 
         end

  otherwise 
    error('Invalid argument value.')      
    
 end
return

