function [ D1 ] = drr_fk_pocs(D,D0,mask,thr,niter,type)
%DRR FK POCS: FK POCS interpolation (POCS + IST + WPOCS)
%  IN   D:   	intput observed data 
%       D0:   	initial data
%               mask: masking operator
%               thr: percentile threshold
%               niter: number of iterations
%               type: POCS or IST or WPOCS
%               'pocs': POCS  (Abma, 2006)
%               'ist': IST (Chen et al., 2015, Seismic data interpolation using nonlinear shaping regularization, Journal of Seismic Exploration, 24, 327-342)
%               'wpocs': Weighted POCS (Gao et al., 2012, Convergence improvement and noise attenuation
%               considerations for beyond alias projection onto convex sets reconstruction)
%      
%  OUT   D1:  	output data
% 
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
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

