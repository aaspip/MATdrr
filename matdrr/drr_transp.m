function [dout]=db_transp(din,plane)
% db_transp: Transpose two axes in a dataset
% by Yangkang Chen, Dec 18, 2019
% Modified on Jan, 2020
% 
% INPUT
% din: input dataset
% plane: Two-digit number with axes to transpose. The default is 12
% OUTPUT
% dout: output dataset
%
% DEMO:
% a=magic(3);b=reshape(a,3,1,3);c=yc_transp(b,23);norm(a-c)
% 
% MIT License
% 
% Copyright (C) 2020 Yangkang Chen
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

if nargin==1
    plane=12;
end

[n1,n2,n3,n4,n5]=size(din);

switch plane
    
    case 12
        dout=zeros(n2,n1,n3,n4,n5);
        
        for i5=1:n5
            for i4=1:n4
                for i3=1:n3
                    dout(:,:,i3,i4,i5)=din(:,:,i3,i4,i5).';
                end
            end
        end
        
    case 23
        dout=zeros(n1,n3,n2,n4,n5);
        for i5=1:n5
            for i4=1:n4
                for i1=1:n1
                    dout(i1,:,:,i4,i5)=squeeze(din(i1,:,:,i4,i5)).';
                end
            end
        end
    case 13
        dout=zeros(n3,n2,n1,n4,n5);
        for i5=1:n5
            for i4=1:n4
                for i2=1:n2
                    dout(:,i2,:,i4,i5)=squeeze(din(:,i2,:,i4,i5)).';
                end
            end
        end
        
     case 14
        dout=zeros(n4,n2,n3,n1,n5);
        for i5=1:n5
            for i3=1:n3
                for i2=1:n2
                    dout(:,i2,i3,:,i5)=squeeze(din(:,i2,i3,:,i5)).';
                end
            end
        end
        
     case 15
        dout=zeros(n5,n2,n3,n4,n1);
        for i4=1:n4
            for i3=1:n3
                for i2=1:n2
                    dout(:,i2,i3,i4,:)=squeeze(din(:,i2,i3,i4,:)).';
                end
            end
        end        
        
    otherwise
        error('Invalid argument value.');
end

return