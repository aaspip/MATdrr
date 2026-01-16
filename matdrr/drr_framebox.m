function [] = drr_framebox(x1,x2,y1,y2,c,lw)
%
% for drawing a frame box
% 
% Jan, 5, 2023, By Yangkang Chen
% 
% x1,x2,y1,y2: intuitive
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

if nargin==4
c='r';
lw=2;
end

hold on;


hold on;
plot([x1,x2],[y1,y1],'-','color',c,'linewidth',lw);
plot([x1,x2],[y2,y2],'-','color',c,'linewidth',lw);
plot([x1,x1],[y1,y2],'-','color',c,'linewidth',lw);
plot([x2,x2],[y1,y2],'-','color',c,'linewidth',lw);


return