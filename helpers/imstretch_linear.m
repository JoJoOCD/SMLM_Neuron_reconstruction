function imout=imstretch_linear(im,low_in,high_in,low_out,high_out)
% imout=imstretch_linear(im,low_in,high_in,low_out,high_out)
% For each element in 'im', set it to high_out if it is large than high_in,
% to low_out if less than low_in, not changed otherwise.
% ? Fang Huang, 2012, fang.huang@yale.edu
imout=(im-low_in)./(high_in-low_in);
imout(imout>1)=1;
imout(imout<0)=0;
imout=(imout).*(high_out-low_out)+low_out;