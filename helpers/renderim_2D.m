function [colorim] = renderim_2D(sz, zm, reconx, recony, ss_xy, Recon_color_highb)
srim = SRreconstructhist(sz,zm,reconx,recony,0);
srim = double(srim);
smoothim = imgaussfilt(srim,[ss_xy ss_xy]);  %    
imstr = imstretch_linear(smoothim,0,Recon_color_highb,0,255);% 修正，将亮度过高和过低的值分别置255和0
colorim = uint8(imstr);
