function [colorim] = renderim(imsz, zm_xy, x, y, c, ss_xy, Recon_color_highb)
[rch, gch, bch]=srhist_color(imsz,zm_xy,x,y,c,64);
rchsm=gaussf(rch,ss_xy);
gchsm=gaussf(gch,ss_xy);
bchsm=gaussf(bch,ss_xy);

rchsmst=imstretch_linear(rchsm,0,Recon_color_highb,0,255);
gchsmst=imstretch_linear(gchsm,0,Recon_color_highb,0,255);
bchsmst=imstretch_linear(bchsm,0,Recon_color_highb,0,255);

colorim=joinchannels('RGB',rchsmst,gchsmst,bchsmst);
