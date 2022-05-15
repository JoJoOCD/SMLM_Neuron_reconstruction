clearvars;
addpath('.\helpers')

% load data
tmp = load(['.\ExampleData\recon_final_backup_3D.mat']);
namei = fields(tmp);
obj = tmp.(namei{1});

% render image
pixelsize = 119; % nm
zm_xy = 20;
Recon_color_highb = 3; % low value for higher contrast
ss_xy = 1;
imsz = 256;
colorim = renderim(imsz, zm_xy, obj.loc_x_f/pixelsize, obj.loc_y_f/pixelsize, obj.loc_z_f, ss_xy, Recon_color_highb);
%% calculate structure size
close all;
h = dipshow(colorim);diptruesize(10);
[xpt,ypt] = getline(h);
hold on;
plot(xpt, ypt, 'r')

% convert unit to nm
ynm = xpt*pixelsize/zm_xy;
xnm = ypt*pixelsize/zm_xy;
%
xp = mean(xnm);
yp = mean(ynm);
% Take localizations on that line
rotate_angle = atan2(diff(xnm),diff(ynm)); % to rotate counter clockwise
TM = [cos(rotate_angle) -sin(rotate_angle);sin(rotate_angle) cos(rotate_angle)];% rotation matrix 2D
loc0 = cat(2, cat(1, obj.loc_x_f, xnm), cat(1, obj.loc_y_f, ynm))';
center = cat(1, repelem(xp, size(loc0,2)), repelem(yp, size(loc0,2)));
loc0_c = loc0 - center;
loc1 = TM*loc0_c;

% crop out the selected line
h = figure; h.Position = [16          42        1166         954];
scatter(loc1(2,:),loc1(1,:),1,'c.');axis equal;set(gca, 'YDir','reverse');%set(gca, 'XDir','reverse');
hold on;
ywidth = 20; %nm <-- change this to select line width for FWHM calculation
masksubx = loc1(2,:) > min(loc1(2,end-1:end)) & loc1(2,:) < max(loc1(2,end-1:end));
masksuby = loc1(1,:) > mean(loc1(1,end-1:end))-ywidth/2 & loc1(1,:) < mean(loc1(1,end-1:end))+ywidth/2;
masksub = masksubx & masksuby;
loc2 = loc1(:,masksub);
scatter(loc2(2,:),loc2(1,:),1,'r.');axis equal;
title('Rotated image (cyan) to crop out the selected line (red)')

% render selected line to Gaussian blurred image
reconx = (loc2(1,:)-min(loc2(1,:)))./pixelsize;
recony = (loc2(2,:)-min(loc2(2,:)))./pixelsize;
subsz = double(ceil(max(max(reconx), max(recony))));
im_line = renderim_2D(subsz, zm_xy, reconx',recony', ss_xy*4, Recon_color_highb);dipshow(im_line);
figure;scatter(recony,reconx,1,'b.');axis equal tight;set(gca, 'YDir','reverse');
title('The selected line to calculate FWHM_x_y')

% show the FWHM
curv_FWHM = sum(double(im_line),1);
data_FWHM = curv_FWHM./max(curv_FWHM);

sx_FWHM = linspace(-subsz/2-1,subsz/2,numel(data_FWHM))*pixelsize;
f = fit(sx_FWHM',data_FWHM','gauss1');
h2=figure;h2.Position=[555         558        1210         350];ha = axes;
bar(sx_FWHM,data_FWHM,0.6,'hist')
hold on
plot(f,'r');
plot(linspace(f.b1-1.665*f.c1/2, f.b1+1.665*f.c1/2,100), ones(1,100).*f.a1/2, '-m')
ha.XLabel.String='lateral distance (nm)';
ha.YLabel.String='relative intensity';
ha.FontSize=20;
ha.XLim = [f.b1-5*f.c1, f.b1+5*f.c1];
text(f.b1+f.c1,f.a1/2,['FWHM',': ',num2str(1.665*f.c1,'%.1f'),' nm'],'fontsize',10)
legend('data', 'fitted curve', 'FWHM line')
A=1.665*f.c1
% crop the selected line with wider linewidth for axial calculation
ywidth_crosssec = 250; %nm
masksubx2 = loc1(2,:) > min(loc1(2,end-1:end)) & loc1(2,:) < max(loc1(2,end-1:end));
masksuby2 = loc1(1,:) > mean(loc1(1,end-1:end))-ywidth_crosssec/2 & loc1(1,:) < mean(loc1(1,end-1:end))+ywidth_crosssec/2;
masksub2 = masksubx2 & masksuby2;
loc3 = loc1(:,masksub2);
h3 = figure; h3.Position = [16          42        1166         954];
scatter(loc1(2,:),loc1(1,:),1,'c.');axis equal;set(gca, 'YDir','reverse');%set(gca, 'XDir','reverse');
hold on
scatter(loc3(2,:),loc3(1,:),1,'r.');axis equal;
title('Rotated image (cyan) to crop out the selected line (red)')

% reconx = (loc3(1,:)-min(loc3(1,:)))./pixelsize;
recony2 = (loc3(2,:)-min(loc3(2,:)))./pixelsize;
reconz = obj.loc_z_f(masksub2);
subsz2 = double(ceil(max(max(reconz./pixelsize), max(recony2))));
reconz(1:2) = [min(obj.loc_z_f(:)),max(obj.loc_z_f(:))]';
reconz = reconz - min(reconz);
im_line_z = renderim(subsz2, zm_xy, recony2', reconz'./pixelsize, reconz', ss_xy*4, Recon_color_highb);dipshow(im_line_z)

figure;scatter(recony2,reconz./pixelsize,1,reconz,'.');colormap('jet');axis equal tight;set(gca, 'YDir','reverse');
title('The selected line to calculate FWHM_z');
ylabel('z (pixel)')
xlabel('lateral (pixel)');

% show the FWHM
curv_FWHM_z = sum(sum(double(im_line_z),3),1);
data_FWHM_z = curv_FWHM_z./max(curv_FWHM_z);

sx_FWHM_z = linspace(-subsz2/2,subsz2/2,numel(data_FWHM_z))*pixelsize;
f2 = fit(sx_FWHM_z',data_FWHM_z','gauss1'); 
h4=figure;h4.Position=[555         558        1210         350];ha = axes;
bar(sx_FWHM_z,data_FWHM_z,0.6,'hist')
hold on
plot(f2,'r');
plot(linspace(f2.b1-1.665*f2.c1/2, f2.b1+1.665*f2.c1/2,100), ones(1,100).*f2.a1/2, '-m')
ha.XLabel.String='axial distance (nm)';
ha.YLabel.String='relative intensity';
ha.FontSize=20;
ha.XLim = [f2.b1-5*f2.c1, f2.b1+5*f2.c1];
text(f2.b1+f2.c1,f2.a1/2,['FWHM',': ',num2str(1.665*f2.c1,'%.1f'),' nm'],'fontsize',10)
legend('data', 'fitted curve', 'FWHM', 'Location', 'NortheastOutside')
B=1.665*f2.c1
cross_section_area=pi*A*B/4
display(cross_section_area)