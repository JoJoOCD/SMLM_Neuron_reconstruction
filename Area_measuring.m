clearvars;
addpath('.\helpers')
which fit -all
% load data
tmp = load(['.\ExampleData\recon_final_backup_2D.mat']);
namei = fields(tmp);
obj = tmp.(namei{1});

% render image
pixelsize = 119; % nm 119
zm_xy = 20;
Recon_color_highb = 6; % <-- adjust this for different image: low value for higher contrast
ss_xy = 1; %filter image
imsz = 256;
colorim = renderim_2D(imsz, zm_xy, obj.loc_x_f/pixelsize, obj.loc_y_f/pixelsize, ss_xy, Recon_color_highb);
%%
% select the line you want to calculate ending with right click
%close all;
h = dipshow(colorim);diptruesize(10);
[xpt,ypt] = getline(h);
hold on;
plot(xpt, ypt, 'r')

% convert pt to nm
ynm = xpt*pixelsize/zm_xy;
xnm = ypt*pixelsize/zm_xy;
sprintf('x (nm): %f', xnm)
sprintf('x1 (nm): %f', xnm(1))
sprintf('x2 (nm): %f', xnm(2))
sprintf('y (nm): %f', ynm)

% selected line length
d = diff([xnm(:) ynm(:)]);
total_length = sum(sqrt(sum(d.*d,2)));
sampling_interval=50;
sprintf('selected line (nm): %f', total_length)

% find line of certain lenght perpendicular to selected line
dist_perp= 500; %1/2 of the full perpendicular line
%Get slope and y int of line AB
angle = atan2(diff(ynm),diff(xnm)); % angle wrong?
sprintf ('angle: %f', angle)
% Choose a point C along line AB at half distance

xval = xnm(1): sampling_interval.*cos(angle):xnm(2);
yval = ynm(1): sampling_interval.*sin(angle):ynm(2);
disp('x1 vs x2')
xnm(1)>xnm(2)
disp('y1 vs y2')
ynm(1)>ynm(2)
disp(cos(angle))
disp(sin(angle))

sampling_sizex=length(xval)
sampling_sizey=length(yval)
%Get slope and all y int of line perpendicular to AB
% perpSlope = -1/slope; 
% perpYint1 = ynm(1) - perpSlope*xnm(1); 
% perpYint2 = ynm(2) - perpSlope*xnm(2); 

%perpYint_all=perpYint1:20:perpYint2 %20 is the sampling interval

% Find the end points of the perpendicular line with length Clen*2

x1 = xval+ones(1,sampling_sizex).*dist_perp.*cos(angle+pi/2);
x2 = xval-ones(1,sampling_sizex).*dist_perp.*cos(angle+pi/2);
y1 = yval+ones(1,sampling_sizey).*dist_perp.*sin(angle+pi/2);
y2 = yval-ones(1,sampling_sizey).*dist_perp.*sin(angle+pi/2);

% Add the perp segment

h = dipshow(colorim);diptruesize(10);
axis equal
hold on
for i =1:sampling_sizex
    
        x=x1(i):-50.*cos(angle+pi/2):x2(i);
    
    %disp('$$$$')
%     disp(x1(i))
%     disp(x2(i))
   
        y=y1(i):-50.*sin(angle+pi/2):y2(i);
  
    %disp('^^^^')
%     disp(y1(i))
%     disp(y2(i))
    diff(xnm)/diff(ynm).*diff(x)/diff(y)
    plot(y/pixelsize*zm_xy,x/pixelsize*zm_xy,'r');
end
plot (xpt,ypt,'g');

figure();
axis equal
hold on
for i =1:sampling_sizex
    x=x1(i):-50.*cos(angle+pi/2):x2(i);
    %disp('$$$$')
%     disp(x1(i))
%     disp(x2(i))
    y=y1(i):-50.*sin(angle+pi/2):y2(i);
    %disp('^^^^')
%     disp(y1(i))
%     disp(y2(i))
    diff(xnm)/diff(ynm).*diff(x)/diff(y)
    plot(x,y,'r');%/pixelsize*zm_xy
end

plot (xnm,ynm,'g')


disp('/////')

% crop out the selected line
rotate_angle = atan2(diff(xnm),diff(ynm)); % to rotate counter clockwise
TM = [cos(rotate_angle) -sin(rotate_angle);sin(rotate_angle) cos(rotate_angle)];% rotation matrix 2D
xp = mean(xnm); % pivot for rotation
yp = mean(ynm); % pivot for rotation
loc0 = cat(2, cat(1, obj.loc_x_f, xnm), cat(1, obj.loc_y_f, ynm))';
center = cat(1, repelem(xp, size(loc0,2)), repelem(yp, size(loc0,2)));
loc0_c = loc0 - center;
loc1 = TM*loc0_c;

h = figure; h.Position = [16          42        1166         954];
scatter(loc1(2,:),loc1(1,:),1,'b.');axis equal;set(gca, 'YDir','reverse');%set(gca, 'XDir','reverse');
hold on;
ywidth = 20; %nm <-- change this to select line width for FWHM calculation
masksubx = loc1(2,:) > min(loc1(2,end-1:end)) & loc1(2,:) < max(loc1(2,end-1:end));
masksuby = loc1(1,:) > mean(loc1(1,end-1:end))-ywidth/2 & loc1(1,:) < mean(loc1(1,end-1:end))+ywidth/2;
masksub = masksubx & masksuby;
loc2 = loc1(:,masksub);
scatter(loc2(2,:),loc2(1,:),1,'r.');axis equal;
title('crop out the selected line (red) in rotated image (blue)')

% render selected line to Gaussian blurred image
reconx = (loc2(1,1:end-2)-min(loc2(1,1:end-2)))./pixelsize;
recony = (loc2(2,1:end-2)-min(loc2(2,1:end-2)))./pixelsize;
subsz = double(ceil(max(max(reconx), max(recony))));
im_line = renderim_2D(subsz, zm_xy, reconx',recony', ss_xy*4, Recon_color_highb);dipshow(im_line)
figure;scatter(recony,reconx,1,'b.');axis equal tight;set(gca, 'YDir','reverse');
title('The selected line to calculate FWHM')

% show the FWHM
curv_FWHM = sum(double(im_line),1);
data_FWHM = curv_FWHM./max(curv_FWHM);

sx_FWHM = linspace(-subsz/2-1,subsz/2,numel(data_FWHM))*pixelsize;
f = fit(sx_FWHM',data_FWHM','gauss2');
h=figure;h.Position=[555         558        1210         350];ha = axes;
bar(sx_FWHM,data_FWHM,0.6,'hist')
hold on
plot(f,'r');

%plot(linspace(f.b1-1.665*f.c1/2, f.b1+1.665*f.c1/2,100), ones(1,100).*f.a1/2)% '-m'

index1 = find(data_FWHM >= 0.5, 1, 'first');

% Find where the data last rises above half the max.
index2 = find(data_FWHM >= 0.5, 1, 'last');
fwhmx = sx_FWHM(index2) - sx_FWHM(index1);
sprintf('cross section diameter (nm): %f', fwhmx)

plot(linspace(sx_FWHM(index1),sx_FWHM(index2)),ones(1,100).*0.5)

% disp(index1)
% disp(index2)

ha.XLabel.String='distance (nm)';
ha.YLabel.String='relative intensity';
ha.FontSize=20;
ha.XLim = [f.b1-5*f.c1, f.b1+5*f.c1];
%text(f.b1+f.c1,f.a1/2,['FWHM',': ',num2str(1.665*f.c1,'%.1f'),' nm'],'fontsize',10) %FWHW gaussian fit
text(f.b1+f.c1,f.a1/2,['FWHM',': ',num2str(fwhmx,'%.1f'),' nm'],'fontsize',10)

legend('data', 'fitted curve', 'FWHM line')