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

sampling_interval=20;
dist_perp= 500; %1/2 of the full perpendicular line
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

%find line of certain lenght perpendicular to selected line
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
    diff(xnm)/diff(ynm).*diff(x)/diff(y);
    plot(y/pixelsize*zm_xy,x/pixelsize*zm_xy,'r');
end
plot (xpt,ypt,'g');

% figure();
% axis equal
% hold on

sampling_diameter=[];
for i =1:sampling_sizex
    x_prep=[x1(i),x2(i)];
    y_prep=[y1(i),y2(i)];
    diff(xnm)/diff(ynm).*diff(x)/diff(y);
    sprintf('this is %f iteration', i);
    %plot(x_prep,y_prep,'r');%/pixelsize*zm_xy
    value=cut_out_square(x_prep',y_prep',obj,i,pixelsize,zm_xy,ss_xy,Recon_color_highb,sampling_interval);
    sampling_diameter=[sampling_diameter;value];
end

figure()
total_length_x=[1:sampling_interval:total_length];
plot (total_length_x, sampling_diameter/sampling_interval)
disp('/////')

%plot (xnm,ynm,'g');
area=sum(sampling_diameter);
area_um=area/(10^6);

sprintf('selected line (nm): %f', total_length)
sprintf('selected area (nm^2): %f',area)
sprintf('selected area (um^2): %f',area_um)

