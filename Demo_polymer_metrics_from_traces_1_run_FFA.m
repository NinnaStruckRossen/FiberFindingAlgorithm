clear
close all
delete('./psf_scan.mat')

%% User variables

psf_density = [24 24 100];  % voxel resolution of psf image [x y z] in nm/pixel
data_density = [82 82 200]; % voxel resolution of data image [x y z] in nm/pixel
stop_factor = 0.4;          % lies within [0:1] and can be optimized in stop-factor demo. 0: estimated noise floor, 1: estimated signal
angle_max = 15;             % 1/2 max cone angle change between each step (in degrees, [10:180[).
l_prio = [15 10 5];         % array allowing prioritized search for lenghts

it_max = -1;                % if set to '-1', it_max is default value. max step number per object (effective trace max length)
step_max = -1;              % If value is '-1', the algorithm will use a step_max based on the PSF size.

%% Input measured actual point-spread-function (psf) or theoretical 3D gaussian approximation

FileTif= './Data_ourPSF.tif';
% FileTif= './Data_theoreticalPSF.tif';
% change "Data_ourPSF" to the measured PSF image of your imaging system and
% objective above if you are analysing your own small dataset, or to
% "Data_theoreticalPSF" if your actual measured PSF is unavailable.
myImInfo = imfinfo(FileTif,'tiff');
z_range = 1:length(myImInfo);
x_el_needed = myImInfo(1).Height ;
y_el_needed = myImInfo(1).Width ;
psf = zeros(x_el_needed,y_el_needed,numel(z_range),'double') ;
for i=z_range
    psf(:,:,i)=double(imread(FileTif,'Index',i));
end


%% Input data

FileTif= './SampleData_CollagenTypeI_1mgpermL_23C.tif';
% change "SampleData_CollagenTypeI_1mgpermL_23C" to your acquired data above
myImInfo = imfinfo(FileTif,'tiff');
z_range = 1:length(myImInfo);
x_el_needed = myImInfo(1).Height ;
y_el_needed = myImInfo(1).Width ;
data = zeros(x_el_needed,y_el_needed,numel(z_range),'double') ;
for i=z_range
    data(:,:,i)=double(imread(FileTif,'Index',i));
end

%% Finding traces

[corr,traces,data_applied,corr_backup] = fibtracer(data,psf,data_density,psf_density, ...
    it_max,stop_factor,angle_max,step_max,l_prio);

%% Figure

cmap_fig_blues =  [0,0,0.5;0,0,0.6875;0,0,0.8750;0,0.0625,1;0,0.2500,1;0,0.4375,1];

figure(1)

subplot(2,2,1)
title('Rotatable 3D traces overlaid on maximum intensity projection')
hold on;
imagesc(max(data_applied(:,:,:),[],3))
for k=1:length(traces)
    L=traces{k};
    plot3(L(:,2),L(:,1),L(:,3),'-','Color',cmap_fig_blues(mod(k,length(cmap_fig_blues))+1,:),'LineWidth',3);
end
colormap(gray)
hold off
daspect([1 1 1])
xlim([1 size(data_applied,2)])
ylim([1 size(data_applied,1)])
%%%%%%%%%%%%%%%%%

subplot(2,2,3)
title('Rotatable 3D traces overlaid on residuals')
hold on;
imagesc(min(corr_backup,[],3))
for k=1:length(traces)
    L=traces{k};
    plot3(L(:,2),L(:,1),L(:,3),'-','Color',cmap_fig_blues(mod(k,length(cmap_fig_blues))+1,:),'LineWidth',3);
end
colormap(gray)
hold off
daspect([1 1 1])
xlim([1 size(data_applied,2)])
ylim([1 size(data_applied,1)])

%%%%%%%%%%%%%%%%%

subplot(2,2,[2,4])
title('Rotatable 3D volume with traces')
patch(isosurface(data_applied,mean(data_applied(:))+25), ...
    'FaceColor',[1 1 1], ...
    'edgealpha',0,...
    'FaceAlpha',0.4);
daspect([1 1 1])
xlim([1 size(data_applied,2)])
ylim([1 size(data_applied,1)])
zlim([1 size(data_applied,3)])
hold on
for k=1:length(traces)
    L=traces{k};
    plot3(L(:,2),L(:,1),L(:,3),'-','Color',cmap_fig_blues(mod(k,length(cmap_fig_blues))+1,:),'LineWidth',3);
    view([20 20])
    camlight
end
set(gca,'color',[0 0 0])

