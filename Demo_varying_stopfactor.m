clear
close all
delete('./psf_scan.mat')

%% User variables

psf_density = [24 24 100];  % voxel resolution of psf image [x y z] in nm/pixel
data_density = [82 82 200]; % voxel resolution of data image [x y z] in nm/pixel
stop_factor = 0.4;          % lies within [0:1] with % 0: estimated noise floor, 1: estimated signal
% This stop_factor variable is being varied below [0.1:0.1:0.9] in this stop-factor demo to optimize to the noise-level of the acquired image.
angle_max = 15;             % 1/2 max cone angle change between each step (in degrees, [10:180[).
l_prio = [15 10 5];         % array allowing priritized search for lenghts

it_max = -1;                % if set to '-1', it_max is default value. max step number per object (effective trace max length)
step_max = -1;              % -1;% If value is '-1', the algorithm will use a step_max based on the PSF size.

%% Input Point-Spread-Function (psf)

FileTif= './Data_ourPSF.tif';
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

% choose a dataset of the three below or substitute with your own:
FileTif= './SampleData_small.tif';
% FileTif= './SampleData_CollagenTypeI_1mgpermL_23C.tif';
% FileTif= './SampleData_CollagenTypeI_1mgpermL_37C.tif';
myImInfo = imfinfo(FileTif,'tiff');
z_range = 1:length(myImInfo);
x_el_needed = myImInfo(1).Height ;
y_el_needed = myImInfo(1).Width ;
data = zeros(x_el_needed,y_el_needed,numel(z_range),'double') ;
for i=z_range
    data(:,:,i)=double(imread(FileTif,'Index',i));
end

%% Finding traces

for n = 1:9
    stop_factor = n*0.1 ;
    [corr,traces,data_applied,corr_backup] = fibtracer(data,psf,data_density,psf_density, ...
        it_max,stop_factor,angle_max,step_max,l_prio);
    traces_master{n} = traces ;
end

%% Figure

cmap_fig = jet;
cmap_fig_reds = [1,0.5000,0;1,0.4375,0;1,0.3750,0;1,0.3125,0;1,0.2500,0;...
    1,0.1875,0;1,0.1250,0;1,0.0625,0;1,0,0;0.9375,0,0;...
    0.8750,0,0;0.8125,0,0;0.7500,0,0;0.6875,0,0;...
    0.6250,0,0;0.5625,0,0;0.5000,0,0];
cmap_fig_blues =  [0,0,0.5;0,0,0.5625;0,0,0.6250;0,0,0.6875;0,0,0.7500;...
    0,0,0.8125;0,0,0.8750;0,0,0.9375;0,0,1;0,0.0625,1;...
    0,0.1250,1;0,0.1875,1;0,0.2500,1;0,0.3125,1;...
    0,0.3750,1;0,0.4375,1;0,0.5000,1];

figure(1)
for n=1:9
    subplot(3,3,n)
    title(strcat('Stopfactor: 0.',num2str(n)))
    traces = traces_master{n} ;
    hold on;
    imagesc(max(data_applied(:,:,:),[],3))
    for k=1:length(traces)
        L=traces{k};
        % choose a color scheme of the three below:
        plot3(L(:,2),L(:,1),L(:,3),'-','Color',cmap_fig(round(rand*63+1),:),'LineWidth',1); % for random colors
        % plot3(L(:,2),L(:,1),L(:,3),'-','Color',cmap_fig_reds(mod(k,length(cmap_fig_reds))+1,:),'LineWidth',1); % for red colors
        % plot3(L(:,2),L(:,1),L(:,3),'-','Color',cmap_fig_blues(mod(k,length(cmap_fig_blues))+1,:),'LineWidth',1); % for blue colors
    end
    colormap(gray)
    hold off
    daspect([1 1 1])
    xlim([1 size(data_applied,2)])
    ylim([1 size(data_applied,1)])
end
