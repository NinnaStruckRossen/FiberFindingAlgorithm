clear
close all
delete('./psf_scan.mat')

%% User variables

psf_density = [24 24 100];  % voxel resolution of psf image [x y z] in nm/pixel
data_density = [82 82 200]; % voxel resolution of data image [x y z] in nm/pixel
stop_factor = 0.4;          % lies within [0:1] and can be optimized in stop-factor demo. 0: estimated noise floor, 1: estimated signal
% angle_max altered between:
angle_max_narrow = 15;      % 1/2 max cone angle change between each step (in degrees, [10:180[).
angle_max_wide = 45;
angle_max_wider = 90;
% l_array altered between:
l_un_prio = 2 ;             % Un-prioritized search
l_prio = [15 10 5];         % array allowing prioritized search for lenghts
l_high_prio = [30 15 10];   % array allowing higher prioritized search for lenghts

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

%% Finding traces without and with prioritizing length

for a = 1:3
    for n = 1:3
        
        if a == 1
            angle_max = angle_max_wider ;
        elseif a == 2
            angle_max = angle_max_wide;
        else
            angle_max = angle_max_narrow;
        end
        if n == 1
            l_array = l_un_prio ;
        elseif n == 2
            l_array = l_prio;
        else
            l_array = l_high_prio;
        end
        [corr,traces,data_applied,corr_backup] = fibtracer(data,psf,data_density,psf_density, ...
            it_max,stop_factor,angle_max,step_max,l_array);
        traces_master{n+(a-1)*3} = traces ;
    end
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
for a=1:3
    for n=1:3
        subplot(3,3,n+(a-1)*3)
        if n == 1 && a == 1
            title('Un-prioritized search within 180^o cone angle')
        elseif n == 2 && a == 1
            title('Length prioritized search [15 10 5] within 180^o cone angle')
        elseif n == 3  && a == 1
            title('Length highly prioritized search [30 15 10] within 180^o cone angle')
        elseif n == 1 && a == 2
            title('Un-prioritized search within 90^o cone angle')
        elseif n == 2 && a == 2
            title('Length prioritized search [15 10 5] within 90^o cone angle')
        elseif n == 3  && a == 2
            title('Length highly prioritized search [30 15 10] within 90^o cone angle')
        elseif n == 1 && a == 3
            title('Un-prioritized search within 30^o cone angle')
        elseif n == 2 && a == 3
            title('Length prioritized search [15 10 5] within 30^o cone angle')
        elseif n == 3  && a == 3
            title('Length highly prioritized search [30 15 10] within 30^o cone angle')
        end
        traces = traces_master{n+(a-1)*3} ;
        hold on;
        imagesc(max(data_applied(:,:,:),[],3))
        for k=1:length(traces)
            L=traces{k};
            % choose a color scheme of the three below:
            plot3(L(:,2),L(:,1),L(:,3),'-','Color',cmap_fig(round(rand*63+1),:),'LineWidth',1); % for random colors
%             plot3(L(:,2),L(:,1),L(:,3),'-','Color',cmap_fig_blues(mod(k,length(cmap_fig_blues))+1,:),'LineWidth',1); % for blue colors
%             plot3(L(:,2),L(:,1),L(:,3),'-','Color',cmap_fig_reds(mod(k,length(cmap_fig_reds))+1,:),'LineWidth',1); % for red colors
        end
        colormap(gray)
        hold off
        daspect([1 1 1])
        xlim([1 size(data_applied,2)])
        ylim([1 size(data_applied,1)])
    end
end

