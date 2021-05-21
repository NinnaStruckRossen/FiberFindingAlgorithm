% Copyright 2014 Anders Kyrsting
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [corr,traces,data,corr_backup] = ...
    fibtracer(data,psf_raw,data_density,psf_density, ...
    it_max,stop_factor,angle_max,step_max,l_prio)
%% input checks
input_flag = 1;
if ndims(data) == 3
    volume_mode = 1;
elseif min(size(data)) > 1
    volume_mode = 0;
else
    fprintf('Input data not 2D or 3D')
    pause(3)
    input_flag = 0;
end
if ndims(psf_raw) ~= ndims(data)
    fprintf('Dimensional mismatch of input data and PSF')
    pause(3)
    input_flag = 0;
end
if ndims(psf_raw) ~= size(psf_density,2)
    fprintf('PSF density does not match PSF dimensions')
    pause(3)
    input_flag = 0;
end
if ndims(data) ~= size(data_density,2)
    fprintf('Input data density does not match input data dimensions')
    pause(3)
    input_flag = 0;
end
if stop_factor < 0 || stop_factor > 1
    fprintf('Stop_factor outside normal interval, proceeding..')
    pause(3)
end
if (step_max < 200) && (step_max ~= -1)
    fprintf('Step_limit_factor below 200 nm, proceeding..')
    pause(3)
end
if (step_max <= 0) && (step_max ~= -1)
    fprintf('Step_limit_factor invalid value')
    pause(3)
    input_flag = 0;
end
if numel(l_prio) < 1
    fprintf('Length_priorities empty, must contain at least one value')
    pause(3)
    input_flag = 0;
end

%% Static
psf_thr = 1/exp(1)^2; % crop psf from this fraction of max
blot_thr = psf_thr ; % crop psf from this fraction of max
blot_expansion = 2; % times larger than psf

%%
if input_flag == 1
    if (exist('./psf_scan.mat') == 2)
        %% case when Scan has been performed
        disp('Loading saved PSF scan')
        load('psf_scan.mat')
    else
        % rescale psf to match data density voxel ratio
        ratio_psf_data_xy = data_density(1)/psf_density(1);
        if volume_mode == 1
            ratio_psf_data_z = data_density(3)/psf_density(3);
            [Yi,Xi,Zi] = meshgrid(1:ratio_psf_data_xy:size(psf_raw,1), ...
                1:ratio_psf_data_xy:size(psf_raw,2), ...
                1:ratio_psf_data_z:size(psf_raw,3));
            psf = interp3(psf_raw,Xi,Yi,Zi,'bicubic');
        elseif volume_mode == 0
            [Yi,Xi] = meshgrid(1:ratio_psf_data_xy:size(psf_raw,1), ...
                1:ratio_psf_data_xy:size(psf_raw,2));
            psf = interp2(psf_raw,Xi,Yi,'bicubic');
        end
        
        % put in uint8 space
        psf = psf - min(psf(:));
        psf = psf./unique(max(psf(:)));
        psf = psf.*255;
        
        % psf crop
        psf_crop_thr = max(psf(:)) * psf_thr;
        psf_crop_elements = psf > psf_crop_thr;
        psf_crop_box = regionprops(psf_crop_elements,'Boundingbox');
        psf_crop_box = psf_crop_box.BoundingBox;
        
        %cropping
        if volume_mode == 1
            psf = psf(ceil(psf_crop_box(2)):floor(psf_crop_box(5))+ceil(psf_crop_box(2)-1), ...
                ceil(psf_crop_box(1)):floor(psf_crop_box(4))+ceil(psf_crop_box(1)-1), ...
                ceil(psf_crop_box(3)):floor(psf_crop_box(6))+ceil(psf_crop_box(3)-1));
        elseif volume_mode == 0
            psf = psf(ceil(psf_crop_box(2)):floor(psf_crop_box(4))+ceil(psf_crop_box(2)-1), ...
                ceil(psf_crop_box(1)):floor(psf_crop_box(3))+ceil(psf_crop_box(1)-1));
        end
        
        %% Generation of blotting psf
        if volume_mode == 1
            [Yi,Xi,Zi] = meshgrid(...
                1:ratio_psf_data_xy/blot_expansion:size(psf_raw,1), ...
                1:ratio_psf_data_xy/blot_expansion:size(psf_raw,2), ...
                1:ratio_psf_data_z/blot_expansion:size(psf_raw,3));
            blot = interp3(psf_raw,Xi,Yi,Zi,'bicubic');
        elseif volume_mode == 0
            [Yi,Xi] = meshgrid(...
                1:ratio_psf_data_xy/blot_expansion:size(psf_raw,1), ...
                1:ratio_psf_data_xy/blot_expansion:size(psf_raw,2));
            blot = interp2(psf_raw,Xi,Yi,'bicubic');
        end
        
        % put in uint8 space
        blot = blot - min(blot(:));
        blot = blot./unique(max(blot(:)));
        blot = blot.*255;
        
        % blot crop
        blot_crop_thr = max(psf(:)) * blot_thr;
        blot_crop_elements = blot > blot_crop_thr;
        blot_crop_box = regionprops(blot_crop_elements,'Boundingbox');
        blot_crop_box = blot_crop_box.BoundingBox;
        
        % cropping
        if volume_mode == 1
            blot = blot(ceil(blot_crop_box(2)):floor(blot_crop_box(5))+ceil(blot_crop_box(2)-1), ...
                ceil(blot_crop_box(1)):floor(blot_crop_box(4))+ceil(blot_crop_box(1)-1), ...
                ceil(blot_crop_box(3)):floor(blot_crop_box(6))+ceil(blot_crop_box(3)-1));
        elseif volume_mode == 0
            blot = blot(ceil(blot_crop_box(2)):floor(blot_crop_box(4))+ceil(blot_crop_box(2)-1), ...
                ceil(blot_crop_box(1)):floor(blot_crop_box(3))+ceil(blot_crop_box(1)-1));
        end
        
        
        %% blot and psf center
        blot_size = size(blot);
        psf_size = size(psf);
        if volume_mode == 1
            [blot_cen(1),...
                blot_cen(2),...
                blot_cen(3)] = ...
                ind2sub([size(blot,1),...
                size(blot,2),...
                size(blot,3)], ...
                find(blot == max(blot(:))));
            [psf_cen(1),psf_cen(2),psf_cen(3)] = ...
                ind2sub([size(psf,1),size(psf,2),size(psf,3)], ...
                find(psf == max(psf(:))));
        elseif volume_mode == 0
            [blot_cen(1),...
                blot_cen(2)] = ...
                ind2sub([size(blot,1),...
                size(blot,2)], ...
                find(blot == max(blot(:))));
            [psf_cen(1),psf_cen(2)] = ...
                ind2sub([size(psf,1),size(psf,2)], ...
                find(psf == max(psf(:))));
        end
        % blot binarization
        blot(blot < max(blot(:)) * blot_thr ) = 0;
        blot = logical(blot);
        
        %% data range remapping to [0:255]
        data = data - poissfit(data(:));
        data = (data./max(data(:)).*255);
        
        %% pad volume to allow psf scanning to original boundary
        if volume_mode == 0
            data_temp(...
                size(data,1)+blot_size(1)*2,...
                size(data,2)+blot_size(2)*2) = 0;
            data_temp(...
                blot_size(1):size(data,1)+ ...
                blot_size(1)-1,...
                blot_size(2):size(data,2)+ ...
                blot_size(2)-1) = data;
            data = data_temp;
        elseif volume_mode == 1
            data_temp(...
                size(data,1)+blot_size(1)*2,...
                size(data,2)+blot_size(2)*2,...
                size(data,3)+blot_size(3)*2) = 0;
            data_temp(...
                blot_size(1):size(data,1)+...
                blot_size(1)-1,...
                blot_size(2):size(data,2)...
                +blot_size(2)-1,...
                blot_size(3):size(data,3)...
                +blot_size(3)-1) = data;
            data = data_temp;
        end
        clear data_temp
        
        %% execute cross-correlation of data with psf
        corr = psfcorr(data,psf,psf_cen);
        
        % saving psf scan(time consuming to run)
        save('psf_scan.mat','psf','blot','blot_cen',...
            'blot_size','psf_cen','psf_size','volume_mode',...
            'data','corr',...
            'psf_density','data_density')
    end
    % set backup
    corr_backup = corr;
    
    %% User-variables check
    % if user value not set
    if step_max == -1
        step_max = max([size(psf,1) size(psf,2)]) * data_density(1) * 2; 
    end
    if it_max == -1
        it_max = 200;
    end
    
    %% stop criterium creation
    % de-pad
    if volume_mode == 0
        corr_hist = corr(...
            blot_size(1):end-blot_size(1)-1, ...
            blot_size(2):end-blot_size(2)-1);
        corr_hist(isnan(corr_hist)) = [];
        % signal level
        corr_signal = mean(min(corr_hist));
        % black level
        corr_hist_sorted = sort(corr_hist,1) ;
        corr_hist_sorted = corr_hist_sorted(floor(size(corr_hist_sorted,1)/2):end,:) ;
        corr_black = poissfit(corr_hist_sorted(:)) ;
    elseif volume_mode == 1
        corr_hist = corr(...
            blot_size(1):end-blot_size(1)-1, ...
            blot_size(2):end-blot_size(2)-1, ...
            blot_size(3):end-blot_size(3)-1);
        corr_hist(isnan(corr_hist)) = [];
        corr_hist_minproj = min(corr_hist,[],3);
        % signal level
        corr_signal = mean(min(corr_hist_minproj));
        % black level
        corr_black = poissfit(corr_hist(:)) ;
    end
    
    % stop criterion for object trace
    corr_limit = corr_black - (corr_black - corr_signal) * stop_factor;
    
    
    %% startingpoints grid enumeration
    init_pt_spacing = 1;
    init_pt = [];
    if volume_mode == 0
        for a = psf_size(1):psf_size(1)*init_pt_spacing:...
                size(corr,1)-psf_size(1)
            for b = psf_size(2):psf_size(2)*init_pt_spacing:...
                    size(corr,2)-psf_size(2)
                % generate 1 each size step
                init_pt = [init_pt ; a b];
            end
        end
        clear a b
    elseif volume_mode == 1
        for a = psf_size(1):psf_size(1)*init_pt_spacing:...
                size(corr,1)-psf_size(1)
            for b = psf_size(2):psf_size(2)*init_pt_spacing:...
                    size(corr,2)-psf_size(2)
                for c = psf_size(3):psf_size(3)*init_pt_spacing:...
                        size(corr,3)-psf_size(3);
                    init_pt = [init_pt ; a b c];
                end
            end
        end
        clear a b c
    end
    
    init_pt_counter = size(init_pt,1);
    
    % setup for step length limit
    if volume_mode == 0
        pad = [round(step_max/data_density(1)) ...
            round(step_max/data_density(2))];
    elseif volume_mode == 1
        pad = [round(step_max/data_density(1)) ...
            round(step_max/data_density(2)) ...
            round(step_max/data_density(3))];
    end
    
    %%  master loop
    blot_val = corr_black*2;
    status_cnt = 0;
    init_pt_max = init_pt_counter;
    blot_l_thr = l_prio(1);
    l_prio_no = 1;
        
    while init_pt_counter > 0
        % status display
        i = round(100 - init_pt_counter/init_pt_max*100);
        if i > status_cnt
            clc
            fprintf(1,'PSF correlation progress:  100%% \nTracing progress:  ');
            fprintf(1,'\b%d',i);
            fprintf(1,'%%',i);
            status_cnt = i;
        end
        
        % loop counters
        init_pt_counter = init_pt_counter-1;
        
        % picking next startingpoint
        step_prev = init_pt(init_pt_counter+1,:);
        
        % length prioritization loop
        if init_pt_counter == 0
            if l_prio_no < numel(l_prio)
                l_prio_no = l_prio_no+1;
                blot_l_thr = l_prio(l_prio_no);
                init_pt_counter = init_pt_max;
            end
        end
        
        % inits
        step_allowed = 1;
        step_num = 0;
        step_direction = 1;
        
        %% main line tracing by PSF
        while (step_allowed == 1) && (step_num < it_max)
            % counter
            step_num = step_num+1;
            
            %% creating subvolume from ranges
            if volume_mode == 0
                % test if subvolume will fall outside the global boundary
                % if subvolume would have exceeded boundary, save the
                % offset xy assumption that it will not violate both lower
                % and upper boundary
                xrange = step_prev(1) - pad(1) :...
                    step_prev(1) + pad(1);
                yrange = step_prev(2) - pad(2) :...
                    step_prev(2) + pad(2);
                
                % offset test
                offset(1) = + numel(xrange(xrange < 1 ));
                offset(2) = + numel(yrange(yrange < 1 ));
                
                % removing elements breaking boundary (set to NaN)
                xrange(xrange > size(corr,1) ) = NaN;
                xrange(xrange < 1 ) = NaN;
                xrange(isnan(xrange)) = [];
                
                yrange(yrange > size(corr,2) ) = NaN;
                yrange(yrange < 1 ) = NaN;
                yrange(isnan(yrange)) = [];
                
                sub_corr = corr(xrange,yrange);
                
            elseif volume_mode == 1
                % test if subvolume will fall outside the global boundary
                % if subvolume would have exceeded boundary, save the
                % offset xy assumption that it will not violate both lower
                % and upper boundary
                xrange = step_prev(1) - pad(1) :...
                    step_prev(1) + pad(1);
                yrange = step_prev(2) - pad(2) :...
                    step_prev(2) + pad(2);
                zrange = step_prev(3) - pad(3) :...
                    step_prev(3) + pad(3);
                
                % offset test
                offset(1) = + numel(xrange(xrange < 1 ));
                offset(2) = + numel(yrange(yrange < 1 ));
                offset(3) = + numel(zrange(zrange < 1 ));
                
                % removing elements breaking boundary (set to NaN)
                xrange(xrange > size(corr,1) ) = NaN;
                xrange(xrange < 1 ) = NaN;
                xrange(isnan(xrange)) = [];
                yrange(yrange > size(corr,2) ) = NaN;
                yrange(yrange < 1 ) = NaN;
                yrange(isnan(yrange)) = [];
                zrange(zrange > size(corr,3) ) = NaN;
                zrange(zrange < 1 ) = NaN;
                zrange(isnan(zrange)) = [];
                
                sub_corr = corr(xrange,yrange,zrange);
            end
            
            %% blot by psf to avoid selecting 'step_prev' as 'step_curr'.
            % getting step_prev in sub_volume indices
            sub_pt = pad - offset;
            if step_num > 1
                %cont. on errors
                if volume_mode == 0
                    try
                        blot_vol =  sub_corr(...
                            sub_pt(1)-blot_cen(1):...
                            sub_pt(1)+(blot_size(1)...
                            -blot_cen(1)-1), ...
                            sub_pt(2)-blot_cen(2):...
                            sub_pt(2)+(blot_size(2)...
                            -blot_cen(2)-1));
                        blot_vol(blot) = blot_val;
                        sub_corr(...
                            sub_pt(1)-blot_cen(1):...
                            sub_pt(1)+(blot_size(1)...
                            -blot_cen(1)-1), ...
                            sub_pt(2)-blot_cen(2):...
                            sub_pt(2)+(blot_size(2)...
                            -blot_cen(2)-1)) = blot_vol;
                    end
                elseif volume_mode == 1
                    try
                        blot_vol =  sub_corr(...
                            sub_pt(1)-blot_cen(1):...
                            sub_pt(1)+(blot_size(1)...
                            -blot_cen(1)-1), ...
                            sub_pt(2)-blot_cen(2):...
                            sub_pt(2)+(blot_size(2)...
                            -blot_cen(2)-1), ...
                            sub_pt(3)-blot_cen(3):...
                            sub_pt(3)+(blot_size(3)...
                            -blot_cen(3)-1));
                        blot_vol(blot) = blot_val;
                        sub_corr(...
                            sub_pt(1)-blot_cen(1):...
                            sub_pt(1)+(blot_size(1)...
                            -blot_cen(1)-1), ...
                            sub_pt(2)-blot_cen(2):...
                            sub_pt(2)+(blot_size(2)...
                            -blot_cen(2)-1), ...
                            sub_pt(3)-blot_cen(3):...
                            sub_pt(3)+(blot_size(3)...
                            -blot_cen(3)-1)) = blot_vol;
                    end
                end
            end
            %% Angle change cost enforcement on correlation map
            %validity only when trace is present ( min 2 points)
            if step_allowed == 1
                % getting vector
                if step_num > 2
                    % get last step
                    vct = trace_found(end,:) - trace_found(end-1,:);
                    
                    %% running angle map, if valid step_prev
                    if (unique(min(sub_pt >= 1)) == 1) && ...
                            (unique(min(sub_pt <= size(sub_corr))) == 1)
                        %% Determining best step
                        step_allowed = 0;
                        best_step_found = 0;
                        step_loop = 1;
                        sub_size = size(sub_corr);
                        while (best_step_found == 0) && (step_loop < length(sub_corr(:)))
                            step_loop = step_loop+1;
                            [best_step_val,best_step_id] = min(sub_corr(:));
                            best_step_val = best_step_val(1);
                            best_step_id = best_step_id(1);
                            if best_step_val > corr_limit
                                break
                            end
                            if volume_mode == 0
                                [best_step_coor(1),best_step_coor(2)] = ...
                                    ind2sub(sub_size,best_step_id);
                                theta =  acosd( ( ...
                                    (best_step_coor(1)-sub_pt(1))*vct(1) + ...
                                    (best_step_coor(2)-sub_pt(2))*vct(2)) / ...
                                    (norm(best_step_coor - ...
                                    [sub_pt(1) sub_pt(2)] )* ...
                                    norm([vct(1) vct(2)] )) );
                                % avoid outside of cone
                                if theta < angle_max
                                    %test euclidean distance
                                    eucl_dist = sqrt(((best_step_coor(1)-sub_pt(1))*data_density(1))^2 + ...
                                        ((best_step_coor(2)-sub_pt(2))*data_density(2))^2 );
                                    if eucl_dist < step_max
                                        step_allowed = 1;
                                        best_step_found = 1;
                                    else
                                        sub_corr(best_step_coor(1),best_step_coor(2)) = NaN;
                                    end
                                else
                                    sub_corr(best_step_coor(1),best_step_coor(2)) = NaN;
                                end
                            elseif volume_mode == 1
                                [best_step_coor(1),best_step_coor(2),best_step_coor(3)] = ...
                                    ind2sub(sub_size,best_step_id);
                                theta =  acosd( ( ...
                                    (best_step_coor(1)-sub_pt(1))*vct(1) + ...
                                    (best_step_coor(2)-sub_pt(2))*vct(2) + ...
                                    (best_step_coor(3)-sub_pt(3))*vct(3)) / ...
                                    (norm(best_step_coor - ...
                                    [sub_pt(1) sub_pt(2) sub_pt(3)] )* ...
                                    norm([vct(1) vct(2) vct(3)] )) );
                                % avoid outside of cone
                                if theta < angle_max
                                    %test euclidean distance
                                    eucl_dist = sqrt(((best_step_coor(1)-sub_pt(1))*data_density(1))^2 + ...
                                        ((best_step_coor(2)-sub_pt(2))*data_density(2))^2 + ...
                                        ((best_step_coor(3)-sub_pt(3))*data_density(2)*data_density(2)/data_density(3) )^2 );
                                    if eucl_dist < step_max
                                        step_allowed = 1;
                                        best_step_found = 1;
                                    else
                                        sub_corr(best_step_coor(1),best_step_coor(2),best_step_coor(3)) = NaN;
                                    end
                                else
                                    sub_corr(best_step_coor(1),best_step_coor(2),best_step_coor(3)) = NaN;
                                end
                            end
                        end
                    else
                        % step_prev is beyond test data (near boundary),
                        % cancel step
                        best_step_val = corr_limit*2;
                    end
                else
                    % unconstrained first step
                    % blot by psf to avoid reselecting
                    % 'step_prev' as 'step_curr'.
                    if volume_mode == 0
                        if step_num == 2
                            try
                                sub_corr(...
                                    sub_pt(1)-psf_cen(1):sub_pt(1)+(psf_size(1)-psf_cen(1)-1), ...
                                    sub_pt(2)-psf_cen(2):sub_pt(2)+(psf_size(2)-psf_cen(2)-1)) = NaN;
                            end
                        end
                        best_step_val = unique(min(sub_corr(:)));
                        [best_step_coor(1),best_step_coor(2)] = ...
                            ind2sub(size(sub_corr),find(sub_corr == best_step_val,1));
                    elseif volume_mode == 1
                        if step_num == 2
                            try
                                sub_corr(...
                                    sub_pt(1)-psf_cen(1):sub_pt(1)+(psf_size(1)-psf_cen(1)-1), ...
                                    sub_pt(2)-psf_cen(2):sub_pt(2)+(psf_size(2)-psf_cen(2)-1), ...
                                    sub_pt(3)-psf_cen(3):sub_pt(3)+(psf_size(3)-psf_cen(3)-1)) = NaN;
                            end
                        end
                        best_step_val = unique(min(sub_corr(:)));
                        [best_step_coor(1),best_step_coor(2),best_step_coor(3)] = ...
                            ind2sub(size(sub_corr),find(sub_corr == best_step_val,1));
                    end
                    
                end
            end
            
            
            %% verify that best step does not result in NaN
            if step_allowed == 1
                if isnan(best_step_val) == 1
                    step_allowed = 0;
                end
            end
            
            %% determining best step, and checking for acceptable correlation
            if step_allowed == 1
                % checking against stop limit
                if best_step_val >= corr_limit
                    step_allowed = 0;
                else
                    % if step is good, pick coords for it
                    % correct for offset
                    if volume_mode == 0
                        best_step_coor(1) = best_step_coor(1) + min(xrange);
                        best_step_coor(2) = best_step_coor(2) + min(yrange);
                    elseif volume_mode == 1
                        best_step_coor(1) = best_step_coor(1) + min(xrange);
                        best_step_coor(2) = best_step_coor(2) + min(yrange);
                        best_step_coor(3) = best_step_coor(3) + min(zrange);
                    end
                end
            end
            
            %% testing for trace procession
            if step_allowed == 1
                if best_step_coor == step_prev
                    step_allowed = 0;
                end
            end
            
            %% saving step into trace
            if step_allowed == 1
                if step_num == 1
                    trace_found = best_step_coor;
                else
                    trace_found(size(trace_found,1)+1,:) = best_step_coor;
                end
            end
            
            %% blotting trace from search map, only blot if trace is at least n steps long.
            %  This avoids odd detection from ruining proper traces.
            if step_allowed == 1
                % when trace is accepted for blotting, blot previous
                % steps in trace
                if size(trace_found,1) == blot_l_thr
                    blot_coords = trace_found;
                elseif size(trace_found,1) > blot_l_thr
                    blot_coords = trace_found(end-1:end,:);
                end
                if size(trace_found,1) >= blot_l_thr
                    for m = 1:size(blot_coords,1)-1
                        if volume_mode == 0
                            % creating bilinear points between current step
                            % and step_prev
                            blot_numel = max(abs( blot_coords(m,:) - blot_coords(m+1,:) ));
                            
                            % trace coordinates in simple line
                            blot_x = floor(linspace(blot_coords(m+1,1),blot_coords(m,1),blot_numel));
                            blot_y = floor(linspace(blot_coords(m+1,2),blot_coords(m,2),blot_numel));
                            
                            % blotting trace
                            for k = 1:length(blot_x)
                                try
                                    % creating blotting block
                                    blot_vol = corr(...
                                        blot_x(k)-blot_cen(1):blot_x(k)+(blot_size(1)-blot_cen(1)-1), ...
                                        blot_y(k)-blot_cen(2):blot_y(k)+(blot_size(2)-blot_cen(2)-1));
                                    
                                    % blotting by psf spot
                                    blot_vol(blot) = blot_val;
                                    corr(...
                                        blot_x(k)-blot_cen(1):blot_x(k)+(blot_size(1)-blot_cen(1)-1), ...
                                        blot_y(k)-blot_cen(2):blot_y(k)+(blot_size(2)-blot_cen(2)-1)) = blot_vol;
                                end
                            end
                        elseif volume_mode == 1
                            % creating bilinear points between current step
                            % and step_prev
                            blot_numel = max(abs( blot_coords(m,:) - blot_coords(m+1,:) ));
                            
                            % trace coordinates in simple line
                            blot_x = floor(linspace(blot_coords(m+1,1),blot_coords(m,1),blot_numel));
                            blot_y = floor(linspace(blot_coords(m+1,2),blot_coords(m,2),blot_numel));
                            blot_z = floor(linspace(blot_coords(m+1,3),blot_coords(m,3),blot_numel));
                            
                            % blotting trace
                            for k = 1:length(blot_x)
                                try
                                    % creating blotting block
                                    blot_vol = corr(...
                                        blot_x(k)-blot_cen(1):blot_x(k)+(blot_size(1)-blot_cen(1)-1), ...
                                        blot_y(k)-blot_cen(2):blot_y(k)+(blot_size(2)-blot_cen(2)-1), ...
                                        blot_z(k)-blot_cen(3):blot_z(k)+(blot_size(3)-blot_cen(3)-1));
                                    
                                    % blotting by psf spot
                                    blot_vol(blot) = blot_val;
                                    corr(...
                                        blot_x(k)-blot_cen(1):blot_x(k)+(blot_size(1)-blot_cen(1)-1), ...
                                        blot_y(k)-blot_cen(2):blot_y(k)+(blot_size(2)-blot_cen(2)-1), ...
                                        blot_z(k)-blot_cen(3):blot_z(k)+(blot_size(3)-blot_cen(3)-1)) = blot_vol;
                                end
                            end
                        end
                    end
                end
            end
            
            %% saving states for next iteration
            if step_allowed == 1
                step_prev = best_step_coor;
            end
            %% When trace ends in first direction, start over from the other end
            if (step_allowed == 0) && (step_direction == 1)
                if step_num > 2
                    % switch to reverse tracking from startpoint
                    step_direction = 2;
                    
                    % append in right order to trace
                    trace_found = flipud(trace_found);
                    
                    % attempting restart from startpoint
                    step_prev = trace_found(end,:);
                    
                    % forcing loop to go on
                    step_allowed = 1;
                end
            end
        end
        
        %% saving trace
        % check that startpoint proceeded to join
        if exist('trace_found') == 1 && size(trace_found,1) > blot_l_thr-1
            % correct for space-mismatch
            trace_found = trace_found-1;
            if exist('traces') == 0
                traces{1,1} = trace_found;
            else
                traces{1,size(traces,2)+1} = trace_found;
            end
        end
        %% clear vars
        clear trace_found
    end
end

if exist('traces') == 0
    fprintf('No traces found using user defined parameters.')
    pause(3)
    traces = NaN;
end

%% correlation of input volume 2 or 3D, with a same dimensionality psf
function [residuals] = psfcorr(volume_in,psf,offset)
vol_size = size(volume_in);
psf_size = size(psf);
residuals = NaN(vol_size);
% 2D data
if length(vol_size) == 2
    % loop setup
    x_range = 1:vol_size(1)-psf_size(1)+1;
    y_range = 1:vol_size(2)-psf_size(2)+1;
    for x = x_range
        % status display
        clc
        fprintf(1,'PSF correlation progress:  ');
        i = round(x/max(x_range)*100);
        fprintf(1,'\b%d',i);

        for y = y_range
            residual_sq = (volume_in( x:x+psf_size(1)-1, ...
                y:y+psf_size(2)-1) -psf).^2;
            % save this correlation with PSF, coord in offsetter
            residuals(x+offset(1)-1,y+offset(2)-1) = sum(residual_sq(:));
        end %end y
    end %end x
    % 3D data
elseif length(vol_size) == 3
    % loop setup
    x_range = 1:vol_size(1)-psf_size(1)+1;
    y_range = 1:vol_size(2)-psf_size(2)+1;
    z_range = 1:vol_size(3)-psf_size(3)+1;
    for x = x_range
        % status display
        clc
        fprintf(1,'PSF correlation progress:  ');
        i = round(x/max(x_range)*100);
        fprintf(1,'\b%d',i);
        fprintf(1,'%%',i);

        for y = y_range
            for z = z_range
                residual_sq = (volume_in( x:x+psf_size(1)-1, ...
                                          y:y+psf_size(2)-1, ...
                                          z:z+psf_size(3)-1)-psf).^2;
                % save this correlation with PSF, coord in offsetter
                residuals(x+offset(1)-1,y+offset(2)-1,z+offset(3)-1) = sum(residual_sq(:));
            end %end z
        end %end y
    end %end x
end