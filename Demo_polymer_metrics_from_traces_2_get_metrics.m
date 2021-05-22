%% FIRST RUN THE FIBER FINDING ALGORITHM ON A 3D IMAGES

% (Then either continue running this or save the workspace and load it
%  before running this.)

% This demo calculates the following metrics of interest for the found
% fibers (bio-polymers):
% 1) End-to-end length and contour length for each individual polymer and
%    Accumulated contour length for all the polymers in the sample
% 2) Persistence length for the sample
% 3) 3D Mesh size using all pixels in the voids for the sample
% 4) 3D Mesh size using only the pixels with max distance to nearest fiber in each void for the sample

    
%% 1. Calculating end-to-end length (and contour lengths) of the individual identified fibers

L_end2end = NaN(1,length(traces));
L_contour = zeros(1,length(traces));
for k=1:length(traces)
    L=traces{k};
    L_micron = L.*(data_density*10^(-3));
    L_end2end(1,k) = norm(L_micron(1,:)-L_micron(end,:)); % sqrt(sum((L_micron(1,:)-L_micron(end,:)).^2));
    for i=1:(length(L_micron)-1)
        L_contour(1,k) = L_contour(1,k) + norm(L_micron(i,:)-L_micron(i+1,:));
    end
end

% Calculating mean and standard deviation of contour length 
% by fitting the distribution of lengths to a lognormal distribution
% (other distributions can be used if appropriate)

% prepare data for fit
[x_contour, yData] = prepareCurveData(L_contour,L_end2end.^2);
% x_contour = 0:.1:ceil(max(xData,[],1));
% calculating mean from data 
Contour_length_mean_data = log(mean(exp(x_contour)));

% fit lognormal distribution to contour lengths
pd_contour_lognormal = fitdist(x_contour,'lognormal');
pdf_contour_lognormal = pdf(pd_contour_lognormal,x_contour);

% Calculating mean and standard deviation from distribution parameters
Contour_length_mean_lognormal = exp(pd_contour_lognormal.mu+(pd_contour_lognormal.sigma^2)/2);
Contour_length_std_lognormal = sqrt(exp(2*pd_contour_lognormal.mu + pd_contour_lognormal.sigma^2)...
                                    *(exp(pd_contour_lognormal.sigma^2)-1));

% Calculating accumulated contour length 
Accumulated_Contour_length = sum(L_contour);

%% 2. Calculating persistence length by fitting L_endtoend^2 vs. L_contour of all identified fibers

% Set function for fit
ft = fittype( '2*e_p*((x))-2*e_p^2*(1-exp(-((x))/e_p))', 'independent', 'x', 'dependent', 'y' );

% Set fit options for either weighted non-weighted or fit:

% FIT 1: weighing longer fibers more
opts1 = fitoptions( ft );
opts1.Display = 'Off';
opts1.Lower = [-Inf];
opts1.StartPoint = [0.35];
opts1.Upper = [Inf];
opts1.Weights = x_contour./max(x_contour);
% fit model to data.
[fitresult_weighted, gof] = fit( x_contour, yData, ft, opts1);

% FIT 2: equal weight to all fibers
opts2 = fitoptions( ft );
opts2.Display = 'Off';
opts2.Lower = [-Inf];
opts2.StartPoint = [0.35];
opts2.Upper = [Inf];
opts2.Weights = ones(size(x_contour));
% fit model to data.
[fitresult_nonweighted, gof] = fit( x_contour, yData, ft, opts2);

% Calculating persistence length
Persistence_length_weighted = fitresult_weighted.e_p;
Persistence_length_nonweighted = fitresult_nonweighted.e_p;

%% 3. Calculating 3D mesh size

% rescale traces to isotropic spacing using the ratio between xy and z
ratio_xy_z = data_density(1,3)/data_density(1,1);
% create new (empty) data matrix of the correct isotropic volume
volume_data = zeros(size(data,1),size(data,2),round(size(data,3)*ratio_xy_z));
% draw in traces using  contour coordinates and Bresenham's line algorithm
% to connect them:
for k=1:length(traces)
    L=traces{k};
    for i=1:(length(L)-1)
        volume_data(L(i,1),L(i,2),round(L(i,3)*ratio_xy_z)) = 1;
        [X_temp, Y_temp, Z_temp] = bresenham_line3d([L(i,1),L(i,2),round(L(i,3)*ratio_xy_z)],[L(i+1,1),L(i+1,2),round(L(i+1,3)*ratio_xy_z)]);
        for j=1:length(X_temp)
            volume_data(X_temp(j),Y_temp(j),Z_temp(j)) = 1;
        end
    end
        volume_data(L(end,:)) = 1;
end

% calculate 3D distance map
distance_map = bwdist(volume_data)*data_density(1,3)*10^(-3);

% fit only non-zero distances
mesh_Data = double(distance_map(distance_map > 0));

% FIT TO A GAMMA DISTRIBUTION:
% PDF(x) = 1/(Gamma(k)*theta^k) * x^(k-1) * exp (-x/theta)
pd_gamma = fitdist(mesh_Data,'gamma');
x_mesh = 0:.1:ceil(max(mesh_Data,[],1));
pdf_gamma = pdf(pd_gamma,x_mesh);

% Calculate mean 3D mesh size from distribution parameters
Mesh_size_mean = pd_gamma.a*pd_gamma.b;
% Calculate stats and goodness-of-fit chi^2 test
[h_gamma,p_gamma,stats_gamma] = chi2gof(mesh_Data,'CDF',pd_gamma);


%% 4. Calculating Max 3D mesh size based on the maximum distance 
% to the nearest fiber within each void:
distance_map_localmax = imregionalmax(distance_map);
MES_radii = distance_map(distance_map_localmax);

max_traced_fiber = round(sqrt((sqrt(size(data_applied,1)^2+size(data_applied,2)^2)*data_density(1,1)*10^(-3))^2 + (size(data_applied,3)*data_density(1,3)*10^(-3))^2));
[index] = find(MES_radii < max_traced_fiber);
MES_radii_cutoff = double(MES_radii(index));

[xData_MES, yData_MES] = prepareCurveData(MES_radii_cutoff(:),MES_radii_cutoff(:));
x_MES = 0:.1:ceil(max(xData_MES,[],1));
pd_MES = fitdist(xData_MES,'gamma');
pdf_MES = pdf(pd_MES,x_MES);

% Calculate mean and standard deviation 3D mesh size from distribution parameters
Max_mesh_size_mean = pd_MES.a*pd_MES.b;
Max_mesh_size_std = sqrt(pd_MES.a*pd_MES.b^2);


%% Figures

figure(2)
subplot(2,2,1)
h = histogram(L_contour,20,'Normalization','pdf')
hold on
[x_contour_sort,index_sort] = sort(x_contour);
plot(x_contour_sort,pdf_contour_lognormal(index_sort),'LineWidth',2,'color',[0,0,0.4])
hold off
legend('Fiber contour lengths','Lognormal fitted distribution','location','northeast')
xlabel( 'Contour length [micron]');
ylabel( 'PDF normalized fiber Counts');
xlim([0 max_traced_fiber])
title('Histogram of contour lengths')

%%%%%%%%%%%%%

figure(2)
subplot(2,2,3)
plot(L_contour,L_end2end.^2,'.b')
L_range = [0:1:max(L_contour,[],2)];
hold on
plot(L_range,L_range.^2,'-','Color',[0.5,0.5,0.5])
plot(L_range,2*fitresult_weighted.e_p*(L_range)-2*fitresult_weighted.e_p^2*(1-exp(-(L_range)/fitresult_weighted.e_p)),':','Color',[0,0.5,0.5])
plot(L_range,2*fitresult_nonweighted.e_p*(L_range)-2*fitresult_nonweighted.e_p^2*(1-exp(-(L_range)/fitresult_nonweighted.e_p)),'--','Color',[0,0,1])
hold off
legend('Identified fibers','Completely straight fibers','Persistence length (weighted fit)','Persistence length (non-weighted fit)','location','northwest')
xlabel( 'Contour length [micron]');
ylabel( 'Squared end-to-end length [micron^2]');
xlim([0 max_traced_fiber*2/3])
ylim([0 (max_traced_fiber*2/3)^2])
title('Calculation of persistence length, P')

%%%%%%%%%%%%%

figure(2)
subplot(2,2,2)
h = histogram(mesh_Data,40,'Normalization','pdf')
hold on
plot(x_mesh,pdf_gamma,'LineWidth',2,'color',[0,0,0.4])
hold off
legend('3D Mesh size','Gamma fitted distribution','location','northeast')
xlabel('Ddistance to nearest fiber [micron]')
ylabel('PDF normalized counts')
xlim([0 max_traced_fiber])
ylim([0 0.25])
title('Histogram of 3D mesh size')

%%%%%%%%%%%%%

figure(2)
subplot(2,2,4)
h_MES = histogram(xData_MES,30,'Normalization','pdf');
hold on
plot(x_MES,pdf_MES,'LineWidth',2,'color',[0,0,0.4])
hold off
legend(' Max 3D Mesh size ',' Gamma fitted distribution','location','northeast')
xlabel('Max distance to nearest fiber [micron]')
ylabel('PDF normalized counts')
xlim([0 max_traced_fiber])
ylim([0 0.25])
title('Histogram of Max 3D mesh size')