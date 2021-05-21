%% Input measured actual point-spread-function (psf) or theoretical 3D gaussian approximation

% Measured actual psf:
FileTif= './Data_ourPSF.tif';
% change "Data_ourPSF" to the measured PSF image of your imaging system and
% objective above if you are analysing your own small dataset.
myImInfo = imfinfo(FileTif,'tiff');
z_range = 1:length(myImInfo);
x_el_needed = myImInfo(1).Height ;
y_el_needed = myImInfo(1).Width ;
psf_actual = zeros(x_el_needed,y_el_needed,numel(z_range),'double') ;
for i=z_range
    psf_actual(:,:,i)=double(imread(FileTif,'Index',i));
end

% Theoretical 3D gaussian approximation of psf
point_emitter = zeros(46,71,31);
point_emitter(23,36,16) = [255];
psf_theoretical = imgaussfilt3(point_emitter,5);
psf_theoretical = double(round(psf_theoretical*255*5));

for i=1:size(psf_theoretical,3)
    imwrite(uint8(psf_theoretical(:,:,i)),'Data_theoreticalPSF.tif','WriteMode','append')
end

figure(2)
for i=1:6
subplot(6,3,i*3-2)
imagesc(psf_actual(:,:,1+(5*(i-1))))
caxis([0 200])
axis equal
title(strcat('actual PSF z=',num2str(1+(5*(i-1)))))

subplot(6,3,i*3-1)
imagesc(psf_theoretical(:,:,1+(5*(i-1))))
caxis([0 200])
axis equal
title(strcat('theoretical PSF z=',num2str(1+(5*(i-1)))))

subplot(6,3,i*3)
imagesc(point_emitter(:,:,1+(5*(i-1))))
caxis([0 200])
axis equal
title(strcat('point emitter z=',num2str(1+(5*(i-1)))))
end