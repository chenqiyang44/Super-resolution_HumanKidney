clear
size_z = 1600;
size_x = 1024;
min_length = 5;
ii = 1;
for i = 1:3
    %ii = 1;
    filename = ['2_500ensemble_corr_thresh01_bloc' num2str(i) '.mat'];
    load(filename);
    n_tracks = numel(adjacency_tracks);
    all_points = vertcat(SR_Localizations{:});
    size(all_points,1)
    for i_track = 1 : n_tracks  
        tmp = adjacency_tracks{i_track};
        if length(tmp)>min_length
            Tracks{ii} = all_points(tmp, :);  
            ii=ii+1;
        end
    end
end

% Apply Kalman Filter
for nn = 1:length(Tracks)
    Traj_raw = Tracks{nn};
    Tracks_post{nn} = Kalman_func(Traj_raw);
end
% Interpolation of the tracks
n_tracks = numel(Tracks);
for i_track=1:n_tracks
    z_track=Tracks_post{i_track}(:,1);
    x_track=Tracks_post{i_track}(:,2);    
    z_track_interp = interp1(1:length(z_track),z_track,1:0.1:length(z_track));
    x_track_interp = interp1(1:length(x_track),x_track,1:0.1:length(x_track));
    x_final=round(x_track_interp(1:end));
    z_final=round(z_track_interp(1:end));    
    [~,ixu]=unique(x_final);[~,izu]=unique(z_final);
    ifin=union(izu,ixu);
    Tracks_Final{i_track}(:,1)=z_final(ifin);
    Tracks_Final{i_track}(:,2)=x_final(ifin);
end

% Reconstruction of tracked image
FinalPoints=cell2mat(Tracks_Final');
FinalPoints(FinalPoints<1) = 1;
for nn = 1:size(FinalPoints,1)
    if FinalPoints(nn,1)>1600
        FinalPoints(nn,1) = 1600;
    end
    if FinalPoints(nn,2)>1024
        FinalPoints(nn,2) = 1024;
    end
end
ind=sub2ind([size_z, size_x],FinalPoints(:,1),FinalPoints(:,2));
Imfinal=zeros([size_z, size_x]);
for pp=1:length(ind)
    Imfinal(ind(pp))=Imfinal(ind(pp))+1;
end

figure
imagesc(Imfinal.^0.5)
colormap hot,caxis([0 4])
axis image,axis off
title('tracking 2 1 15')


index_y_SR1 = index_y_SR + radiusOfCurvature;
[THETA,R] = meshgrid(theta_SR',index_y_SR1);
[X,Y] = pol2cart(THETA,R);
figure;p=pcolor(Y,X-radiusOfCurvature,Imfinal.^0.5);
set(p,'EdgeColor','none');
view(0,90);colorbar;axis image;colormap(hot);caxis([0 ])
%xlabel('lateral [mm]');ylabel('axial [mm]');
xlim([-50 50]);ylim([0 90]);
set(gca,'Ydir','reverse');backColor=[0 0 0];set(gca,'color',backColor);