addpath(genpath('SimpleTracker'))

%% Pair MB localizations with simpletracker
max_linking_distance=12;
max_gap_closing=1;
[tracks,adjacency_tracks] = simpletracker(SR_Localizations,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', true);

%% Reconstruction of tracks 
min_length=10;
n_tracks = numel(adjacency_tracks);
colors = hsv(n_tracks);
all_points = vertcat(SR_Localizations{:});
ii=1;

figure;hold on
for i_track = 1 : n_tracks  
    tmp = adjacency_tracks{i_track};
    if length(tmp)>min_length
        Tracks{ii} = all_points(tmp, :);  
        plot(Tracks{ii}(:,2), Tracks{ii}(:, 1), 'Color', colors(i_track, :))
        ii=ii+1
        %pause(0.001);
    end
end
set(gca,'YDir','reverse')
title('Microbuble tracks')
%% Reconstruction of tracked image
FinalPoints=cell2mat(Tracks');
ind=sub2ind([size_z, size_x],round(FinalPoints(:,1)),round(FinalPoints(:,2)));
Imfinal=zeros([size_z, size_x]);
for pp=1:length(ind)
    Imfinal(ind(pp))=Imfinal(ind(pp))+1;
end

figure
imagesc(Imfinal)
colormap gray
axis image,axis off
caxis([0 2])

for ii = 1:length(Tracks)
    Traj_raw = Tracks{ii};
    Tracks_post{ii} = Kalman_func(Traj_raw);
end

%% Interpolation of the tracks
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

%% Reconstruction of tracked image
FinalPoints=cell2mat(Tracks_Final');
FinalPoints(FinalPoints<1) = 1;
for ii = 1:size(FinalPoints,1)
    if FinalPoints(ii,1)>1600
        FinalPoints(ii,1) = 1600;
    end
    if FinalPoints(ii,2)>1024
        FinalPoints(ii,2) = 1024;
    end
end
ind=sub2ind([size_z, size_x],FinalPoints(:,1),FinalPoints(:,2));
Imfinal=zeros([size_z, size_x]);
for pp=1:length(ind)
    Imfinal(ind(pp))=Imfinal(ind(pp))+1;
end

figure
imagesc(Imfinal.^0.5)
colormap hot
axis image,axis off
title('650 frames, threshold 0.25 12 2 10')

%save('removed_focus_bloc1.mat','Tracks','Tracks_Final');
%caxis([0 10])
% figure
% imagesc(20*log10(Imfinal/max(Imfinal(:))))
% colormap hot
% axis equal,axis off
% title('linking distance 8,minlength 5')
%caxis([0 10])