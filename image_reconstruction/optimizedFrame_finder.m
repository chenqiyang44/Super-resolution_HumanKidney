%% parameter initiation
clear
load('param');load('MyColormap');
f0 = Trans.frequency;
fs = 4*f0;
SoundSpeed=1540;
dataDepth = 1600;
axial = ((1:(dataDepth))*SoundSpeed/(fs*1e6)/2+Receive(1).startDepth/Trans.frequency*SoundSpeed/1e6)*1e3; %mm % US Pulse-echo case
theta = Trans.ElementPos(:,4);

index_y_SR=linspace(axial(1),axial(end),length(axial));
theta_SR=linspace(theta(1),theta(end),8*length(theta));
z_slope = mean(diff(index_y_SR));
z_bias = index_y_SR(1);
x_slope = mean(diff(theta_SR));
x_bias = theta_SR(1);
radiusOfCurvature = Trans.radiusMm;

index_y_SR1 = index_y_SR + radiusOfCurvature;
[THETA,R] = meshgrid(theta_SR',index_y_SR1);
[X,Y] = pol2cart(THETA,R);

size_x = 1024;
size_z = 1600;

min_length=10;

%% loop to find the optimized number of frames
for no_blocs = 1:11
    ii = 1;
    for blocs = 1:no_blocs
        filename_load = ['Trial4_017_link2_bloc' num2str(blocs) '.mat'];
        load(filename_load);
        n_tracks = numel(adjacency_tracks);
        all_points = vertcat(SR_Localizations_temp{:});

        for i_track = 1 : n_tracks  
            tmp = adjacency_tracks{i_track};
            if length(tmp)>min_length
                Tracks{ii} = all_points(tmp, :);  
                ii=ii+1;
            end
        end
    end
    
    % Apply Kalman Filter
    for ii = 1:length(Tracks)
        Traj_raw = Tracks{ii};
        Tracks_post{ii} = Kalman_func(Traj_raw);
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
    FinalPoints(FinalPoints<1)=1;
    for n = 1:size(FinalPoints,1)
        if FinalPoints(n,2)>1024
            FinalPoints(n,2) = 1024;
        end
    end     
    ind=sub2ind([size_z, size_x],FinalPoints(:,1),FinalPoints(:,2));
    Imfinal=zeros([size_z, size_x]);
    for pp=1:length(ind)
        Imfinal(ind(pp))=Imfinal(ind(pp))+1;
    end

    figure;p=pcolor(Y,X-radiusOfCurvature,Imfinal.^0.5);
    set(p,'EdgeColor','none');
    view(0,90);colorbar;axis image;colormap(hot);%axis off
    %xlabel('lateral [mm]');ylabel('axial [mm]');
    set(gca,'Ydir','reverse');backColor=[0 0 0];set(gca,'color',backColor);
    xlim([-50 50]);ylim([0 90]);%caxis([0 5])
    title([num2str(300+(no_blocs-1)*100) 'frames']); 
    set(gcf, 'WindowState', 'maximized')
    fig_name = [num2str(300+(no_blocs-1)*100) '.tif'];
    saveas(gcf,fig_name)
    no_blocs
end