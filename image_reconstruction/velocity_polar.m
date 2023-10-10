%%parameter initiation
load('param');load('MyColormap');
f0 = Trans.frequency;
fs = 4*f0;
SoundSpeed=1540;
dataDepth = 1600;
axial = ((1:(dataDepth))*SoundSpeed/(fs*1e6)/2+Receive(1).startDepth/Trans.frequency*SoundSpeed/1e6)*1e3; %mm % US Pulse-echo case
%lat = Trans.ElementPos(:,1)'*(SoundSpeed*1e3)/(Trans.frequency*1e6); % Center 128 channel
theta = Trans.ElementPos(:,4);

%index_x_SR=linspace(lat(1),lat(end),5*length(lat));  %psf radius: 0.4312mm, lat: 0.7451mm,22 indices in lateral, 10 indices in axial
index_y_SR=linspace(axial(1),axial(end),length(axial));
theta_SR=linspace(theta(1),theta(end),8*length(theta));
z_slope = mean(diff(index_y_SR));
z_bias = index_y_SR(1);
x_slope = mean(diff(theta_SR));
x_bias = theta_SR(1);
radiusOfCurvature = Trans.radiusMm;
%% Interpolation of the tracks
VelocityMap_final = zeros(1600,1024);
V_index = zeros(1600,1024);
pixel_value = 1540/3.125/1e6/8*1e3; %mm
n_tracks = numel(Tracks);
nn = 1;
for i_track=1:n_tracks
    z_track=Tracks_post{i_track}(:,1);
    x_track=Tracks_post{i_track}(:,2);
    velocity = zeros(length(z_track)-1,1);
    for j = 1:length(z_track)-1
        %calculation velocity
%         line1=z_track(j+1) - z_track(j);
%         line2=x_track(j+1) - x_track(j);
%         flow_direction=z_track(j+1)-z_track(j);
%         dis=sqrt(line1^2+line2^2) * pixel_value;
%         velocity(j)=dis/(1/500)*sign(-flow_direction);
        %velocity(j)=dis/(1/40)*sign(-flow_direction);
        
        %Velocity in polar coordinate
        y_pre = (z_track(j)-1)*z_slope + z_bias;
        y_cur = (z_track(j+1)-1)*z_slope + z_bias;
        x_pre = (x_track(j)-1)*x_slope + x_bias;
        x_cur = (x_track(j+1)-1)*x_slope + x_bias;
        
        line1=(y_pre+radiusOfCurvature)*sin(x_pre-x_cur);
        line2=y_cur+radiusOfCurvature-(y_pre+radiusOfCurvature)*cos(x_pre-x_cur);
        flag_direction=(y_pre+radiusOfCurvature)*cos(x_pre)-(y_cur+radiusOfCurvature)*cos(x_cur);
        dis=sqrt(line1^2+line2^2);
        velocity(j)=dis/(1/500)*sign(flag_direction);      

    end
    z_track_interp = interp1(1:length(z_track),z_track,1:0.1:length(z_track));
    x_track_interp = interp1(1:length(x_track),x_track,1:0.1:length(x_track));
    x_final=round(x_track_interp(1:end));
    z_final=round(z_track_interp(1:end));    
    [~,ixu]=unique(x_final);[~,izu]=unique(z_final);
    ifin=union(izu,ixu);
    for j = 1:length(ifin)
        if x_final(ifin(j))<1 
            x_final(ifin(j)) = 1;
        end
        if x_final(ifin(j))>1024 
            x_final(ifin(j)) = 1024;
        end
        if z_final(ifin(j))> 1600 
            z_final(ifin(j)) = 1600;
        end
        if z_final(ifin(j))< 1 
            z_final(ifin(j)) = 1;
        end
        VelocityMap_final(z_final(ifin(j)),x_final(ifin(j))) = velocity(floor(abs((ifin(j)-2))/10)+1);
        %VelocityMap(z_final(ifin(j)),x_final(ifin(j))) = VelocityMap(z_final(ifin(j)),x_final(ifin(j))) + velocity(floor(abs((ifin(j)-2))/10)+1);
        %V_index(z_final(ifin(j)),x_final(ifin(j))) = V_index(z_final(ifin(j)),x_final(ifin(j))) + 1;
    end
end
%VelocityMap_final = VelocityMap ./ V_index;
%% Reconstruction of tracked image
index_y_SR1 = index_y_SR + radiusOfCurvature;
[THETA,R] = meshgrid(theta_SR',index_y_SR1);
[X,Y] = pol2cart(THETA,R);

VelocityMap_final(isnan(VelocityMap_final))=0;
figure;p=pcolor(Y,X-radiusOfCurvature,VelocityMap_final);
set(p,'EdgeColor','none');
view(0,90);colorbar;axis image;colormap(cmap);caxis([-350 350])
%xlabel('lateral [mm]');ylabel('axial [mm]');
xlim([-50 50]);ylim([0 90]);
set(gca,'Ydir','reverse');backColor=[0 0 0];set(gca,'color',backColor);
%title('trial 4 thresh 0.20 2.1 0.7 1600 frames 0.80 size 30 pcolor');

