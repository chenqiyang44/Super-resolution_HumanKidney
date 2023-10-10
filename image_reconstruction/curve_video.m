load('param.mat');
%IQ_corr = BmodeData(1:1300,:,:);
SoundSpeed = 1540;%Resource.Parameters.speedOfSound;
dataDepth = size(BmodeData,1);
fs = Receive(1).samplesPerWave * Trans.frequency;    % In MHz
lat = Trans.ElementPos(:,1)'*(SoundSpeed*1e3)/(Trans.frequency*1e6); % Center 128 channel
axial = ((1:(dataDepth))*SoundSpeed/(fs*1e6)/2+Receive(1).startDepth/Trans.frequency*SoundSpeed/1e6)*1e3; %mm % US Pulse-echo case
% polar coordinate
radius=axial(1:1600);
radiusOfCurvature = Trans.radiusMm;
radius = radius + radiusOfCurvature;
theta = linspace(Trans.ElementPos(1,4),Trans.ElementPos(end,4),size(Trans.ElementPos,1));
%theta = Trans.ElementPos(:,4);
[THETA,R] = meshgrid(theta,radius);
[X,Y] = pol2cart(THETA,R);
writerObj = VideoWriter('Trial1.avi');
writerObj.FrameRate = 500;
% set the seconds per image
% open the video writer
open(writerObj);
for i=1:3000
    a=abs(BmodeData(1:1600,:,i));
    figure(1);surf(Y,X-radiusOfCurvature,20*log10(a/max(a(:))),'edgecolor','none'),title(['Bmode Frame' num2str(i)]);
    %figure(1);surf(Y,X-radiusOfCurvature,a,'edgecolor','none'),title([num2str(i)]);
    view(0,90);
    %caxis([0.05 0.2])
    caxis([-40 0]);
    xlabel('lateral [mm]');
    ylabel('axial [mm]');
    colormap gray;
    %colorbar;
    axis image;
    set(gca,'Ydir','reverse');
    backColor=[0 0 0];
    set(gca,'color',backColor);
    %axis off
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end
close(writerObj);

