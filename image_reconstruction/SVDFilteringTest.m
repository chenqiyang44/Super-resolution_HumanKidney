% clear;
% close all;
% load('tmpCompen_dbl');
%load('tmpEst_dbl');


%%
%cut100 = compenAll(:,:,1:100);
%cut300 = compenAll(:,:,1:300);
%cut300_2x = compenAll(:,:,1:2:600);
%cut500 = compenRawAll(:,:,1:500);
%cut500_2x = compenRawAll(:,:,1:2:1000);
%cut1000 = compenRawAll(279:end,:,1:1000);
cut400_2x= IQ_corr(1:1600,:,1:500);
%cut500_2x = estImageAll(:,:,1:2:1000);

[size_y size_x ~] = size(cut400_2x);
%%

BaseBandData_Align = CasoratiMtx(cut400_2x);

%for cutOff_Eigenvector = 1:10:200
%disp(num2str(cutOff_Eigenvector));
end_Eigenvector = size(BaseBandData_Align,2);
disp('SVD Filtering');
[SpatialBasis_U EigenValue_S TemporalBasis_V] = svd(BaseBandData_Align,'econ');
%temp = SpatialBasis_U(:,cutOff_Eigenvector:end_Eigenvector) * EigenValue_S(cutOff_Eigenvector:end_Eigenvector,cutOff_Eigenvector:end_Eigenvector)' * TemporalBasis_V(:,cutOff_Eigenvector:end_Eigenvector)';
enS = diag(EigenValue_S);
figure;plot(20*log10(enS));




%%
    nCnt = 0;

for tissueCutOffEnd = 1:4:490
    nCnt=nCnt+1;
    tissueCutOffStart = 1;
    tissueCutOffEnd
    uBCutOffStart = tissueCutOffEnd+1;
    uBCutOffEnd=500;
    noiseCutOffStart = uBCutOffEnd+1;
    %noiseCutOffEnd = size(BaseBandData_Align,2);
    
    tempTissue = SpatialBasis_U(:,tissueCutOffStart:tissueCutOffEnd) * EigenValue_S(tissueCutOffStart:tissueCutOffEnd,tissueCutOffStart:tissueCutOffEnd)' * TemporalBasis_V(:,tissueCutOffStart:tissueCutOffEnd)';
    tempBlood = SpatialBasis_U(:,uBCutOffStart:uBCutOffEnd) * EigenValue_S(uBCutOffStart:uBCutOffEnd,uBCutOffStart:uBCutOffEnd)' * TemporalBasis_V(:,uBCutOffStart:uBCutOffEnd)';
    tempNoise = SpatialBasis_U(:,noiseCutOffStart:end_Eigenvector) * EigenValue_S(noiseCutOffStart:end_Eigenvector,noiseCutOffStart:end_Eigenvector)' * TemporalBasis_V(:,noiseCutOffStart:end_Eigenvector)';
    %
    disp('Reshaping');
    for k = 1:100 %size(tempTissue,2)
        tempImg1 = squeeze(tempTissue(:,k));
        FilteredTissue(:,:,k) = reshape(tempImg1, size_y, size_x);
        
        tempImg1 = squeeze(tempBlood(:,k));
        FiltereduB(:,:,k) = reshape(tempImg1, size_y, size_x);
        
        tempImg1 = squeeze(tempNoise(:,k));
        FilteredNoise(:,:,k) = reshape(tempImg1, size_y, size_x);
    end
    
    FilteredTissueMean(:,:,nCnt) = mean(abs(FilteredTissue),3);
    %FilteredTissueMean(:,:,nCnt) = FilteredTissueMean1/max(FilteredTissueMean1(:));

    FiltereduBMean(:,:,nCnt) = mean(abs(FiltereduB),3);%% <- CHECK THIS PART!!!!!!!!
    %FiltereduBMean(:,:,nCnt) = FiltereduBMean1/max(FiltereduBMean1(:));

    %FilteredNoiseMean = mean(abs(FilteredNoise),3);
%     FilteredNoiseMean = FilteredNoiseMean/max(FilteredNoiseMean(:));
end
%%
load('param.mat');
SoundSpeed = 1540;%Resource.Parameters.speedOfSound;
dataDepth = size(IQ_corr,1);
fs = Receive(1).samplesPerWave * Trans.frequency;    % In MHz
lat = Trans.ElementPos(:,1)'*(SoundSpeed*1e3)/(Trans.frequency*1e6); % Center 128 channel
axial = ((1:(dataDepth))*SoundSpeed/(fs*1e6)/2+Receive(1).startDepth/Trans.frequency*SoundSpeed/1e6)*1e3; %mm % US Pulse-echo case
% polar coordinate
radius=axial(1:1600);
radiusOfCurvature = Trans.radiusMm;
radius = radius + radiusOfCurvature;
theta = Trans.ElementPos(:,4);
[THETA,R] = meshgrid(theta,radius);
[X,Y] = pol2cart(THETA,R);
writerObj = VideoWriter('Trial1_filteredData.avi');
writerObj.FrameRate = 30;
open(writerObj);
%%
for k = 1:300
    figure(1);surf(Y,X-radiusOfCurvature,FiltereduBMean(1:1600,:,k),'edgecolor','none');
    view(0,90);colormap hot;colorbar;axis image;    
    %caxis([3e4 5e5]);
    set(gca,'Ydir','reverse');backColor=[0 0 0];set(gca,'color',backColor);
    title(num2str(((k-1)*4+1)))
    pause(0.001);
    k
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end
close(writerObj);
