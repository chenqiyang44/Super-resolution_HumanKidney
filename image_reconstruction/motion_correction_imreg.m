%% Draw singular value curve
Beamformed_IQ = BmodeData(1:1600,:,1:3000);
[size_z,size_x,N_frames]=size(Beamformed_IQ); 
Casorati = reshape(Beamformed_IQ,[size_z*size_x,N_frames]); 
[U,S,V]=svd(Casorati,'econ');
enS = diag(S);
figure;plot(20*log10(enS));
%% Extract tissue signal for tracking
Cutoff=5;
S([Cutoff:end],[Cutoff:end])=0;
Casorati_f=U*S*V';
Bmode_tissue=reshape(Casorati_f,[size_z,size_x,N_frames]);

%check the images
for ii=1:size(Beamformed_IQ,3)
    imagesc(abs(Bmode_tissue(:,:,ii)))
    title(ii)
    axis off,colormap gray
    caxis([0 0.5e4]),drawnow,pause(0.001)
end

%% Plot Correlation Coefficient with a reference frame
Bmode_samp = Beamformed_IQ(820:1260,40:90,:);
for ii=1:size(Bmode_tissue,3)
   CorrCoeff(ii)=corr2(abs(Bmode_samp(:,:,ii)),abs(Bmode_samp(:,:,100)));
end

figure
plot(CorrCoeff,'LineWidth',2)
axis([0 size(Bmode_tissue,3) 0.2 1]),drawnow
xlabel('Frame Number')
ylabel('Correlation Coefficient with frame 100')
set(gca,'FontSize',14)
grid on
%% Estimate translation
[optimizer,metric]=imregconfig('monomodal');
registration.optimizer.MazximumIterations = 200;
registration.optimizer.MazximumStepLength = 0.001;
for ii=1:size(Beamformed_IQ,3)
   tform=imregtform(abs(Bmode_samp(:,:,ii)),abs(Bmode_samp(:,:,44)),'translation',optimizer,metric);
   X(ii)=tform.T(3,1);
   Z(ii)=tform.T(3,2);
   IQ_corr(:,:,ii)=imwarp(Beamformed_IQ(:,:,ii),tform,'OutputView',imref2d([size(Beamformed_IQ,1) size(Beamformed_IQ,2)]));
   ii
end

% tform = affine2d([ ...
%    cosd(theta) sind(theta) 0;...
%   -sind(theta) cosd(theta) 0; ...
%    x y 1])

figure,hold on
h1=plot(X,'LineWidth',2);
h2=plot(Z,'LineWidth',2);
axis([0 650 -5 5])
grid on
xlabel('Frame Number')
ylabel('Displacement (pixels)')
set(gca,'FontSize',14)
legend([h1,h2],{'lateral displacement','axial displacement'})


% for ii=1:300 %size(Bmode_tissue,3)
%    tform=imregtform(abs(Bmode_samp(:,:,ii)),abs(Bmode_samp(:,:,40)),'similarity',optimizer,metric);
%    IQ_corr(:,:,ii)=imwarp(Beamformed_IQ(:,:,ii),tform,'OutputView',imref2d([size(Beamformed_IQ,1) size(Beamformed_IQ,2)]));
%    ii
% end
% 
% 
%%  plot the comparison 
for ii=1:size(Bmode_tissue,3)
    figure(4);
    subplot(1,2,1)
    imagesc(abs(IQ_corr(:,:,ii)))
    title(ii)
    colormap gray, axis square
    caxis([0 1e4]),drawnow
    subplot(1,2,2)
    imagesc(abs(Beamformed_IQ(:,:,ii)))
    title(ii)    %axis equal,axis off,
    colormap gray,axis square
    caxis([0 1e4]),drawnow
    %pause(0.0001)
end
% 
