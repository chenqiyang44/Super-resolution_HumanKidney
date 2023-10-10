clear
Threshold=0.1;
m=linspace(1,1600,1600);
n0=linspace(1,128,128);
n=linspace(1,128,1024);
addpath(genpath('SimpleTracker'))
%% loacalization and tracking
for No_blocs = 1:6
    tic
    filename_load = ['2_500ensemble_corr_bloc' num2str(No_blocs) '_filtered'];
    load(filename_load);
    for i = 1:size(FilteredData,3)
        prediction(:,:,i) = abs(interp2(n0,m',FilteredData(:,:,i),n,m','spline'));
    end
    [size_z,size_x,N_frames]=size(prediction(:,:,:)); 
    [Pix_X,Pix_Z]=meshgrid([-2:2],[-2:2]);
    for nn=1:N_frames
        Im=squeeze(prediction(:,:,nn));
        %Im = imgaussfilt(Im,1,'filtersize',3);
        Im = Im/max(Im(:));
        Reg_max=Im.*imregionalmax(Im);
        Reg_max(Reg_max<Threshold)=0;
        [Ind_z,Ind_x]=ind2sub(size(Reg_max),find(Reg_max));    
        kk = 1;
        for i_loc=1:length(Ind_z)
            Z_init=Ind_z(i_loc);
            X_init=Ind_x(i_loc);
            if Z_init>2 && Z_init<1599 && X_init>2 && X_init<1023
                Im_MB=Im(Z_init-2:Z_init+2,X_init-2:X_init+2);
                dz=sum(sum(Im_MB.*Pix_Z))/sum(sum(Im_MB));
                dx=sum(sum(Im_MB.*Pix_X))/sum(sum(Im_MB));
                Z_waverage=Z_init+dz+0.5;
                X_waverage=X_init+dx+0.5;
                SR_Localizations{nn}(kk,:)=[Z_waverage X_waverage];
                kk = kk+1;
            end
        end
    end


    Points=cell2mat(SR_Localizations');
    ind=sub2ind([size_z,size_x],round(Points(:,1)),round(Points(:,2)));
    Imfinal=zeros([size_z, size_x]);
    for pp=1:length(ind)
        Imfinal(ind(pp))=Imfinal(ind(pp))+1;
    end

    figure
    imagesc(Imfinal)
    colormap gray
    axis image,axis off
    caxis([0 2]),drawnow
    
    % Tracking
    max_linking_distance=12;
    max_gap_closing=1;
    [tracks,adjacency_tracks] = simpletracker(SR_Localizations,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', false);

    filename_save = ['2_500ensemble_corr_thresh01_bloc' num2str(No_blocs) '.mat'];
    save(filename_save,'adjacency_tracks','SR_Localizations');
    clear SR_Localizations prediction
    No_blocs
    toc
end

