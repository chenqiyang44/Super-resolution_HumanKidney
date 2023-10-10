%%
Beamformed_IQ = IQ_corr(1:1600,:,1:500);
[size_z,size_x,N_frames]=size(Beamformed_IQ); 
Casorati = reshape(Beamformed_IQ,[size_z*size_x,N_frames]); 
[U,S,V]=svd(Casorati,'econ');
enS = diag(S);
figure;plot(20*log10(enS));
%%
Cutoff=65;
Cutoff_noise = 500;
S([1:Cutoff],[1:Cutoff])=0;
%S([Cutoff_noise:N_frames],[Cutoff_noise:N_frames])=0;
Casorati_f=U*S*V';
FilteredData=reshape(Casorati_f,[size_z,size_x,N_frames]);

for ii=1:size(FilteredData,3)
    figure(19)
    imagesc(abs(FilteredData(:,:,ii)))
    title(ii)
    colormap hot;axis off
    caxis([256 1134]),drawnow,pause(0.01)
end
figure;
imagesc(mean(abs(FilteredData),3));colormap hot
%%
for i=1:2:300
    a=abs(FilteredData(1:1500,:,i));
    a = a/max(a(:));
    a(a<0.3) = 0;
    figure(1);
    imagesc(a); colormap(hot);
    drawnow
    pause(0.1)
end

%%
% n = 1;
% 
% for i = 1:size(FilteredData,3)
%      if CorrCoeff(i) > 0.9
%          temp = abs(FilteredData(:,:,i));
%          temp = temp/max(temp(:));
%          temp(temp<0.05) = 0;
%          testdata(:,:,n) = temp;
%          n = n + 1;
%      end
% end

for i = 1:size(FilteredData,3)

    temp = abs(FilteredData(:,:,i));
    %temp = temp/max(temp(:));
    temp(temp<400) = 0;
    temp = temp/2000;
    temp(temp>1) = 1;
    %temp = temp/max(temp(:));
    testdata(:,:,i) = temp;
%     figure(1);
%     imagesc(temp); colormap(hot);%caxis([0 2000]);
%     drawnow
%     pause(0.01)
    i
end
save(['0402_humankidney_filtered_400globalnorm.mat'],'testdata');

