thresh = [65 70 122 100 81 88];
for i = 1:6
    Beamformed_IQ = IQ_corr(1:1600,:,(i-1)*500+1:(i-1)*500+500);
    [size_z,size_x,N_frames]=size(Beamformed_IQ); 
    Casorati = reshape(Beamformed_IQ,[size_z*size_x,N_frames]); 
    [U,S,V]=svd(Casorati,'econ');
    enS = diag(S);
    figure;plot(20*log10(enS));
    %%
    Cutoff=thresh(i);
    S([1:Cutoff],[1:Cutoff])=0;
    %S([300:end],[300:end])=0;
    Casorati_f=U*S*V';
    FilteredData=reshape(Casorati_f,[size_z,size_x,N_frames]);
    
    fname =  ['1_500ensemble_corr_bloc' num2str(i) '_filtered'];
    save(fname,'FilteredData');
    i
end
