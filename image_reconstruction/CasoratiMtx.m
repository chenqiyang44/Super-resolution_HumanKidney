function [ outputMtx ] = CasoratiMtx( inputMtx )

% inputMtx (Y, X, t);
[size_y size_x size_t] = size(inputMtx);
for k = 1:size_t
    temp1 = squeeze(inputMtx(:,:,k));
    outputMtx(:,k) = reshape(temp1, size_x*size_y, 1);
end

