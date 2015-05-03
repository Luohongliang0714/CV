function [ M,corner_count ] = find_N_features( tile,N,M,i,j,corner_count,tile_size )
    %find_N_features return the N best features and appends it to M
    %Sort the non-zero entries of the tile
    [sortedValues,sortIndex] = sort(tile(:),'descend');
    if size(sortIndex,1)>=N % more entries than N found
        maxIndex = sortIndex(1:N);
        maxIndex(sortedValues(1:N)==0)=[]; % get the N largest, get rid of 0 elements
    else % less than N found
        maxIndex = sortIndex(1:end);
        maxIndex(sortedValues(1:end)==0)=[];
    end
    [row,col] = ind2sub(size(tile), maxIndex);
    row=row+(i-1)*tile_size(2);
    col=col+(j-1)*tile_size(1);
    if size(maxIndex,1)~=0 
        % add features to result matrix
        M(:,corner_count+1:corner_count+size(row,1))=[col,row]';
        corner_count=corner_count+size(row,1);
    end
end

