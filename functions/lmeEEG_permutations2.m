function [rperms] = lmeEEG_permutations2(nperms,SS,Item)
% lmeEEG_permutations: Compute unique permutations within subjects and items
% This function is for fully-crossed designs  
% [Inputs]
% - SS: Nominal variable of subjects' indeces
% - Item: Nominal variable of items' indeces
% - nperms: number of permutations
% [Output]
% -rperms: matrix with unique permutations
disp('Compute permutation indeces')
subj =unique(SS);
it = unique(Item);
rperms = nan(length(SS),nperms);

for id = 1:length(subj)
    for itx = 1:length(it)
        idx = find(SS==subj(id)&Item==it(itx));
        idx2 = repmat(idx,1,nperms);
        for i =1:nperms
            goon = true;
            while goon
                tmp = idx2(randperm(size(idx2,1)),i);
                if isempty(find(sum(tmp==idx2)==size(idx2,1)))
                    idx2(:,i)=tmp;
                    goon = false;
                end
            end
        end
        progressbar(id/length(subj));
        rperms(idx,:)=idx2;
    end
end

