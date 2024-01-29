function [rperms] = lmeEEG_permutations(SS,nperms)
% lmeEEG_permutations: Compute unique permutations within subjects
% This function is for stimuli-within-condition designs
% [Inputs]
% - SS: Nominal variable of subjects' indeces
% - nperms: number of permutations
% [Output]
% -rperms: matrix with unique permutations

subj =unique(SS);
disp('Compute permutation indeces')

rperms = nan(length(SS),nperms);
for id = 1:length(subj)
    idx = find(SS==subj(id));
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

