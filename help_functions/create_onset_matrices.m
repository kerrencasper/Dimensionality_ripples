function [onsetmat_ret] = create_onset_matrices(data,onsets, RT)

   % for retrieval onsets ending at RT.

    onsetmat_ret = zeros(numel(onsets),5001);
    for itrl=1:numel(onsets)
        to_add = onsets(1,itrl):onsets(1,itrl)+((RT(itrl,1)*data.fsample));
        onsetmat_ret(itrl,1:size(to_add,2)) = to_add;
        clear to_add
    end
        
end