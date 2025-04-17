function [artitime_ret] = extract_artifacts(data,alldat,onsets, onsetmat_ret)

mat_arti = zeros(numel(alldat),numel(data.time{1}));

for ichannel = 1:numel(alldat)
    ichannel
    
    % mark all artifact including padding period
    for iartifact = 1:size(data.artifact{ichannel},1)
        idx  = data.artifact{ichannel}(iartifact,1)+alldat{ichannel}.param.artfctPad(1)*data.fsample : data.artifact{ichannel}(iartifact,2)+alldat{ichannel}.param.artfctPad(2)*data.fsample;
        idx = idx(idx>0 & idx <= numel(mat_arti(ichannel,:)));
        
        mat_arti(ichannel,idx) = 1;
    end
    
    % mark all clean segments shorter than minimum ripple duration as
    % unavailable
    transitions      = [mat_arti(ichannel,1) diff(mat_arti(ichannel,:))];
    clean_beginnings = find(transitions == 1);
    clean_endings    = find(transitions == -1);
    
    for iblock = 1:numel(clean_beginnings)-1
        
        if clean_endings(iblock)-clean_beginnings(iblock) < alldat{ichannel}.criterion.len(1) * data.fsample
            mat_arti(ichannel,clean_beginnings(iblock):clean_endings(iblock)) = 1;
        end
    end
    
end  
    
        %% Find when the artifact happened by extracting the time points
    for ichannel = 1:numel(alldat)
        tmp = mat_arti(ichannel,:);
        allarti_samples(ichannel,1:size(find(tmp),2)) = find(tmp);
        clear tmp
    end

    %% for retrieval
    artitime_ret_tmp = [];
    for ichannel = 1:numel(alldat)
        tmp = allarti_samples(ichannel,:);
        tmp(tmp==0)=NaN;
        for itrl=1:numel(onsets)
            [~,b] = find(ismember(tmp,onsetmat_ret(itrl,:)));
            when = tmp(1,b)-onsets(1,itrl);
            %     when = allripples_samples(1,b);
            [artitime_ret_tmp(itrl,1:numel(when))] = when;
        end
        artitime_ret_tmp(artitime_ret_tmp==0)=NaN; % NaN zeros
        if isempty(artitime_ret_tmp)
            artitime_ret_tmp = NaN(size(onsetmat_ret,1),1);
        end
        artitime_ret{ichannel,:} = artitime_ret_tmp;
        clear artitime_ret_tmp when b tmp
    end

    
    
end