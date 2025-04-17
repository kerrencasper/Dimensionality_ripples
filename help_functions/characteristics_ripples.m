function [rippletime_ret, rippledur_ret, rippleamp_ret, ripplefreq_ret] = characteristics_ripples(alldat, allripples_samples, onsetmat_ret, onsets)

    % for retrieval

    rippletime_ret_tmp  = [];
    rippledur_ret_tmp   = [];
    rippleamp_ret_tmp   = [];
    ripplefreq_ret_tmp  = [];
    rippletime_ret      = cell(numel(alldat)-1,1);
    rippledur_ret       = cell(numel(alldat)-1,1);
    rippleamp_ret       = cell(numel(alldat)-1,1);
    ripplefreq_ret      = cell(numel(alldat)-1,1);

    for ichannel = 1:numel(alldat)
        ichannel
        tmp = allripples_samples(ichannel,:); % All timestamps per channel
        tmp(tmp==0)=[];
        for nucol = 1:size(tmp,2)
            time_rip = tmp(1,nucol); % single timestamp
            for itrl=1:numel(onsets)
                [~,b] = find(onsetmat_ret(itrl,:)==time_rip); % Find trial of timestamp
                if isempty(b)
                    time_stamp  = NaN; % If not that trial, NaN it.
                    dur         = NaN;
                    amp         = NaN;
                    freq        = NaN;
                else
                    [~,tmp_b]   = find(alldat{ichannel}.evtIndiv.maxTime == time_rip); % If a trial, find the idx of that ripple in the data.
                    time_stamp  = alldat{ichannel}.evtIndiv.maxTime(1,tmp_b); % Add timestamp of ripple.
                    dur         = alldat{ichannel}.evtIndiv.duration(1,tmp_b); % Use idx to find the duration of that particular ripple
                    amp         = alldat{ichannel}.evtIndiv.maxAmp(1,tmp_b); % Use idx to find the amplitude of that particular ripple
                    freq        = alldat{ichannel}.evtIndiv.freq(1,tmp_b); % Use idx to find the frequency of that particular ripple
                end
                [rippletime_ret_tmp(itrl,:)]    = time_stamp-onsets(1,itrl);
                [rippledur_ret_tmp(itrl,:)]     = dur;
                [rippleamp_ret_tmp(itrl,:)]     = amp;
                [ripplefreq_ret_tmp(itrl,:)]    = freq;
                clear dur b time_stamp amp freq
            end
            rippletime_ret_col(:,nucol)     = rippletime_ret_tmp;
            rippledur_ret_col(:,nucol)      = rippledur_ret_tmp;
            rippleamp_ret_col(:,nucol)      = rippleamp_ret_tmp;
            ripplefreq_ret_col(:,nucol)     = ripplefreq_ret_tmp;
            clear rippletime_ret_tmp rippledur_ret_tmp...
                time_rip rippleamp_ret_tmp ripplefreq_ret_tmp
        end
        clear tmp
        rippletime_col  = [];
        rippledur_col   = [];
        rippleamp_col   = [];
        ripplefreq_col  = [];
        for nurows = 1:size(rippletime_ret_col,1) % Reduce matrix
            [~,b] = find(rippletime_ret_col(nurows,:)>0);
            if isempty(b)
                [rippletime_col(nurows,1)]  = NaN;
                [rippledur_col(nurows,1)]   = NaN;
                [rippleamp_col(nurows,1)]   = NaN;
                [ripplefreq_col(nurows,1)]  = NaN;
            else
                [rippletime_col(nurows,1:size(b,2))]    = rippletime_ret_col(nurows,b);
                [rippledur_col(nurows,1:size(b,2))]     = rippledur_ret_col(nurows,b);
                [rippleamp_col(nurows,1:size(b,2))]     = rippleamp_ret_col(nurows,b);
                [ripplefreq_col(nurows,1:size(b,2))]    = ripplefreq_ret_col(nurows,b);
            end
            clear b
        end
        clear rippletime_ret_col rippledur_ret_col rippleamp_ret_col ripplefreq_ret_col
        rippletime_col(rippletime_col==0)   = NaN;
        rippletime_ret{ichannel,:}          = rippletime_col;
        rippletime_col                      = [];
        rippledur_col(rippledur_col==0)     = NaN;
        rippledur_ret{ichannel,:}           = rippledur_col;
        rippledur_col                       = [];
        rippleamp_col(rippleamp_col==0)     = NaN;
        rippleamp_ret{ichannel,:}           = rippleamp_col;
        rippleamp_col                       = [];
        ripplefreq_col(ripplefreq_col==0)   = NaN;
        ripplefreq_ret{ichannel,:}          = ripplefreq_col;
        ripplefreq_col                      = [];
        clear tmp

    end
    
end