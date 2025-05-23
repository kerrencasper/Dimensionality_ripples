function dpca_plot_iapsRec_4s_2d_2color(data, time, yspan, explVar, compNum, events, signif, marg)

% Modify this function to adjust how components are plotted.
%
% Parameters are as follows:
%   data      - data matrix, size(data,1)=1 because it's only one component
%   time      - time axis
%   yspan     - y-axis spab
%   explVar   - variance of this component
%   compNum   - component number
%   events    - time events to be marked on the time axis
%   signif    - marks time-point where component is significant
%   marg      - marginalization number


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displaying legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(data, 'legend')
    
    % if there is only time and no other parameter - do nothing
    if length(time) == 2
        return

    % if there is one parameter
    elseif length(time) == 3

        numOfStimuli = time(2); % time is used to pass size(data) for legend
        colors = lines(numOfStimuli);
        hold on
        
        for f = 1:numOfStimuli
            plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
            text(1.2, f, ['Stimulus ' num2str(f)])
        end
        axis([0 3 -1 1.5+numOfStimuli])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return        

    % two parameters: stimulus and decision (decision can only have two
    % values)

    % 1=emo, 2=neutral
% for dec = 0:1 % decision, 1=old, 2=new
    elseif length(time) == 4 && time(3) == 2
        numOfStimuli = time(2); % time is used to pass size(data) for legend
        colors = lines(numOfStimuli);
        hold on
        
%         outs = {'CorrRem','CorrFam','Missed','CorrNotRem','FalseRem','FalseFam'};
        outs = {'RHit','Miss','CR','Know'};
        styles = {'-','--',':','-.'};
        for f = 1:numOfStimuli
%         f = 1; %
            plot([0.5 1], [f*2 f*2], 'LineStyle', styles(f), 'color', [0.5 0.5 0.5], 'LineWidth', 2)
            text(1.2, f*2, outs{f})
        end

%         for f = 1:numOfStimuli
% %         f = 1; %
%             plot([0.5 1], [f*2 f*2], 'color', colors(f,:), 'LineWidth', 2)
%             text(1.2, f*2, outs{f})
%         end

%         f = 2;
%         plot([0.5 1], [f f], 'color', colors(f,:), 'LineWidth', 2)
%         text(1.2, f, 'Neutral')
%         end

        plot([0.5 1], [-2 -2], 'Color', rgb('FireBrick'), 'LineWidth', 2)
        plot([0.5 1], [-4 -4], 'Color', rgb('Navy'), 'LineWidth', 2)
        text(1.2, -2, 'Emotional')
        text(1.2, -4, 'Neutral')

%         plot([0.5 1], [-2 -2], 'k', 'LineWidth', 2)
%         plot([0.5 1], [-4 -4], 'k--', 'LineWidth', 2)
%         text(1.2, -2, 'Negative')
%         text(1.2, -4, 'Neutral')
        
        axis([0 3 -4.5 1.5+(numOfStimuli*2)])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        set(gca,'Visible','off')
        return
        
    % other cases - do nothing
    else
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting up the subplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(time)
    time = 1:size(data, ndims(data));
end
axis([time(1) time(end) yspan])
hold on

if ~isempty(explVar)
    title(['dPC ' num2str(compNum) ' [' num2str(explVar,'%.1f') '%]'],'FontWeight','normal')
else
    title(['dPC ' num2str(compNum)])
end

if ~isempty(events)
    plot([events; events], yspan, 'Color', [0.6 0.6 0.6])
end

if ~isempty(signif)
    signif(signif==0) = nan;
%     plot(time, signif + yspan(1) + (yspan(2)-yspan(1))*0.05, 'k', 'LineWidth', 3) % DF changed to move line low on 23Jan2023
    plot(time, signif + yspan(1) + (yspan(2)-yspan(1))*0.0005 - 0.5, 'k', 'LineWidth', 1.5)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data) == 2
    % only time - plot it
    plot(time, squeeze(data(1, :)), 'k', 'LineWidth', 2)

elseif ndims(data) == 3
    % different stimuli in different colours
    numOfStimuli = size(data, 2);
    plot(time, squeeze(data(1,:,:)), 'LineWidth', 2)    

elseif ndims(data) == 4 && size(data,3)==2
    % different stimuli in different colours and binary condition as
    % solid/dashed
%     numOfStimuli = size(data, 2);
%     colors = lines(numOfStimuli);
% 
%     for f=1:numOfStimuli 
%         plot(time, squeeze(data(1, f, 1, :)), 'color', colors(f,:), 'LineWidth', 1.25)
%         plot(time, squeeze(data(1, f, 2, :)), '--', 'color', colors(f,:), 'LineWidth', 1.25)
%     end

    numOfDecisions = size(data, 3);
    colors = [rgb('FireBrick');rgb('Navy')];
%     colors = lines(numOfStimuli);

    for f=1:numOfDecisions 
        plot(time, squeeze(data(1, 1, f, :)), 'color', colors(f,:), 'LineWidth', 1.25)
        plot(time, squeeze(data(1, 2, f, :)), '--', 'color', colors(f,:), 'LineWidth', 1.25)
        plot(time, squeeze(data(1, 3, f, :)), ':', 'color', colors(f,:), 'LineWidth', 1.5)
        try
            plot(time, squeeze(data(1, 4, f, :)), '-.', 'color', colors(f,:), 'LineWidth', 1.5)
        catch
            1;
        end
    end
else
    % in all other cases pool all conditions and plot them in different
    % colours
    data = squeeze(data);
    dims = size(data);
    data = permute(data, [numel(dims) 1:numel(dims)-1]);
    data = reshape(data, size(data,1), []);
    data = data';
    
    plot(time, data, 'LineWidth', 2)    
end