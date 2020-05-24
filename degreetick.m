function degreetick(ax)
% degreetick places degree symbols on tick labels. Labels 
% dynamically update after zooming. 
% 
%% Syntax 
% 
% degreetick 'x'
% degreetick 'y'
% degreetick 'z'
% 
%% Description 
% 
% degreetick 'x' places degree symbols on x tick labels. 
%
% degreetick 'y' places degree symbols on y tick labels. 
%
% degreetick 'z' places degree symbols on z tick labels. 
% 
%% Example 
% 
% plot(1:10,(1:10).^2)
% degreetick 'y'
% box off
% ylabel 'temperature'
% 
%% Author Info
% This function was written by Chad A. Greene of the University
% of Texas at Austin's Institute for Geophysics, May 2015. 
% http://www.chadagreene.com.

narginchk(1,1)
assert(isnumeric(ax)==0,'Specified axis must be x,y, or z.') 

z = zoom(gcf);
p = pan(gcf);

switch lower(ax(1))
    case 'x'
        SetNewXLabels
        set(z,'ActionPostCallback',@SetNewXLabels);
        set(p,'ActionPostCallback',@SetNewXLabels);
        
    case 'y'
        SetNewYLabels
        set(z,'ActionPostCallback',@SetNewYLabels);
        set(p,'ActionPostCallback',@SetNewYLabels);
        
    case 'z'
        SetNewXLabels
        set(z,'ActionPostCallback',@SetNewZLabels);
        set(p,'ActionPostCallback',@SetNewZLabels);
        
    otherwise
        error('Specified axis must be x, y, or z.')
end


end

function SetNewXLabels(hObject, evt)
    % Get tick values: 
    tickvals = get(gca,'xtick');

    % Preallocate empty array for new tick labels: 
    NewTickLabels = cell(1,length(tickvals)); 

    % Create new tick labels with degree symbol:
    for k=1:numel(tickvals);
        NewTickLabels{k}=sprintf('%g°',tickvals(k));
    end

    % set new tick labels: 
    set(gca,'xticklabel',NewTickLabels);
end

function SetNewYLabels(hObject, evt)
    % Get tick values: 
    tickvals = get(gca,'Ytick');

    % Preallocate empty array for new tick labels: 
    NewTickLabels = cell(1,length(tickvals)); 

    % Create new tick labels with degree symbol:
    for k=1:numel(tickvals);
        NewTickLabels{k}=sprintf('%g°',tickvals(k));
    end

    % set new tick labels: 
    set(gca,'Yticklabel',NewTickLabels);
end

function SetNewZLabels(hObject, evt)
    % Get tick values: 
    tickvals = get(gca,'Ztick');

    % Preallocate empty array for new tick labels: 
    NewTickLabels = cell(1,length(tickvals)); 

    % Create new tick labels with degree symbol:
    for k=1:numel(tickvals);
        NewTickLabels{k}=sprintf('%g°',tickvals(k));
    end

    % set new tick labels: 
    set(gca,'Zticklabel',NewTickLabels);
end



    


