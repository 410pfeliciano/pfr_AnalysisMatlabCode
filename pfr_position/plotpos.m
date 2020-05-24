function plotpos(Pos1, Pos2, Pos3, Pos4, Pos5);
% Position have to be a matrix with 2 columns. Every column represents a
% coordinate.

if nargin < 6
    hold on
    if exist('Pos1', 'var')
    plot(Pos1(:,1), Pos1(:,2), '.','Color',[64/255,224/255,208/255],'MarkerSize',15);
    var1 = inputname(1);
    else
     Pos1 = [];
     var1 = [];
    end
    if exist('Pos2', 'var')
    plot(Pos2(:,1),Pos2(:,2),'.','Color',[255/255,0/255,127/255],'MarkerSize',15);
    var2 = inputname(2);
    else
        Pos2 = [];
        var2 = [];
    end
    
    if exist('Pos3', 'var')
    plot(Pos3(:,1),Pos3(:,2),'.','Color',[51/255,255/255,153/255],'MarkerSize',15);
    var3 = inputname(3);
    else
        Pos3 = [];
        var3 = [];
    end
    
    if exist('Pos4', 'var')
    plot(Pos4(:,1),Pos4(:,2),'.','Color',[127/255,0/255,255/255],'MarkerSize',15);
    var4 = inputname(4);
    else
        Pos4 = [];
        var4 = [];
    end
    
    if exist('Pos5', 'var')
    plot(Pos5(:,1),Pos5(:,2),'.','Color',[255/255,128/255,0/255],'MarkerSize',15);
    var5 = inputname(5);
    else
        Pos5 = [];
        var5 = [];
    end
    xlabel('x coordinates(pixels)')
    ylabel('y coordinates(pixels)')
    % List all possible legend handles and their strings
    legHands = {var1 var2 var3 var4, var5}; 
    % Determine which handles aren't empty
    hasData = ~cellfun(@isempty, legHands); 
    % Produce legend, ignore empties
    legend(legHands{hasData}, 'NumColumns', 3,'Box',...
        'off','Location','northoutside')
    hold off
else
    
end
end