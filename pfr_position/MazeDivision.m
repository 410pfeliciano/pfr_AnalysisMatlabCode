function sites = MazeDivision(position, parts)
% By drawing a polygon you can divide the maze into regions. The function
% will give you the indices of the position coordinates when the animal was
% in a region of interest. 
%Input arguments:
% 1- Position = vector with [x y] coordinates of the animal position
% 2- Parts = integer of the number of divisions that you are interested
%Output
% 1- Struct of cell arays with the indices of each region
          figure
          plot(position(:,1),position(:,2),'.','Color',...
              [178/255,0/255,29/255],'MarkerSize',10);
          hold on
     for ii = 1 :parts-1
          plotline = plot(1,1,'b');
          button = 1;
          [x y] = deal([]);
          while button == 1
             [xt, yt, button] = ginput(1);
             x = [x , xt];
             y = [y yt];
             set(plotline, 'XData', x, 'YData', y)
          end
    
          dotsin = inpolygon(position(:,1),position(:,2), x, y);
          sites.regions{ii} = find(dotsin);
     end
      lastregion = ones(length(position),1);
      allregions = vertcat(sites.regions{1:parts-1});
      lastregion(allregions,1)= 0;
      sites.regions{parts} = find(lastregion);
      col=['r', 'b', 'k','m', 'c'];
      
      for kk = 1:parts
      plot(position(sites.regions{kk},1),position(sites.regions{kk},2),...
          '.','Color',col(kk),'MarkerSize',10)
      end
      hold off
      save('MazeRegions.mat', 'sites');
end