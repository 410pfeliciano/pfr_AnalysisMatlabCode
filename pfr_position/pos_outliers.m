
function [Part1A, Part2A, Part3A, Part4A, Part5A] = pos_outliers(Part1,Part2, Part3, Part4, Part5)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 6
      if exist('Part1', 'var')
          figure
          plot(Part1(:,1),Part1(:,2),'.','Color',[75/255,0/255,130/255],'MarkerSize',15);
          var1 = inputname(1);
          title([var1,' ', 'Coordinates in Pixels'])
          hold on
          plotline = plot(1,1,'b');
          button = 1;
          [x y] = deal([]);
          while button == 1
             [xt, yt, button] = ginput(1);
             x = [x , xt];
             y = [y yt];
             set(plotline, 'XData', x, 'YData', y)
          end
          hold off
          dotsin = inpolygon(Part1(:,1),Part1(:,2), x, y);
          Part1A = Part1;
          Part1A(dotsin,:) = NaN;
      else
          Part1 = [];
          
      end
      
      if exist('Part2', 'var')
          figure
          plot(Part2(:,1),Part2(:,2),'.','Color',[255/255,0/255,127/255],'MarkerSize',15)
          var1 = inputname(2);
          title([var1,' ', 'Coordinates in Pixels'])
          hold on
          plotline = plot(1,1,'b');
          button = 1;
          [x y] = deal([]);
          while button == 1
             [xt, yt, button] = ginput(1);
             x = [x , xt];
             y = [y yt];
             set(plotline, 'XData', x, 'YData', y)
          end
          hold off
          dotsin = inpolygon(Part2(:,1),Part2(:,2), x, y);
          Part2A = Part2;
          Part2A(dotsin,:) = NaN;
      else
          Part2 = [];
      end
      
      if exist('Part3', 'var')
          figure
          plot(Part3(:,1),Part3(:,2),'.','Color',[51/255,255/255,153/255],'MarkerSize',15)
          var1 = inputname(3);
          title([var1,' ', 'Coordinates in Pixels'])
          hold on
          plotline = plot(1,1,'b');
          button = 1;
          [x y] = deal([]);
          while button == 1
             [xt, yt, button] = ginput(1);
             x = [x , xt];
             y = [y yt];
             set(plotline, 'XData', x, 'YData', y)
          end
          hold off
          dotsin = inpolygon(Part3(:,1),Part3(:,2), x, y);
          Part3A = Part3;
          Part3A(dotsin,:) = NaN;
      else
          Part3 = [];
      end
      
      if exist('Part4', 'var')
          figure
          plot(Part4(:,1),Part4(:,2),'.','Color',[127/255,0/255,255/255],'MarkerSize',15)
          var1 = inputname(4);
          title([var1,' ', 'Coordinates in Pixels'])
          hold on
          plotline = plot(1,1,'b');
          button = 1;
          [x y] = deal([]);
          while button == 1
             [xt, yt, button] = ginput(1);
             x = [x , xt];
             y = [y yt];
             set(plotline, 'XData', x, 'YData', y)
          end
          hold off
          dotsin = inpolygon(Part4(:,1),Part4(:,2), x, y);
          Part4A = Part4;
          Part4A(dotsin,:) = NaN;
      else
          Part4 = [];
      end
      
      if exist('Part5', 'var')
          figure
          plot(Part5(:,1),Part5(:,2),'.','Color',[255/255,128/255,0/255],'MarkerSize',15)
          var1 = inputname(5);
          title([var1,' ', 'Coordinates in Pixels'])
          hold on
          plotline = plot(1,1,'b');
          button = 1;
          [x y] = deal([]);
          while button == 1
             [xt, yt, button] = ginput(1);
             x = [x , xt];
             y = [y yt];
             set(plotline, 'XData', x, 'YData', y)
          end
          hold off
          dotsin = inpolygon(Part5(:,1),Part5(:,2), x, y);
          Part5A = Part5;
          Part5A(dotsin,:) = NaN;
      else
          Part5 = [];
      end
else
   
end
end

