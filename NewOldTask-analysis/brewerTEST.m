z = peaks(100);
 
% choosing a good colormap is especially important for diverging data, that
% is, data that is centred at zero and has minima and maxima below and
% above. cbrewer has a range of nice colormaps.
colors = cbrewer2('div', 'RdBu', 64);
colors = flipud(colors); % puts red on top, blue at the bottom
colormap(colors);
 
% when the data are sequential (eg. only going from 0 to positive, use for
% example colors = cbrewer('seq', 'YlOrRd', 64); or the default parula.
 
subplot(3,3,6); % take a bit more space here because the colorbar also needs to fit in
imagesc(z);
 
% note that imagesc cannot handle unevenly spaced axes. if you want eg. a
% logarithmically scaled colormap, see uimagesc.m from the file exchange
% (also included in fieldtrip)
 
% imagesc automatically flips the y-axis so that the smallest values go on
% top. Set this right if we want the origin to be in the left bottom
% corner.
set(gca, 'ydir', 'normal');
axis square;
 
% add the colorbar, make it prettier
handles = colorbar;
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Label.String = '% change';
drawnow;
 
% this looks okay, but the colorbar is very wide. Let's change that!
% get original axes
axpos = get(gca, 'Position');
cpos = handles.Position;
cpos(3) = 0.5*cpos(3);
handles.Position = cpos;
drawnow;
 
% restore axis pos
set(gca, 'position', axpos);
drawnow;
 
xlabel('SomethingX'); ylabel('SomethingY');
set(gca, 'xtick', 25:25:75, 'ytick', [25:25:75]);