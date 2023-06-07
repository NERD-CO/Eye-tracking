time = 0:0.01:10; % seconds, sampled at 100 Hz
data(:, :, 1) = bsxfun(@plus, sin(time), randn(100, length(time)));
data(:, :, 2) = bsxfun(@plus, cos(time), randn(100, length(time)));
 
colors = cbrewer2('qual', 'Set2', 8);
 
subplot(4,4,[13 14]);  % plot across two subplots
hold on;
bl = boundedline(time, mean(data(:, :, 1)), std(data(:, :, 1)), ...
    time, mean(data(:, :, 2)), std(data(:, :, 2)), ...
    'cmap', colors);
% boundedline has an 'alpha' option, which makes the errorbars transparent
% (so it's nice when they overlap). However, when saving to pdf this makes
% the files HUGE, so better to keep your hands off alpha and make the final
% figure transparant in illustrator
 
xlim([-0.4 max(time)]); xlabel('Time (s)'); ylabel('Signal');
 
% instead of a legend, show colored text
lh = legend(bl);
legnames = {'sin', 'cos'};
for i = 1:length(legnames),
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
 
% move a bit closer
lpos = lh.Position;
lpos(1) = lpos(1) + 0.15;
lh.Position = lpos;
 
