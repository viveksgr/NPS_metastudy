function plotGlobalMeans(fmridat)
% PLOTGLOBALMEANS   Plot the global mean ±1 SD of each image in fmridat
%
%   plotGlobalMeans(fmridat)
%
% Input:
%   fmridat.dat  – [nVoxels×nImages] data matrix
%
% Produces a single figure showing each image’s mean intensity ±1 SD.

    % compute per-image means and SDs
    globalMean = mean( fmridat.dat, 1, 'omitnan' );   % 1×nImages
    globalSD   = std(  fmridat.dat, 0, 1, 'omitnan' );

    nImages = numel(globalMean);
    x       = 1:nImages;

    % figure('Name','Global Mean ± 1 SD','NumberTitle','off');
    hold on;

    % plot the mean line
    plot(x, globalMean, '-k', 'LineWidth', 1.5);

    % plot mean ± 1 SD as error bars
    errorbar(x, globalMean, globalSD, 'ok', ...
        'MarkerFaceColor','k', ...
        'MarkerSize',6, ...
        'LineWidth',1);

    % xlabel('Image index');
    % ylabel('Global mean intensity');
    % title('Global Mean \pm 1 SD across voxels');
    grid on;
    box on;

    hold off;
end
