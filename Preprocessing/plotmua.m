%% Plot some basic figures for the first 10 trials.
if ~exist('plotrange')
	display('No "plotrange" given. Default plot range (1~10) automatically selected.')
	display('Help: plotrange = [min max]')
	plotrange = [1 10];
end

% Determine what ylim should be like by taking maximum of Raw, MUA, and MUAHT.
xmax = 1;
for i=1:num_trial
	if max(size(neural_raw_tr{1,i})) > xmax
		xmax = max(size(neural_raw_tr{1,i}));
	end
end

% Plot Raw, MUA, and Hilbert Transform of MUA in an array.
figure
numrow = ceil((max(plotrange)-min(plotrange)+1)/2);
for i=min(plotrange):min(max(plotrange),max(size(neural_raw_tr)))
    subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+1);
    plot(eval(sprintf('neural_raw_tr{1,%i}',i)));
    hold on;
    plot(eval(sprintf('neural_mua_tr{1,%i}',i)),'g');
	xlim([0 xmax])
	ylim([-150 150])
    legend('Raw','MUA');
    title(sprintf('Neural Signal Tr%i',i));
    subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+2);
    plot(abs(eval(sprintf('neural_mua_ht_tr{1,%i}',i))));
	xlim([0 xmax])
	ylim([0 50])
    title('Hilbert Transform');
    subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+3);
    periodogram(eval(sprintf('neural_mua_tr{1,%i}',i)));
end