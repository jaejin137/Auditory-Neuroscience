if ~exist('plotrange')
	display('No "plotrange" given. Default plot range (1~10) automatically selected.')
	display('Help: plotrange = [min max]')
	plotrange = [1 10];
end

% Determine what ylim should be like by taking maximum of Raw, MUA, and MUAHT.
NumTrial=max(size(NeuralRawTr));
XMax = 1;
YMax = 1;
YMaxHT = 1;
YLim = 1.0e+3;
YLimHT = 1.0e+2;
for i=1:NumTrial
	if max(size(NeuralRawTr{1,i})) > XMax
		XMax = max(size(NeuralRawTr{1,i}));
	end
	if max(NeuralRawTr{1,i}) > YMax & max(NeuralRawTr{1,i}) < YLim
		YMax = max(NeuralRawTr{1,i});
	end
	if max(NeuralMUAHTTr{1,i}) > YMaxHT & max(NeuralMUAHTTr{1,i}) < YLimHT
		YMaxHT = max(NeuralMUAHTTr{1,i});
	end
end

% Plot Raw, MUA, and Hilbert Transform of MUA in an array.
figure
numrow = ceil((max(plotrange)-min(plotrange)+1)/2);
for i=min(plotrange):min(max(plotrange),max(size(NeuralRawTr)))
    subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+1);
    plot(eval(sprintf('NeuralRawTr{1,%i}',i)));
    hold on;
    plot(eval(sprintf('NeuralMUATr{1,%i}',i)),'g');
	xlim([0 XMax])
	if max(NeuralRawTr{1,i}) <= YLim
		ylim([-YMax YMax])
	end
    legend('Raw','MUA');
    title(sprintf('NeuralTr%i',i));
    subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+2);
    plot(abs(eval(sprintf('NeuralMUAHTTr{1,%i}',i))));
	xlim([0 XMax])
	if max(NeuralMUAHTTr{1,i}) <= YLimHT
		ylim([0 YMaxHT])
	end
    title('Hilbert Transform');
    subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+3);
    periodogram(eval(sprintf('NeuralMUATr{1,%i}',i)));
end
