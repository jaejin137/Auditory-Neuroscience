% Allocate 3-dimensional array for calibration map.
calibmap = zeros(23,5,9);

% Assign frequency values to the first column of each 2D array.
for k=1:9
	for i=1:23
		calibmap(i,1,k) = cTDT(1).freqMappings(i).value;
	end
end

% Assign voltage values for each speaker, frequency, and dB value.
for k=1:9
	for j=2:5
		for i=1:23
			calibmap(i,j,k) = cTDT(k).freqMappings(i).dBMappings(j-1,2);
		end
	end
end

% Generate 2D plot of frequency responses of speakers
figure(1)
for k=1:9
	subplot(3,3,k)
	hfig1 = plot(calibmap(:,2:5,k));
	xlim([1 23])
	xlabel('Frequency [kHz]'); ylabel('Voltage [V]');
	title(sprintf('Speaker %i',k))
end

% Generate surface plot of voltages vs dB numbers vs frequency numbers for each speaker.
figure(2)
for k=1:9
	subplot(3,3,k)
	hfig2 = surf(calibmap(:,2:5,k));
	set(hfig2,'facecolor','interp');
	xlim([1 4]); ylim([1 23]);
	xlabel('dB No.'); ylabel('Freq. No.'); zlabel('voltage')
   	title(sprintf('Speaker %i',k))
end

