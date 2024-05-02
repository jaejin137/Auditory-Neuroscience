function plot_RI_par(clustNum,idx_clust,STRFSigData)

	%% Plot the result of RI
	binSize = .02;
	%nRow = 3;
	%nCol = numel(clustNum);
	nRow = ceil(sqrt(numel(clustNum)));
	nCol = ceil(sqrt(numel(clustNum)));
	%h_RI = figure;
	figure
	%currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),nCol*200,nRow*200])
	%suptitle(sprintf('%s-%s', monkeyName, sessionDate))
	%for s = 1:numel(clustNum)
	for s = 1
		%for i = 1:3
		for i = 3
			if i == 1
				subplot(nRow,nCol,nCol*(i-1)+s)
				histogram(STRFSigData.R1,[-1:binSize:1],'normalization','probability')
				hold on
				histogram(STRFSigData.R1sh,[-1:binSize:1],'normalization','probability')
				xlim([-.5 .5])
				ylim([0 1])
				xlabel('Reliability Index')
				ylabel('Probability')
				legend('non-Shuffled','Shuffled')
				legend('boxoff')
				legend('Location','NorthWest')
				title(sprintf('cluster #%i\nR1, R1shuf',clustNum(s)))
				%drawnow
			elseif i == 2
				subplot(nRow,nCol,nCol*(i-1)+s)
				histogram(STRFSigData.R2,[-1:binSize:1],'normalization','probability')
				hold on
				histogram(STRFSigData.R2sh,[-1:binSize:1],'normalization','probability')
				xlim([-.5 .5])
				ylim([0 1])
				xlabel('Reliability Index')
				ylabel('Probability')
				legend('non-Shuffled','Shuffled')
				legend('boxoff')
				legend('Location','NorthWest')
				title('R2, R2shuf')
				%drawnow
			elseif i == 3
				%subplot(nRow,nCol,nCol*(i-1)+s)
				%subplot(nRow,nCol,idx_clust)
				histogram(STRFSigData.R12,[-1:binSize:1],'normalization','probability')
				hold on
				histogram(STRFSigData.R12sh,[-1:binSize:1],'normalization','probability')
				xlim([-.5 .5])
				ylim([0 1])
				xlabel('Reliability Index')
				ylabel('Probability')
				legend('non-Shuffled','Shuffled')
				legend('boxoff')
				legend('Location','NorthWest')
				%title('R12, R12shuf')
				title(sprintf('Cluster #%i',clustNum(idx_clust)))
				%drawnow
			end
		end
		drawnow
	end
	%tightfig;
	%print(sprintf('RI_%s-%s_%s.png',monkeyName,sessionDate,driveID),'-dpng')


end		% End of function definition
