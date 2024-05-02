%% Plot the result of RI
binSize = .02;

clustNum = [STRFSig.clustNum]';
nRow = floor(sqrt(numel(clustNum)));
nCol = ceil(sqrt(numel(clustNum)));
if nRow*nCol < numel(clustNum)
	nRow = ceil(sqrt(numel(clustNum)));
end
	
h_RI = figure;
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),nCol*200,nRow*200])
%suptitle(sprintf('%s-%s', monkeyName, sessionDate))
for s = 1:length(clustNum)
%for s = 1
	try
		%for i = 1:3
		for i = 3
			if i == 1
				%STRFSigTemp = STRFSigData{s,4};
				STRFSigTemp = STRFSig(s).STRFSigData;
				%[N,X] = hist(STRFSigTemp.R1,[-1:.02:1]);
				%[Nshuf,X] = hist(STRFSigTemp.R1sh,[-1:.02:1]);
				subplot(nRow,nCol,nCol*(i-1)+s)
				%bar(X,N/NB);
				histogram(STRFSigTemp.R1,[-1:binSize:1],'normalization','probability')
				hold on
				%bar(X,Nshuf/NB,'r')
				histogram(STRFSigTemp.R1sh,[-1:binSize:1],'normalization','probability')
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
				%STRFSigTemp = STRFSigData{s,4};
				STRFSigTemp = STRFSig(s).STRFSigData;
				%[N,X] = hist(STRFSigTemp.R2,[-1:.02:1]);
				%[Nshuf,X] = hist(STRFSigTemp.R2sh,[-1:.02:1]);
				subplot(nRow,nCol,nCol*(i-1)+s)
				%bar(X,N/NB);
				histogram(STRFSigTemp.R2,[-1:binSize:1],'normalization','probability')
				hold on
				%bar(X,Nshuf/NB,'r')
				histogram(STRFSigTemp.R2sh,[-1:binSize:1],'normalization','probability')
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
				%STRFSigTemp = STRFSigData{s,4};
				STRFSigTemp = STRFSig(s).STRFSigData;
				%[N,X] = hist(STRFSigTemp.R12,[-1:.02:1]);
				%[Nshuf,X] = hist(STRFSigTemp.R12sh,[-1:.02:1]);
				%subplot(nRow,nCol,nCol*(i-1)+s)
				subplot(nRow,nCol,s)
				%bar(X,N/NB);
				histogram(STRFSigTemp.R12,[-1:binSize:1],'normalization','probability')
				hold on
				%bar(X,Nshuf/NB,'r')
				histogram(STRFSigTemp.R12sh,[-1:binSize:1],'normalization','probability')
				xlim([-.5 .5])
				ylim([0 1])
				xlabel('Reliability Index')
				ylabel('Probability')
				legend('non-Shuffled','Shuffled')
				legend('boxoff')
				legend('Location','NorthWest')
				%title('R12, R12shuf')
				title(sprintf('Cluster #%i (p=%.2d)',clustNum(s),STRFSig(s).p))
				%drawnow
			end
		end
		drawnow
	catch
	end
end
%tightfig;
%print(sprintf('RI_%s-%s_%s.png',monkeyName,sessionDate,driveID),'-dpng')
