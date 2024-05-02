%% This script takes a human psychophysics data and makes a plot.

DataFile = input('Enter the name of data file to process (without .mat): ','s');
load(sprintf('%s.mat',DataFile));
BehavDat = eval(DataFile);

Params = cell(length(BehavDat),6);
for i=1:length(BehavDat)
	Params(i,:) = strsplit(BehavDat(i).CurrentParam,'.');
end



Results = zeros(5,3);

for i=1:length(Params)
	if str2num(cell2mat(Params(i,2))) == 99
		if BehavDat(i).Error(1) == 0
			Results(1,1) = Results(1,1) +1;
		elseif BehavDat(i).Error(1) == 1
			Results(1,2) = Results(1,2) +1;
		end
	elseif str2num(cell2mat(Params(i,2))) == 1
		if BehavDat(i).Error(1) == 0
			Results(2,1) = Results(2,1) +1;
		elseif BehavDat(i).Error(1) == 1
			Results(2,2) = Results(2,2) +1;
		end
	elseif str2num(cell2mat(Params(i,2))) == 5
		if BehavDat(i).Error(1) == 0
			Results(3,1) = Results(3,1) +1;
		elseif BehavDat(i).Error(1) == 1
			Results(3,2) = Results(3,2) +1;
		end
	elseif str2num(cell2mat(Params(i,2))) == 8
		if BehavDat(i).Error(1) == 0
			Results(4,1) = Results(4,1) +1;
		elseif BehavDat(i).Error(1) == 1
			Results(4,2) = Results(4,2) +1;
		end
	elseif str2num(cell2mat(Params(i,2))) == 12
		if BehavDat(i).Error(1) == 0
			Results(5,1) = Results(5,1) +1;
		elseif BehavDat(i).Error(1) == 1
			Results(5,2) = Results(5,2) +1;
		end
	end
end

Results(:,3) = Results(:,1)./(Results(:,1)+Results(:,2));

figure
Xticks = [1,2,3,4,5];
XtickLabels = {'.5','1','5','8','12'};
plot(Results(:,3),'o-')
ylim([0 1])
set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
xlabel('Semitone Diffefence')
ylabel('Success Rate')
%title(DataFile)
