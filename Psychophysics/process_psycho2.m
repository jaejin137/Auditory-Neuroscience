%% This script takes a human psychophysics data and makes a plot.

DataFile = input('Enter the name of data file to process (without .mat): ','s');
load(sprintf('%s.mat',DataFile));
BehavDat = eval(DataFile);

Params = cell(length(BehavDat),6);
ReactTime = zeros(length(BehavDat),1);
for i=1:length(BehavDat)
	Params(i,:) = strsplit(BehavDat(i).CurrentParam,'.');
	ReactTime(i,1) = BehavDat(i).LeverAcquired(2) - BehavDat(i).TimeStamp(1);
end

% Trial numbers for each cases.
trialnum_p5_succ = [];
trialnum_p5_fail = [];
trialnum_01_succ = [];
trialnum_01_fail = [];
trialnum_05_succ = [];
trialnum_05_fail = [];
trialnum_08_succ = [];
trialnum_08_fail = [];
trialnum_12_succ = [];
trialnum_12_fail = [];

for i=1:length(Params)
	if str2num(cell2mat(Params(i,2))) == 99
		if BehavDat(i).Error(1) == 0
			trialnum_p5_succ(end+1) = i;
		elseif BehavDat(i).Error(1) == 1
			trialnum_p5_fail(end+1) = i;
		end
	elseif str2num(cell2mat(Params(i,2))) == 1
		if BehavDat(i).Error(1) == 0
			trialnum_01_succ(end+1) = i;
		elseif BehavDat(i).Error(1) == 1
			trialnum_01_fail(end+1) = i;
		end
	elseif str2num(cell2mat(Params(i,2))) == 5
		if BehavDat(i).Error(1) == 0
			trialnum_05_succ(end+1) = i;
		elseif BehavDat(i).Error(1) == 1
			trialnum_05_fail(end+1) = i;
		end
	elseif str2num(cell2mat(Params(i,2))) == 8
		if BehavDat(i).Error(1) == 0
			trialnum_08_succ(end+1) = i;
		elseif BehavDat(i).Error(1) == 1
			trialnum_08_fail(end+1) = i;
		end
	elseif str2num(cell2mat(Params(i,2))) == 12
		if BehavDat(i).Error(1) == 0
			trialnum_12_succ(end+1) = i;
		elseif BehavDat(i).Error(1) == 1
			trialnum_12_fail(end+1) = i;
		end
	end
end

% Put all trial numbers together in Results cell array.
% Format of Results: successful trial numbers, failed trial numbers, success rate,
%					 , average reaction time, and s.d. of reaction time in each column,
%					 different semitone differences in each row
Results = cell(5,4);

% Trial numbers for each cases.
Results{1,1} = trialnum_p5_succ;
Results{1,2} = trialnum_p5_fail;
Results{2,1} = trialnum_01_succ;
Results{2,2} = trialnum_01_fail;
Results{3,1} = trialnum_05_succ;
Results{3,2} = trialnum_05_fail;
Results{4,1} = trialnum_08_succ;
Results{4,2} = trialnum_08_fail;
Results{5,1} = trialnum_12_succ;
Results{5,2} = trialnum_12_fail;

for i=1:5
	% Percent success rate
	Results{i,3} = length(Results{i,1})./(length(Results{i,1})+length(Results{i,2}));
	% Average of reaction time
	Results{i,4} = mean(ReactTime(Results{i,1}));
	% s.d. of reaction time
	Results{i,5} = std(ReactTime(Results{i,1}));
end

% Plot percent success rate vs semitone differences.
figure(1)
Xticks = [1,2,3,4,5];
XtickLabels = {'.5','1','5','8','12'};
plot(cell2mat(Results(:,3)),'o-')
xlim([.9 5.1])
ylim([0 1])
set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
xlabel('Semitone Diffefence')
ylabel('Success Rate')
hold on

% Plot reaction time vs semitone differences.
figure(2)
Xticks = [1,2,3,4,5];
XtickLabels = {'.5','1','5','8','12'};
errorbar(cell2mat(Results(:,4)),cell2mat(Results(:,5)),'o-')
xlim([.9 5.1])
ylim([0 max(cell2mat(Results(:,4)))+max(cell2mat(Results(:,5)))*1.1])
set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
xlabel('Semitone Diffefence')
ylabel('Reaction Time')
hold on


