%% This script takes a human psychophysics data and makes a plot.

%% Clear
clear
%close all

%DataFile = input('Enter the name of data file to process (without .mat): ','s');
%load(sprintf('%s.mat',DataFile));
%BehavDat = eval(DataFile);

key_common = input('Enter common keyword: ','s');
cmd_list = sprintf('ls %s*',key_common);
[status,list] = system(cmd_list);
target = sort(strsplit(list));
target = target(2:end);
for k=1:length(target)
	load(target{k})
	dataname = strrep(target{k},'.mat','');
	if k==1
		BehavDat = eval(dataname);
	else
		BehavDat = [BehavDat,eval(dataname)];
	end
end

%% Define a cell for Summary of all the trials
%% column 1-6: parameters
%% column 7: result (h: hit, r: correct rejection, m: miss, f: false alarm)
%% column 8: reaction time.
Summary = cell(length(BehavDat),8);

% parameter indices.
dfInd = {'91','92','93','95','96','97','99','01','02','03','05','07','09','10','12'};
angInd = [0, 15, 30, 45, 60];
% Define matrices for Results: hit, correct rejection, misses, and false alarm.
% dimension: 8 by 5 for 7 tone differences(.5,1,3,5,7,9,12) and 5 angular separation(0,15,30,45,60deg). 
numDf = length(dfInd);
numAng = length(angInd);
Hit = zeros(numDf,numAng);
Rej = zeros(numDf,numAng);
Mis = zeros(numDf,numAng);
FA = zeros(numDf,numAng);
RT = cell(numDf,numAng);

for i=1:length(BehavDat)
	Summary(i,1:6) = strsplit(BehavDat(i).CurrentParam,'.');
	% Find indices for Results cell.
	if Summary{i,2} == '91'		% .25 semitone
		l = 1;
	elseif Summary{i,2} == '92'	% .5 semitone
		l = 2;
	elseif Summary{i,2} == '93'	% .75 semitone
		l = 3;
	elseif Summary{i,2} == '01'	% 1 semitone
		l = 4;
	%elseif Summary{i,2} == '95' % 1.25 semitone
	%	l = 5;
	elseif Summary{i,2} == '96' % 1.5 semitone
		l = 5;
	%elseif Summary{i,2} == '97' % 1.75 semitone
	%	l = 7;
	elseif Summary{i,2} == '02'	% 2 semitone
		l = 6;
	elseif Summary{i,2} == '99'	% 2.5 semitone
		l = 7;
	elseif Summary{i,2} == '03'	% 3 semitone
		l = 8;
	elseif Summary{i,2} == '05'	% 5 semitone
		l = 9;
	elseif Summary{i,2} == '07'	% 7 semitone
		l = 10;
	elseif Summary{i,2} == '09'	% 9 semitone
		l = 11;
	elseif Summary{i,2} == '10' % 10 semitone
		l = 12;
	elseif Summary{i,2} == '12'	% 12 semitone
		l = 13;
	else
		fprintf('Error found on Df values: %s\n', Summary{i,2})
	end
	m = 1+abs(str2num(Summary{i,4})-str2num(Summary{i,5}));
	% Count the number of occurances for each cases and assign some basic values.
	if BehavDat(i).Error == [0 0]
		if Summary{i,6} == '0'
			Summary{i,7} = 'h';
			Summary{i,8} = BehavDat(i).JoystickAcquired(2) - BehavDat(i).TimeStamp(1);
			Hit(l,m) = Hit(l,m)+1;
		else
			Summary{i,7} = 'r';
			Rej(l,m) = Rej(l,m)+1;
		end
	else
		if Summary{i,6} == '0'
			Summary{i,7} = 'm';
			Mis(l,m) = Mis(l,m)+1;
		else
			Summary{i,7} = 'f';
			FA(l,m) = FA(l,m)+1;
		end
	end
end


%% Basic statistics for each parameters.
% Third dimension: success rate, mean reaction time, s.d. of reaction time, ratio of false alarm to catch trials,  
Stat = nan(numDf,numAng,4);
for j=1:numAng
	for i=1:numDf
		% Percent success rate
		if Hit(i,j)+Rej(i,j)+Mis(i,j)+FA(i,j) ~= 0
			Stat(i,j,1) = (Hit(i,j)+Rej(i,j))/(Hit(i,j)+Rej(i,j)+Mis(i,j)+FA(i,j));
		end
		% Average of reaction time
		if Hit(i,j)+Rej(i,j)+Mis(i,j)+FA(i,j) ~= 0 && Hit(i,j) ~= 0
			Stat(i,j,2) = mean(cell2mat(Summary(logical((~cellfun(@isempty,strfind(Summary(:,7),'h')).*~cellfun(@isempty,strfind(Summary(:,2),dfInd{i}))).*(abs(cell2mat(Summary(:,4))-cell2mat(Summary(:,5)))+1 == j)),8)));
			Stat(i,j,3) = std(cell2mat(Summary(logical((~cellfun(@isempty,strfind(Summary(:,7),'h')).*~cellfun(@isempty,strfind(Summary(:,2),dfInd{i}))).*(abs(cell2mat(Summary(:,4))-cell2mat(Summary(:,5)))+1 == j)),8)));
		end
		% Ratio of # of false alarms to # of catch trials
		if Rej(i,j)+FA(i,j) ~= 0
			Stat(i,j,4) = FA(i,j)/(Rej(i,j)+FA(i,j));
		end
	end
end

%% Plot the result
if strfind(key_common,'_ton_')
	% Plot percent success rate vs semitone differences.
	figure
	Xticks = 1:numDf;
	XtickLabels = {'.25','.5','.75','1','1.5','2','2.5','3','5','7','9','10','12'};
	h1 = plot(Stat(:,:,1),'o-');
	xlim([1 numDf])
	ylim([0 1])
	set(h1,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Semitone Diffefence')
	ylabel('Success Rate')
	hold on
	% Plot reaction time vs semitone differences.
	figure
	Xticks = 1:numDf;
	XtickLabels = {'.25','.5','.75','1','1.5','2','2.5','3','5','7','9','10','12'};
	h2 = errorbar(Stat(:,:,2),Stat(:,:,3),'o-');
	xlim([1 numDf])
	ylim([0 max(max(Stat(:,:,2)))+max(max(Stat(:,:,3)))*1.1])
	set(h2,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Semitone Diffefence')
	ylabel('Reaction Time')
	hold on
	% Plot false alarm rate.
	figure
	Xticks = 1:numDf;
	XtickLabels = {'.25','.5','.75','1','1.5','2','2.5','3','5','7','9','10','12'};
	h3 = plot(Stat(:,:,4),'ro-');
	xlim([1 numDf])
	ylim([0 1])
	set(h3,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Semitone Diffefence')
	ylabel('False Alarm Rate')
	hold on
elseif strfind(key_common,'_loc_')
	% Plot percent success rate vs angular separation.
	figure
	Xticks = 1:numAng;
	XtickLabels = angInd;
	h1 = plot(Stat(4,:,1),'o-');
	xlim([1 numAng])
	ylim([0 1])
	set(h1,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Angular Separation')
	ylabel('Success Rate')
	hold on
	% Plot reaction time vs angular separation.
	figure
	Xticks = 1:numAng;
	XtickLabels = angInd;
	h2 = errorbar(Stat(4,:,2),Stat(4,:,3)','o-');
	xlim([1 numAng])
	ylim([0 max(max(Stat(:,:,2)))+max(max(Stat(:,:,3)))*1.1])
	set(h2,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Angular Separation')
	ylabel('Reaction Time')
	hold on
	% Plot false alarm rate.
	figure
	Xticks = 1:numAng;
	XtickLabels = angInd;
	h3 = plot(Stat(4,:,4),'ro-');
	xlim([1 numAng])
	ylim([0 1])
	set(h3,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Angular Separation')
	ylabel('False Alarm Rate')
	hold on
else
	display('Cannot decide what to plot. More information needed. Sorry!')
end

