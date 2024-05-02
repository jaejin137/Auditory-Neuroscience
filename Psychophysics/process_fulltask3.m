%% This script takes a human psychophysics data and makes a plot.

%% Clear
clear
%close all


% Some basic parameters
timeCatch = 4200;
respTimeWindow = 1800;

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
%% column 8: response time.
Summary = cell(length(BehavDat),8);
respTime = nan(length(BehavDat),1);

% parameter indices.
dfInd = {'91','92','93','01','96','02','99','03','04','05','07','09','10','11','12'};
angInd = {'-60','-45','-30','-15','0','15','30','45','60'};
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
	% For semitone differences.
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
	elseif Summary{i,2} == '04'	% 3 semitone
		l = 9;
	elseif Summary{i,2} == '05'	% 5 semitone
		l = 10;
	elseif Summary{i,2} == '07'	% 7 semitone
		l = 11;
	elseif Summary{i,2} == '09'	% 9 semitone
		l = 12;
	elseif Summary{i,2} == '10' % 10 semitone
		l = 13;
	elseif Summary{i,2} == '11' % 10 semitone
		l = 14;
	elseif Summary{i,2} == '12'	% 12 semitone
		l = 15;
	else
		fprintf('Error found on Df values: %s\n', Summary{i,2})
	end
	% For azimuthal anglular separations.
%	m = 1+abs(str2num(Summary{i,4})-str2num(Summary{i,5}));
	m = 5+str2num(Summary{i,5})-str2num(Summary{i,4});
	% Count the number of occurances for each cases and assign some basic values.
	respTime(i,1) = BehavDat(i).JoystickAcquired(2) - BehavDat(i).TimeStamp(1);
	if str2num(Summary{i,3}) < timeCatch
		if respTime(i,1) > 0  && respTime(i,1) <= respTimeWindow
			Summary{i,7} = 'h';
			Summary{i,8} = respTime(i,1);
			Hit(l,m) = Hit(l,m)+1;
		elseif BehavDat(i).JoystickAcquired(2) == -1 || respTime(i,1) > respTimeWindow
			Summary{i,7} = 'm';
			Mis(l,m) = Mis(l,m)+1;
		elseif respTime(i,1) <= 0
			Summary{i,7} = 'f';
			FA(l,m) = FA(l,m)+1;
		else
			fprintf('\nError: Cannot determine joystick movement for a target trial!\n')
		end
	elseif str2num(Summary{i,3}) >= timeCatch
		if BehavDat(i).EndOfBehavior ~= -1 %&& BehavDat(i).JoystickAcquired(2) <= BehavDat(i).EndOfBehavior
			Summary{i,7} = 'r';
			Rej(l,m) = Rej(l,m)+1;
		elseif BehavDat(i).EndOfBehavior == -1 %% respTime(i,1) <= 0
			Summary{i,7} = 'f';
			FA(l,m) = FA(l,m)+1;
		else
			fprintf('\nError: Cannot determine joystick movement for a catch trial!\n') 
		end
	else
		fprintf('\nError: Cannot determine if target or catch trial!\n')
	end
end


%% Basic statistics for each parameters.
% Third dimension: success rate, mean response time, s.d. of response time, ratio of false alarm to catch trials,  
Stat = nan(numDf,numAng,4);
for j=1:numAng
	for i=1:numDf
		% Percent success rate
		if Hit(i,j)+Rej(i,j)+Mis(i,j)+FA(i,j) ~= 0
			Stat(i,j,1) = (Hit(i,j)+Rej(i,j))/(Hit(i,j)+Rej(i,j)+Mis(i,j)+FA(i,j));
		end
		% Average of response time
		if Hit(i,j)+Rej(i,j)+Mis(i,j)+FA(i,j) ~= 0 && Hit(i,j) ~= 0
			Stat(i,j,2) = nanmean(cell2mat(Summary(logical((~cellfun(@isempty,strfind(Summary(:,7),'h')).*~cellfun(@isempty,strfind(Summary(:,2),dfInd{i}))).*(5+cell2mat(Summary(:,5))-cell2mat(Summary(:,4)) == j)),8)));
			Stat(i,j,3) = nanstd(cell2mat(Summary(logical((~cellfun(@isempty,strfind(Summary(:,7),'h')).*~cellfun(@isempty,strfind(Summary(:,2),dfInd{i}))).*(5+cell2mat(Summary(:,5))-cell2mat(Summary(:,4)) == j)),8)));
			%Stat(i,j,2) = mean(cell2mat(Summary(logical((~cellfun(@isempty,strfind(Summary(:,7),'h')).*~cellfun(@isempty,strfind(Summary(:,2),dfInd{i}))).*(abs(cell2mat(Summary(:,4))-cell2mat(Summary(:,5)))+1 == j)),8)));
			%Stat(i,j,3) = std(cell2mat(Summary(logical((~cellfun(@isempty,strfind(Summary(:,7),'h')).*~cellfun(@isempty,strfind(Summary(:,2),dfInd{i}))).*(abs(cell2mat(Summary(:,4))-cell2mat(Summary(:,5)))+1 == j)),8)));
		end
		% Ratio of # of false alarms to # of catch trials
		if Rej(i,j)+FA(i,j) ~= 0
			Stat(i,j,4) = FA(i,j)/(Rej(i,j)+FA(i,j));
		end
	end
end

%% Plot the result
if strfind(key_common,'_ton')
	% Plot percent success rate vs semitone differences.
	figure(81)
	%clf
	Xticks = 1:numDf;
	XtickLabels = {'.25','.5','.75','1','1.5','2','2.5','3','4','5','7','9','10','11','12'};
	h1 = plot(Stat(:,:,1),'o-');
	xlim([1 numDf])
	ylim([0 1])
	set(h1,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Semitone Diffefence')
	ylabel('Success Rate')
	hold on
	% Plot response time vs semitone differences.
	figure(82)
	%clf
	Xticks = 1:numDf;
	h2 = errorbar(Stat(:,:,2),Stat(:,:,3),'o-');
	xlim([1 numDf])
	ylim([0 1500])
	%ylim([0 max(max(Stat(:,:,2)))+max(max(Stat(:,:,3)))*1.1])
	set(h2,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Semitone Diffefence')
	ylabel('Response Time')
	hold on
	% Plot false alarm rate.
	figure(83)
	%clf
	Xticks = 1:numDf;
	h3 = plot(Stat(:,:,4),'ro-');
	xlim([1 numDf])
	ylim([0 1])
	set(h3,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Semitone Diffefence')
	ylabel('False Alarm Rate')
	hold on
	% Plot everything
	figure(89)
	%clf
	Xticks = 1:numDf;
	XtickLabels = {'.25','.5','.75','1','1.5','2','2.5','3','4','5','7','9','10','11','12'};
	Yticks = 1:numAng;
	YtickLabels = {'-60','-45','-30','-15','0','15','30','45','60'};
	%[xn, yn] = meshgrid(1:numAng,1:numDf);
	%h99 = mesh(xn,yn,Stat(:,:,1));
	h89 = waterfall(Stat(:,:,1)');
	hidden off;
	xlim([1 numDf]);
	ylim([1 numAng]);
	zlim([0 1]);
	xlabel('Semitone Difference')
	ylabel('Azimuthal Angular Separation')
	zlabel('Success Rate')
	set(h89,'marker','o','linewidth',2,'edgecolor','interp')
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	set(gca,'ytick',Yticks,'yticklabel',YtickLabels)
	hold on
elseif strfind(key_common,'_loc')
	% Plot percent success rate vs angular separation.
	figure(91)
	%clf
	Xticks = 1:numAng;
	XtickLabels = {'-60','-45','-30','-15','0','15','30','45','60'};
	h91 = plot(Stat(:,:,1)','o-');
	xlim([1 numAng])
	ylim([0 1])
	set(h91,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Angular Separation')
	ylabel('Success Rate')
	hold on
	% Plot response time vs angular separation.
	figure(92)
	%clf
	Xticks = 1:numAng;
	h92 = errorbar(Stat(:,:,2)',Stat(:,:,3)','o-');
	xlim([1 numAng])
	ylim([0 1500])
%	ylim([0 max(max(Stat(:,:,2)))+max(max(Stat(:,:,3)))*1.1])
	set(h92,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Angular Separation')
	ylabel('Response Time')
	hold on
	% Plot false alarm rate.
	figure(93)
	%clf
	Xticks = 1:numAng;
	h93 = plot(Stat(:,:,4)','ro-');
	xlim([1 numAng])
	ylim([0 1])
	set(h93,'linewidth',2)
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	xlabel('Angular Separation')
	ylabel('False Alarm Rate')
	hold on
	% Plot everything
	figure(99)
	%clf
	Yticks = 1:numDf;
	YtickLabels = {'.25','.5','.75','1','1.5','2','2.5','3','4','5','7','9','10','11','12'};
	Xticks = 1:numAng;
	XtickLabels = {'-60','-45','-30','-15','0','15','30','45','60'};
	%[xn, yn] = meshgrid(1:numAng,1:numDf);
	%h99 = mesh(xn,yn,Stat(:,:,1));
	h99 = waterfall(Stat(:,:,1));
	hidden off;
	ylim([1 numDf]);
	xlim([1 numAng]);
	zlim([0 1]);
	ylabel('Semitone Difference')
	xlabel('Azimuthal Angular Separation')
	zlabel('Success Rate')
	set(h99,'marker','o','linewidth',2,'edgecolor','interp')
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	set(gca,'ytick',Yticks,'yticklabel',YtickLabels)
	hold on
else
	% Plot everything
	figure(9)
	%clf
	Xticks = 1:numDf;
	XtickLabels = {'.25','.5','.75','1','1.5','2','2.5','3','4','5','7','9','10','11','12'};
	Yticks = 1:numAng;
	YtickLabels = {'-60','-45','-30','-15','0','15','30','45','60'};
	%[xn, yn] = meshgrid(1:numAng,1:numDf);
	%h99 = mesh(xn,yn,Stat(:,:,1));
	h9 = waterfall(Stat(:,:,1));
	hidden off;
	xlim([1 numDf]);
	ylim([1 numAng]);
	zlim([0 1]);
	xlabel('Semitone Difference')
	ylabel('Azimuthal Angular Separation')
	zlabel('Success Rate')
	set(h9,'marker','o','linewidth',2,'edgecolor','interp')
	set(gca,'xtick',Xticks,'xticklabel',XtickLabels)
	set(gca,'ytick',Yticks,'yticklabel',YtickLabels)
	hold on
	display('Warning: Cannot decide what to plot. More information needed. Sorry!')
end

