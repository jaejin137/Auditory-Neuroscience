%%================== Plotting for parprocess_abu_chorus_v7 ====================

plotImpMap = input('Also plot the impedance map?(y/n) ','s');
display('Generating figures... Please wait.')

%% Plot impedance measure.
if plotImpMap == 'y'
	OpenE_Ztest(datPathImp,1,impThrMax);
end


%------------------------------------------------------------------------------
% Plot % correct vs numch.
figure
sigBand = {'broad','theta','alpha','beta','gamma','hfo'};
for m = 1:numel(chstat_feed)
	if bandProc == 'm'
		subplot(2,3,m)
	end
	errorbar(classAcc{m,2},classAcc{m,3},'k-','LineWidth',2)
	hold on
	errorbar(classAcc{m,6},classAcc{m,7},'g--','LineWidth',2)
	ylim([0 100])
	xlabel('Number of channels')
	ylabel('% Correct')
	Xticks = 1:length(numchpow);
	XtickLabels = numchpow;
	set(gca,'Xtick',Xticks,'XtickLabel',XtickLabels)
	str1 = sprintf('%s',sigBand{m});
	str2 = sprintf('%s_{bstr}',sigBand{m});
	legend(str1,str2)
	legend boxoff
	% Resize figure
	%currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1400,800]);
end


return

%------------------------------------------------------------------------------
% Plot variance vs mean.
if ifFilter == 'n'
	%for i=1:length(condIndUB)
	i = length(condIndUB);
		if numel(classUniq{i}) >= 2
			figure
			title(sprintf('Condition %i',i))
			for k = 1:length(chproc)
				gscatter(chstat_broad(condIndUB{i},1,chproc(k)),chstat_broad(condIndUB{i},2,chproc(k)),cell2mat(behavSummary(condIndUB{i},classLabel)))
				%gscatter(nanmean(chstat_broad(condIndUB{i},1,:),3),nanmean(chstat_broad(condIndUB{i},2,:),3),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_broad(i).SVs(:,1),svmModel_broad(i).SVs(:,2),'ko')
			%for k = chproc
			%	scatter3(chstat_broad(condIndUB{i},1,k),chstat_broad(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)),str2num(cell2mat(behavSummary(condIndUB{i},4))),cell2mat(behavSummary(condIndUB{i},8)),'filled')
			%	drawnow
			%	hold on
			%end
%			xlim([-3 3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			%zlabel('TNR')
			%ZtickLabels = tnrLabel;
			%set(gca,'ZtickLabel',ZtickLabels)
			title('Broad')
			if isnumeric(classUniq{i})
				legend([className{classUniq{i}}])
			else
				legend(classUniq{i})
			end
			%legend boxoff
		end
	%end

else

	for i=1:length(condIndUB)
		if numel(classUniq{i}) >= 2
			figure
			suptitle(sprintf('Condition %i',i))
			% For broad data
			subplot(2,3,1)
			for k = chproc
				gscatter(chstat_broad(condIndUB{i},1,k),chstat_broad(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_broad(i).SVs(:,1),svmModel_broad(i).SVs(:,2),'ko')
%			xlim([-3 3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Broad')
			legend([className{classUniq{i}}])
			legend boxoff

			% For theta band
			subplot(2,3,2)
			for k = chproc
				gscatter(chstat_theta(condIndUB{i},1,k),chstat_theta(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_theta(i).SVs(:,1),svmModel_theta(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Theta')
			legend([className{classUniq{i}}])
			legend boxoff

			% For alpha band
			subplot(2,3,3)
			for k = chproc
				gscatter(chstat_alpha(condIndUB{i},1,k),chstat_alpha(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_alpha(i).SVs(:,1),svmModel_alpha(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Alpha')
			legend([className{classUniq{i}}])
			legend boxoff

			% For beta band
			subplot(2,3,4)
			for k = chproc
				gscatter(chstat_beta(condIndUB{i},1,k),chstat_beta(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_beta(i).SVs(:,1),svmModel_beta(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Beta')
			legend([className{classUniq{i}}])
			legend boxoff

			% For gamma band
			subplot(2,3,5)
			for k = chproc
				gscatter(chstat_gamma(condIndUB{i},1,k),chstat_gamma(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_gamma(i).SVs(:,1),svmModel_gamma(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Gamma')
			legend([className{classUniq{i}}])
			legend boxoff

			% For hfo band
			subplot(2,3,6)
			for k = chproc
				gscatter(chstat_hfo(condIndUB{i},1,k),chstat_hfo(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_hfo(i).SVs(:,1),svmModel_hfo(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('HFO')
			legend([className{classUniq{i}}])
			legend boxoff
			% Resize figure
			%currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1400,800]);
		end
	end
end


%% Plot mean trial by trial activities over all channel grid
%figure
%subplot(1,2,1)
%h1 = imagesc(chstat_broad_mean_grid(:,:,1));
%set(h1,'Alphadata',~isnan(chstat_broad_mean_grid(:,:,1)))
%caxis([0 1])
%colorbar
%title(sprintf('Mean, Tr#%03i',1),'fontsize',8)
%subplot(1,2,2)
%h2 = imagesc(chstat_broad_var_grid(:,:,2));
%set(h2,'Alphadata',~isnan(chstat_broad_var_grid(:,:,1)))
%caxis([0 1])
%colorbar
%title(sprintf('Var, Tr#%03i',1),'fontsize',8)
%%currPos = get(gcf,'Position');
%set(gcf,'Position',[%currPos(1),currPos(2),700,200]);
%for k=2:2*numtr
%	set(h1,'cdata',chstat_broad_mean_grid(:,:,k));
%	set(h2,'cdata',chstat_broad_var_grid(:,:,k));
%	if mod(k,2)
%		subplot(1,2,1)
%		title(sprintf('Mean, Tr#%i',k),'fontsize',8)
%		subplot(1,2,2)
%		title(sprintf('Var, Tr#%i',k),'fontsize',8)
%	end
%	pause(.5)
%	dbroadnow;
%end


