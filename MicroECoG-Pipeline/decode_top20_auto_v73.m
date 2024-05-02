%% This script decodes the behavioral responses (hit, miss, FA, and CR)
%% for the top 20 percentile channels using SVM with bootstrap.
%% External functions/scripts used are
%%
%% decode_single_v73.m
%% find_percentile_v73.m
%% decode_select_v73.m
%% plot_accuracy_top20_v73.m


%% Set up
% Check if continue the previous run
tag_continue = 0;
while tag_continue == 0
	ifContinue = input('Continue the previous run?(y/n) ','s');
	if ifContinue == 'n'
		clear
		close all
		% Get the keyword to filter session data
		decodeWhat = input('Which do you want to decode between (r)esponse and (t)nr? ','s');
		key_common = input('Enter session keyword: ','s');
		key_except = {};
		tag_except = 1;
		ind_except = 1;
		while tag_except==1
		    key_except{ind_except} = input('Any exceptions? (hit Enter if none) ','s');
		    if isempty(key_except{ind_except})
		        key_except = key_except(~cellfun('isempty',key_except));
		        tag_except = 0;
		    else
		        ind_except = ind_except+1;
		    end
		end
		
		% Get the session file list
		cmd_list = sprintf('ls *%s*.mat',key_common);
		[status,list] = system(cmd_list);
		preprocDataName = strsplit(list)';
		preprocDataName = preprocDataName(~cellfun(@isempty,preprocDataName));
		sessKey = cell(length(preprocDataName),1);
		for i = 1:length(preprocDataName)
			splitDataName_tmp = strsplit(preprocDataName{i},'_');
			sessKey{i} = splitDataName_tmp{2};
		end
		
		% Channel mapping
		LUT = [78 250 249 229 230 186 177 165 119 111 100 56 55 39 38 82;
		164 245 241 234 238 180 185 171 117 121 107 48 53 47 42 120;
		174 247 243 232 236 178 175 169 122 113 104 108 59 51 40 112;
		242 244 239 235 176 181 168 172 114 115 106 110 52 49 43 44;
		240 251 231 228 184 187 173 170 116 118 101 102 58 57 41 46;
		70 246 237 233 182 179 183 167 123 103 109 105 50 54 45 204;
		33 158 66 124 31 36 25 37 199 205 197 202 196 206 200 203;
		35 1 7 63 3 61 5 29 207 254 217 252 215 209 211 201;
		62 2 60 6 94 88 98 0 213 140 194 210 198 208 192 212;
		96 90 92 91 95 93 97 193 20 143 139 138 142 132 144 136;
		99 89 65 67 125 71 127 248 27 147 148 151 149 141 145 133;
		64 69 126 160 155 163 157 255 21 75 68 76 72 150 146 154;
		153 161 159 189 135 191 137 162 26 83 87 81 77 79 73 74;
		131 128 129 188 130 190 134 156 16 32 86 30 84 24 85 28;
		216 224 220 226 214 222 218 166 22 11 10 14 8 80 12 34;
		223 219 227 221 225 195 253 152 18 23 19 9 17 13 15 4];
		LUT = LUT + 1;

		errorSession = {};

		sessNum_init = 1;
		tag_continue = 1;

	elseif ifContinue == 'y'
		tag_skip = 0;
		while tag_skip == 0
			ifSkip = input('Skip to the next session data?(y/n) ','s');
			if ifSkip == 'y'
				sessNum_init = sessNum+1;
				tag_skip = 1;
			elseif ifSkip == 'n'
				sessNum_init = sessNum;
				tag_skip = 1;
			else
				return
			end
		end
		tag_continue = 1;
	else
		display('Wrong answer!')
	end
end
	

% Iterate over all session data
for sessNum = sessNum_init:length(preprocDataName)
	%% Load session data file.
	fprintf('\nLoading session data for %s ... ',sessKey{sessNum})
	load(preprocDataName{sessNum})
	fprintf('done!\n')
	
	%% Decode individual channel.
	display('Decoding individual channels')
	decode_single_v73
	
	%% Plot decoding map.
	display('Plotting decoding map')
	figure(1)
	if sessNum == 1
		currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1200,1200]);
	end
	%suptitle('Accuracy Map (Broad Band), Low Imp, Ver.7')
	subplot(6,6,sessNum)	
	if decodeWhat == 'r'
		accumatrix = class_acc_single_resp_broad;
	else
		accumatrix = class_acc_single_tnr_broad;
	end
	accumatrix_new = reshape(accumatrix(LUT,1),16,16);
	h11 = imagesc(accumatrix_new);
	h1_c = colorbar;
	caxis([0 100])
	drawnow
	
	%% Find top 20 percentile and plot accuracies.	
	display('Finding top 20 percentile')
	find_percentile_v73
	figure(2)
	if sessNum == 1
		currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1200,1200]);
	end
	%suptitle('Accuracy - Top and Bottom 20 Percentile')
	subplot(6,6,sessNum)
	h21 = scatter(accumatrix_valid(find(accumatrix_valid(:,2)>=acc_top20),3),accumatrix_valid(find(accumatrix_valid(:,2)>=acc_top20),2),'LineWidth',2);
	ylim([0 100])
	xticklabels('off')
	hold on
	h22 = scatter(accumatrix_valid(find(accumatrix_valid(:,2)<=acc_bot20),3),accumatrix_valid(find(accumatrix_valid(:,2)<=acc_bot20),2),'LineWidth',2);
	ylim([0 100])
	xticklabels('off')
	drawnow

	% Select channels to process
	chproc_select = accumatrix_valid(find(accumatrix_valid(:,2)>=acc_top20),1)'

	%% Calculate decoding accuracies with bootstrap for only top 20 percentile.
	if ~isempty(chproc_select)
		decode_select_v73
		% Save the decoding accuracy in a structure.
		%decAccVar = who('class_acc_select_*');
		%for i=1:length(decAccVar)
		%	DecodingAccuracy.(decAccVar{i}) = eval(decAccVar{i});
		%end
		if decodeWhat == 'r';
			GrandDecodingAccuracyResp(sessNum) = DecodingAccuracyResp;
		else
			GrandDecodingAccuracyTNR(sessNum) = DecodingAccuracyTNR;
		end

		% Save decoding accuracy and all related workspace except chdata into a structure
		fprintf('\nSaving decoding accuracy ... ')
		save(preprocDataName{sessNum})
		fprintf('done\n')

		% Plot decoding accuracies as a function of the number of channels.
		display('Plotting decoding accuracies')
		plot_accuracy_top20_v73
	else
		display('No valid channels to decode found!! Proceeding to next session.')
		errorSession(end+1,1) = {preprocDataName{sessNum}};
		continue
	end


	% Clear variables
	clearvars -except decodeWhat LUT preprocDataName sessKey sessNum GrandDecodingAccuracy*
end

