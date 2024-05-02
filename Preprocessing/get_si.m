%% Get SI(similarity indices) using peak amplitudes for a given df using moving time window.
TMWNDWTH = 800; % width of time window(in msec) for similarity index.
TMWNDMOV = 25; % size of moving of time window(in msec) for similarity index.
fprintf('\nCalculating similarity indices ...\n')
% Convert time to the number of data points.
tw_width = round(TMWNDWTH/1000*fs_neural);
tw_mov = round(TMWNDMOV/1000*fs_neural);
num_bin = floor(pstpmin/tw_mov)+1;
num_df = length(DF);
peak_int_p5 = zeros(num_int_p5,num_bin);
peak_int_3 = zeros(num_int_3,num_bin);
peak_int_5 = zeros(num_int_5,num_bin);
peak_int_12 = zeros(num_int_12,num_bin);
peak_seg_p5 = zeros(num_seg_p5,num_bin);
peak_seg_3 = zeros(num_seg_3,num_bin);
peak_seg_5 = zeros(num_seg_5,num_bin);
peak_seg_12 = zeros(num_seg_12,num_bin);
peak_int_mean = zeros(num_df,num_bin);
peak_seg_mean = zeros(num_df,num_bin);
si_peak = zeros(num_df,num_bin);
for j=1:num_bin
	for i=1:num_int_p5
		peak_int_p5(i,j) = max(lfp_int_p5(i,(j-1)*tw_mov+1:(j-1)*tw_mov+1+tw_width));
	end
	for i=1:num_int_3
		peak_int_3(i,j) = max(lfp_int_3(i,(j-1)*tw_mov+1:(j-1)*tw_mov+1+tw_width));
	end
	for i=1:num_int_5
		peak_int_5(i,j) = max(lfp_int_5(i,(j-1)*tw_mov+1:(j-1)*tw_mov+1+tw_width));
	end
	for i=1:num_int_12
		peak_int_12(i,j) = max(lfp_int_12(i,(j-1)*tw_mov+1:(j-1)*tw_mov+1+tw_width));
	end
	for i=1:num_seg_p5
		peak_seg_p5(i,j) = max(lfp_seg_p5(i,(j-1)*tw_mov+1:(j-1)*tw_mov+1+tw_width));
	end
	for i=1:num_seg_3
		peak_seg_3(i,j) = max(lfp_seg_3(i,(j-1)*tw_mov+1:(j-1)*tw_mov+1+tw_width));
	end
	for i=1:num_seg_5
		peak_seg_5(i,j) = max(lfp_seg_5(i,(j-1)*tw_mov+1:(j-1)*tw_mov+1+tw_width));
	end
	for i=1:num_seg_12
		peak_seg_12(i,j) = max(lfp_seg_12(i,(j-1)*tw_mov+1:(j-1)*tw_mov+1+tw_width));
	end
	peak_int_mean(1,j) = sum(peak_int_p5(:,j),1)/num_int_p5;
	peak_int_mean(2,j) = sum(peak_int_3(:,j),1)/num_int_3;
	peak_int_mean(3,j) = sum(peak_int_5(:,j),1)/num_int_5;
	peak_int_mean(4,j) = sum(peak_int_12(:,j),1)/num_int_12;
	peak_seg_mean(1,j) = sum(peak_seg_p5(:,j),1)/num_seg_p5;
	peak_seg_mean(2,j) = sum(peak_seg_3(:,j),1)/num_seg_3;
	peak_seg_mean(3,j) = sum(peak_seg_5(:,j),1)/num_seg_5;
	peak_seg_mean(4,j) = sum(peak_seg_12(:,j),1)/num_seg_12;
	for k=1:num_df
		si_peak(k,j) = (peak_int_mean(k,j)-peak_seg_mean(k,j))/(peak_int_mean(k,j)+peak_seg_mean(k,j));
	end
end

% Mean and s.d. of SI over all time bins
if nnz(isnan(si_peak))
	cut_nan = floor(min(find(isnan(si_peak)==1)+1)/4);
	si_peak_mean = mean(si_peak(:,1:cut_nan-1),2);
	si_peak_std = std(si_peak(:,1:cut_nan-1),1,2);
else
	si_peak_mean = mean(si_peak,2);
	si_peak_std = std(si_peak,1,2);
end


%% Get ITPC(inter-trial phase coherence) using phase.
phase_int_p5 = zeros(num_int_p5,tw_width,num_bin);
phase_int_3 = zeros(num_int_3,tw_width,num_bin);
phase_int_5 = zeros(num_int_5,tw_width,num_bin);
phase_int_12 = zeros(num_int_12,tw_width,num_bin);
phase_seg_p5 = zeros(num_seg_p5,tw_width,num_bin);
phase_seg_3 = zeros(num_seg_3,tw_width,num_bin);
phase_seg_5 = zeros(num_seg_5,tw_width,num_bin);
phase_seg_12 = zeros(num_seg_12,tw_width,num_bin);
itpc_int_mean = zeros(num_df,num_bin);
itpc_seg_mean = zeros(num_df,num_bin);
si_itpc = zeros(num_df,num_bin);
for j=1:num_bin
	% For .5 semitone intergration trials.
	phase_int_p5(:,:,j) = atan(imag(lfp_ht_int_p5(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width))./real(lfp_ht_int_p5(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width)));
	sum_cos_int_p5 = sum(cos(phase_int_p5(:,:,j)),1);
	sum_sin_int_p5 = sum(sin(phase_int_p5(:,:,j)),1);
	itpc_int_p5 = sqrt(sum_cos_int_p5.^2+sum_sin_int_p5.^2)/num_int_p5;
	itpc_int_mean(1,j) = sum(itpc_int_p5,2)/tw_width;
	% For 3 semitone intergration trials.
	phase_int_3(:,:,j) = atan(imag(lfp_ht_int_3(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width))./real(lfp_ht_int_3(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width)));
	sum_cos_int_3 = sum(cos(phase_int_3(:,:,j)),1);
	sum_sin_int_3 = sum(sin(phase_int_3(:,:,j)),1);
	itpc_int_3 = sqrt(sum_cos_int_3.^2+sum_sin_int_3.^2)/num_int_3;
	itpc_int_mean(2,j) = sum(itpc_int_3,2)/tw_width;
	% For 5 semitone intergration trials.
	phase_int_5(:,:,j) = atan(imag(lfp_ht_int_5(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width))./real(lfp_ht_int_5(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width)));
	sum_cos_int_5 = sum(cos(phase_int_5(:,:,j)),1);
	sum_sin_int_5 = sum(sin(phase_int_5(:,:,j)),1);
	itpc_int_5 = sqrt(sum_cos_int_5.^2+sum_sin_int_5.^2)/num_int_5;
	itpc_int_mean(3,j) = sum(itpc_int_5,2)/tw_width;
	% For 12 semitone intergration trials.
	phase_int_12(:,:,j) = atan(imag(lfp_ht_int_12(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width))./real(lfp_ht_int_12(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width)));
	sum_cos_int_12 = sum(cos(phase_int_12(:,:,j)),1);
	sum_sin_int_12 = sum(sin(phase_int_12(:,:,j)),1);
	itpc_int_12 = sqrt(sum_cos_int_12.^2+sum_sin_int_12.^2)/num_int_12;
	itpc_int_mean(4,j) = sum(itpc_int_12,2)/tw_width;
	% For .5 semitone segregation trials.
	phase_seg_p5(:,:,j) = atan(imag(lfp_ht_seg_p5(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width))./real(lfp_ht_seg_p5(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width)));
	sum_cos_seg_p5 = sum(cos(phase_seg_p5(:,:,j)),1);
	sum_sin_seg_p5 = sum(sin(phase_seg_p5(:,:,j)),1);
	itpc_seg_p5 = sqrt(sum_cos_seg_p5.^2+sum_sin_seg_p5.^2)/num_seg_p5;
	itpc_seg_mean(1,j) = sum(itpc_seg_p5,2)/tw_width;
	% For 3 semitone segregation trials.
	phase_seg_3(:,:,j) = atan(imag(lfp_ht_seg_3(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width))./real(lfp_ht_seg_3(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width)));
	sum_cos_seg_3 = sum(cos(phase_seg_3(:,:,j)),1);
	sum_sin_seg_3 = sum(sin(phase_seg_3(:,:,j)),1);
	itpc_seg_3 = sqrt(sum_cos_seg_3.^2+sum_sin_seg_3.^2)/num_seg_3;
	itpc_seg_mean(2,j) = sum(itpc_seg_3,2)/tw_width;
	% For 5 semitone segregation trials.
	phase_seg_5(:,:,j) = atan(imag(lfp_ht_seg_5(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width))./real(lfp_ht_seg_5(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width)));
	sum_cos_seg_5 = sum(cos(phase_seg_5(:,:,j)),1);
	sum_sin_seg_5 = sum(sin(phase_seg_5(:,:,j)),1);
	itpc_seg_5 = sqrt(sum_cos_seg_5.^2+sum_sin_seg_5.^2)/num_seg_5;
	itpc_seg_mean(3,j) = sum(itpc_seg_5,2)/tw_width;
	% For 12 semitone segregation trials.
	phase_seg_12(:,:,j) = atan(imag(lfp_ht_seg_12(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width))./real(lfp_ht_seg_12(:,(j-1)*tw_mov+1:(j-1)*tw_mov+tw_width)));
	sum_cos_seg_12 = sum(cos(phase_seg_12(:,:,j)),1);
	sum_sin_seg_12 = sum(sin(phase_seg_12(:,:,j)),1);
	itpc_seg_12 = sqrt(sum_cos_seg_12.^2+sum_sin_seg_12.^2)/num_seg_12;
	itpc_seg_mean(4,j) = sum(itpc_seg_12,2)/tw_width;
	% Get SI using temporal means of itpc for integration and segretation.
	for k=1:4
		si_itpc(k,j) = itpc_int_mean(k,j)-itpc_seg_mean(k,j);
	end
end

% Mean and s.d. of SI over all time bins
if nnz(isnan(si_itpc))
	cut_nan = floor(min(find(isnan(si_itpc)==1)+1)/4);
	si_itpc_mean = mean(si_itpc(:,1:cut_nan-1),2);
	si_itpc_std = std(si_itpc(:,1:cut_nan-1),1,2);
else
	si_itpc_mean = mean(si_itpc,2);
	si_itpc_std = std(si_itpc,1,2);
end

fprintf('Done!\n')
