respType = {'Hit','Miss','FA','CR');
for m = 1:sum(~cellfun(@isempty,trind))
	if numel(trind{m,1}) >= 16
		figure
		for i=1:16
			subplot(4,4,i)
			spectrogram(chdata_erp_broad_dns{trind{m,1}(i),1},100,20,[],4e3,'yaxis')
			view(150,45)
			colorbar('off')
			drawnow
		end
		suptitle(sprintf('%s',respType{m}))
	end
end
