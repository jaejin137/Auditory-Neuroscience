% Create reference area map
monkeyName{1} = 'Cassius';
monkeyName{2} = 'Miyagi';
for m = 1:2
	htbl_area = ones(15,15)*-1;
	load(sprintf('coord_ac_%s.mat',monkeyName{m}))
	if strcmp(monkeyName{m},'Cassius')
		coord_core = [coord_ac(2:23,1); coord_ac(2:19,2)];	
		coord_belt = [coord_ac(2:36,3); coord_ac(2:24,4)];
	else
		coord_core = [coord_ac(2:20,1); coord_ac(2:25,2)];	
		coord_belt = [coord_ac(2:31,3); coord_ac(2:16,4)];
	end
	
	for ldx_core = 1:numel(coord_core)
		htbl_area(coord_core{ldx_core}(2)+8,coord_core{ldx_core}(1)+8) = 1;
	end
	for ldx_belt = 1:numel(coord_belt)
		htbl_area(coord_belt{ldx_belt}(2)+8,coord_belt{ldx_belt}(1)+8) = 0;
	end
	
	coordLabels_x = num2str([-7:7]');
	coordLabels_y = flipud(coordLabels_x);
	
	cmap = [0 0 0; 0 1 0; 1 0 0];
	h_refmap(m) = figure;
	heatmap(flipud(htbl_area),'CellLabelColor','none')
	grid off
	colormap(cmap);
	colorbar('off')
	set(gca,'XDisplayLabels',coordLabels_x)
	set(gca,'YDisplayLabels',coordLabels_y)
	if strcmp(monkeyName{m},'Cassius')
		xlabel('Medial <----> Lateral')
	else
		xlabel('Lateral <----> Medial')
	end
	ylabel('Posterior <----> Anterior')
	sgtitle(sprintf('Area Map for %s',monkeyName{m}),'FontSize',24)
	currPos = get(h_refmap(m),'Position'); set(h_refmap(m),'Position',[currPos(1),currPos(2),420,360]);
	saveas(h_refmap(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_RefMap_v4_%s.png',monkeyName{m}))
end
