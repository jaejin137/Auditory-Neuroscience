% Create reference area map
htbl_area = ones(15,15)*-1;
if strcmp(monkeyName{m},'Cassius')
	coord_core = [coord_ac(2:23,1); coord_ac(2:19,2)];	
	coord_belt = [coord_ac(2:36,3); coord_ac(2:24,4)];
else
	coord_core = [coord_ac(2:20,1); coord_ac(2:25,2)];	
	coord_belt = [coord_ac(2:31,3); coord_ac(2:16,4)];
end

for ldx_core = 1:numel(coord_core)
	htbl_area(-coord_core{ldx_core}(2)+8,coord_core{ldx_core}(1)+8) = 1;
end
for ldx_belt = 1:numel(coord_belt)
	htbl_area(-coord_belt{ldx_belt}(2)+8,coord_belt{ldx_belt}(1)+8) = 0;
end

cmap = [0 0 0; 0 1 0; 1 0 0];
h_refmap(m) = figure;
heatmap(htbl_area,'CellLabelColor','none')
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
currPos = get(h_refmap(m),'Position'); set(h_refmap(m),'Position',[currPos(1),currPos(2),394,344]);
saveas(h_refmap(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_RefMap_v5_%s.png',monkeyName{m}))
