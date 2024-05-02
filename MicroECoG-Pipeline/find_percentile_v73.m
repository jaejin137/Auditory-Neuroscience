%accumatrix = accCV_broad;
if decodeWhat=='r'
	accumatrix = class_acc_single_resp_broad;
else
	accumatrix = class_acc_single_resp_tnr;
end
accumatrix(:,2) = accumatrix;
accumatrix(:,1) = 1:length(accumatrix);
accumatrix_valid = accumatrix(find(~isnan(accumatrix(:,2))),:);
accumatrix_valid(:,3) = 1;
acc_top20 = prctile(accumatrix_valid(:,2),80);
acc_bot20 = prctile(accumatrix_valid(:,2),20);
%figure
%scatter(accumatrix_valid(find(accumatrix_valid(:,2)>acc_top20),3),accumatrix_valid(find(accumatrix_valid(:,2)>acc_top20),2),'LineWidth',2)
%ylim([0 100])
%xticklabels('off')
%hold on
%scatter(accumatrix_valid(find(accumatrix_valid(:,2)<=acc_bot20),3),accumatrix_valid(find(accumatrix_valid(:,2)<=acc_bot20),2),'LineWidth',2)
%ylim([0 100])
%xticklabels('off')
%title('Accuracy - Top and Bottom 20 Percentile')
%set(gcf,'Position',[currPos(1),currPos(2),480,400])
