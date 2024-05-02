load abu_2017-11-03_16-51-51_multiband_roi100.mat
datPathImp
OpenE_Ztest(datPathImp,1,impThrMax);
%-- 10/15/18, 4:23 PM --%
load abu_2017-11-03_16-51-51_multiband_roi100.mat
LUT
figure
LUT(1)
HP
Mag(LUT(1))
OpenE_Ztest_Quick
./
OpenE_Ztest_Quick
.
imagesc(reshape(LUT,16,16))
LUT(1)
LUT(2)
LUT
LUTMod = LUT;
figure
LUTSQ = reshape(LUTMod,16,16);
LUTSQ
imagesc(reshape(LUTMod,16,16))
figure
LUTMod
imagesc(reshape(LUT,16,16))
LUTMod = LUT;
imagesc(reshape(LUTMod,16,16))
figure
imagesc(reshape(LUTMod,16,16))
LUTMod(1:128) = 0;
imagesc(reshape(LUTMod,16,16))
LUTMod(1:128) = 1000;
imagesc(reshape(LUTMod,16,16))
colormap(jet)
LUTMod(1:128) = 100;
imagesc(reshape(LUTMod,16,16))
LUTMod(1:128) = 1000;
imagesc(reshape(LUTMod,16,16))
LUTMod = LUT;
imagesc(reshape(LUTMod,16,16))
LUTMod(1:16) = 1000;
imagesc(reshape(LUTMod,16,16))
LUTMod = LUT; i=1; j=4; LUTMod(16*(i-1)+j) = 1000; imagesc(reshape(LUTMod,16,16))
LUTMod = LUT; i=4; j=4; LUTMod(16*(i-1)+j) = 1000; imagesc(reshape(LUTMod,16,16))
size(LUTMod)
size(LUT)
LUTMod(4,4) = 1000; LUTMod(4,20) = 1000; LUTMod(20,4) = 1000; LUTMod(20,20) = 1000;
imagesc(LUTMod)
imagesc(reshape(LUTMod,16,16))
size(LUTMod)
LUTMod = LUT;
LUTMod(4,4) = 1000; LUTMod(4,12) = 1000; LUTMod(12,4) = 1000; LUTMod(12,12) = 1000;
imagesc(LUTMod)
LUTMod = LUT;
LUT(2,2)
LUT(2,6)
LUT(12,4)
LUT(12,12)
chproc
chproc = [246,181,161,77];
decode_select_v73
chproc_select
chproc_select = [246,181,161,77];
decode_select_v73
PN = 0:log2(256)
numch = length(chproc_select);
[2.^PN(2.^PN<numch) numch]
decode_givench_v73
r
decode_givench_v73
t
chproc
decode_single_v73
tnrLabel
load abu_2017-11-28_15-26-26_multiband_roi100.mat
chproc
chproc = [246,181,161,77];
chproc_select = [246,181,161,77];
condIndUB
decode_single_v73
classInd
i
~cellfun(@isempty,classInd(:,i))
classInd(i,~cellfun(@isempty,classInd(:,i)))
decodeWhat = 't';
decodeWhat = 't';
decode_single_v73
minoccur
decodeWhat
decode_single_v73
decodeWhat
chproc
class_acc_single_tnr_broad([77,161,181,246])
figure
imagesc(reshape(class_acc_single_tnr_broad,16,16))
size(reshape(class_acc_single_tnr_broad,16,16))
size(reshape(class_acc_single_tnr_broad(LUT),16,16))
imagesc(reshape(class_acc_single_tnr_broad(LUT),16,16))
class_acc_single_tnr_broad(find(~nan(class_acc_single_tnr_broad)))
