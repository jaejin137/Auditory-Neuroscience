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

% Assign top 20 percentil channels and accuracies.
top20 = accumatrix_valid(find(accumatrix_valid(:,2)>=acc_top20),1);
top20(:,2) = accumatrix_valid(find(accumatrix_valid(:,2)>=acc_top20),2);

% Make maps of channels and accuracies and plot them
top20_map_chnum = nan(16,16);
for i=1:length(top20)
	top20_map_chnum(find(LUT==top20(i,1))) = top20(i,1);
end
figure
imagesc(top20_map_chnum)
grid on
set(gca,'GridColor','w')
title('Top 20 Channels')
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),440,400]);

top20_map_accu = nan(16,16);
for i=1:length(top20)
	top20_map_accu(find(LUT==top20(i,1))) = top20(i,2);
end
figure
imagesc(top20_map_accu)
grid on
set(gca,'GridColor','w')
title('Top 20 Accuracies')
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),440,400]);

