%function [Good_Z_loc] = OpenE_Ztest(foldername,plotswitch)

% ARO Z test 

%Init
%clear all;
%close all;
%clc;

%% User Input
    %folder_path = 'C:\Users\ken_c\Dropbox\ARO data\2017-06-08_16-11-34_013';
    %folder_path = 'C:\Users\ken_c\Dropbox\ARO data\2017-06-10_16-28-44_001';
    %folder_path = '2017-06-22_15-38-14';
    %folder_path = sprintf('/home/jaejin/ArmyProject/open-ephys-gui-master/%s',foldername);
    %filename = 'impedance_measurement';
	%foldername = input('Enter the impdedance folder name: ','s');
	foldername = '.';
	filename = input('Enter the impedance file name: ','s');
    Ch_count = 256;
    EL_design = 5; % 1=Center, 2=Left, 3=Right, 4=Human Project, 5=ARO final mapping
    Thr_range = [1000 1E6]; % [Lower_limit Higher_limit], unit in ohm, -1 = auto
    Red = 1E6;

plotswitch = 1;

%% Getting Z data from proper ZIF connection
    %OpenE_file = xml2struct(strcat(folder_path,'/',filename,'.xml'));
    OpenE_file = xml2struct(strcat(foldername,'/',filename,'.xml'));

    for i = 1:Ch_count
        %[Ch_num{i} Mag{i} Name{i} Phase{i} Stream{i}] = OpenE_file.Children(i*2).Attributes.Value;
        [a b c d e] = OpenE_file.Children(i*2).Attributes.Value;
        Ch_num(i) = str2num(a);
        Mag(i) = str2num(b);
        Name{i} = c;
        Phase(i) = str2num(d);
        Stream(i) = str2num(e);
        clear a b c d e;
    end

    Z = @(a) (Mag(a)+exp(1i*Phase(a)));

    
    lmt = (1:10)*1E5;
    if Thr_range(2) == -1
        temp = find(max(Mag<Red));
        %Thr_range(2) = Mag(temp);
        
        tmp = abs(lmt-Mag(temp));
        [idx idx] = min(tmp); %index of closest value
        closest = lmt(idx); %closest value
        if Mag(temp) > closest
            closest = lmt(idx+1);            
        end
        Thr_range(2) = closest;
    end
    
    
    
    
    %valid_ch = [16:25 4:15 79:90 67:78 142:149 151:153 130:141 205:216 193 194 196:204];
    %valid_ch = sort(valid_ch);
    valid_ch = 1:256;

    % Center EL, v1-1
        Ignore_ch{1} = []; % All working

    % Left EL, v1-2
        Ignore_ch{2} = [16 25 5 87 88 72 146 139 140 194 196 198 200 199 197];
        Ignore_ch{2} = sort(Ignore_ch{2});

    % Right EL, v1-3 
        Ignore_ch{3} = [19:23 25 79 82 73 153 135 132 211 203];
        Ignore_ch{3} = sort(Ignore_ch{3});
    % Human Project
        Ignore_ch{4} = [1 31 63 65 95 127 129 159 191 193 223 255]+1;
        Ignore_ch{4} = sort(Ignore_ch{4});
    % ARO Final
        Ignore_ch{5} = [78 82];  
        Ignore_ch{5} = sort(Ignore_ch{5});

        Final_ch_map = setxor(valid_ch,Ignore_ch{EL_design});    
        Res = real(Z(valid_ch)); %unit in ohm
    
    %For Human projecto only
%     LUT = [1 31 49 39 61 44 58 32 6 20 5 23 27 29 63 65;% ...
%                    95 33 35 53 47 50 36 62 24 10 19 9 13 15 225 127;% ...
%                    93 79 51 37 48 46 60 34 8 22 7 21 25 243 227 241;% ...
%                    91 77 89 55 45 52 38 0 26 12 17 11 247 229 245 231;% ...
%                    87 73 85 75 59 42 56 30 4 18 3 251 237 240 239 253;% ...
%                    69 83 71 81 67 54 40 2 28 14 246 234 244 238 242 236;% ...
%                    84 74 86 76 82 78 41 57 16 233 232 248 230 252 228 250;% ...
%                    70 88 72 90 68 92 80 43 235 249 194 222 192 226 254 224;% ...
%                    96 126 98 64 94 66 121 107 171 208 220 196 218 200 216 198;% ...
%                    122 100 124 102 120 104 105 144 185 169 206 210 204 214 202 212;% ...
%                    108 114 110 116 106 118 142 156 130 168 182 195 209 199 211 197;% ...
%                    125 111 112 109 123 131 146 132 158 184 170 187 203 213 201 215;% ...
%                    103 117 101 119 139 145 140 154 128 166 180 173 183 217 205 219;% ...
%                    113 99 115 153 149 135 150 136 162 188 174 176 165 179 207 221;% ...
%                    129 97 143 141 137 147 138 152 190 164 178 175 181 163 161 159;% ...
%                    191 193 157 155 151 133 148 134 160 186 172 189 167 177 223 255];
   
    %For ARO complete set only            
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

    LUT = LUT + 1; %change to 1 base instead of 0 base
    
 %% Calculating the yield
    EL_ch_count = [92 81 81 244 254];
    Good_Z_loc = find((Res <= Red) & (Res >= Thr_range(1)));
    Good_Z = Mag(Good_Z_loc); %Real only
    Yield = length(Good_Z_loc)/EL_ch_count(EL_design)*100;
    MeanZ = mean(Good_Z/1E3);
    stdev = std(Good_Z/1E3);
    
    %Re-map
    for i = 1:256
        HP(i) = Mag(LUT(i));
    end
    HP_grid = reshape(HP,16,16);

	% Plot only if plotswitch is on
	if plotswitch == 1;
		figure
    	imagesc(HP_grid)
    	%title('Z result - Human Project (face up)');
    	%title('Z result - ARO EL (face up)')
    	h = colorbar; 
    	caxis([0 Thr_range(2)])

    	set(get(h,'title'), 'string', 'Ohms')

    	% set(h,'scale', 'log')
    	ylabel('Rows')
    	xlabel('Columns')
    	
    	s = [num2str(EL_ch_count(EL_design)),'ch (Face-up): ', ' Mean  =  ', num2str(MeanZ),  '  �  ', num2str(stdev), '   kOhms      Yield  = ', num2str(Yield),' %'];
    	%title(s,'fontsize',8);
		currPos = get(gcf,'Position');
		set(gcf,'Position',[currPos(1),currPos(2),400,340]);
		set(gca,'FontSize',16);
	end

