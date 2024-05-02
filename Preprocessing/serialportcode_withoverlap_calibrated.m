



function serialportcode_withoverlap_calibrated

% Modified by KLCL

filePath = 'C:\Users\cohenlab\Documents\MATLAB'

% calibratedvalues.m


% Location of TDT circuit

circuitFile = fullfile('C:\','Users','cohenlab','Documents','MATLAB', 'tonetriplesTDT2.rcx'); 

samplerate=48828.125;
durationofstim=10;
durationofsilence=50;
stimulusnpts=durationofstim*samplerate;
silencenpts=durationofsilence*samplerate;
npts=stimulusnpts+silencenpts; 

% Initial Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Connect to RX6 and load circuit

RP=actxcontrol('RPco.x',[5 5 26 26]);
if RP.ConnectRZ6( 'GB', 1)
    disp('Connected!');
else
    disp('Unable to connect to RZ6');
end

%%%%%End Initial Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















% Initialize Serial Port, make appropriate COM channel
s2 = serial('COM1','BaudRate',9600,'DataBits',7);
set(s2, 'Terminator', 'CR/LF');
fopen(s2);
global s2

endoftrial=0;
while (endoftrial~=1)
  frequencies=[     400
         425
         450
         475
         500
         525
         550
         575
         600
         625
         650
         675
         700
         725
         750
         775
         800
         825
         850
         875
         900
         925
         950
         975
        1000
        1025
        1050
        1075
        1100
        1125
        1150
        1175
        1200
        1225
        1250
        1275
        1300
        1325
        1350
        1375
        1400
        1425
        1450
        1475
        1500
        1525
        1550
        1575
        1600
        1625
        1650
        1675
        1700
        1725
        1750
        1775
        1800
        1825
        1850
        1875
        1900
        1925
        1950
        1975
        2000
        2025
        2050
        2075
        2100
        2125
        2150
        2175
        2200
        2225
        2250
        2275
        2300
        2325
        2350
        2375
        2400
        2425
        2450
        2475
        2500
        2525
        2550
        2575
        2600
        2625
        2650
        2675
        2700
        2725
        2750
        2775
        2800
        2825
        2850
        2875
        2900
        2925
        2950
        2975
        3000
        3025
        3050
        3075
        3100
        3125
        3150
        3175
        3200
        3225
        3250
        3275
        3300
        3325
        3350
        3375
        3400
        3425
        3450
        3475
        3500
        3525
        3550
        3575
        3600
        3625
        3650
        3675
        3700
        3725
        3750
        3775
        3800
        3825
        3850
        3875
        3900
        3925
        3950
        3975
        4000];
    
    
    amplitudes=[ 0.0071
    0.0116
    0.0079
    0.0061
    0.0086
    0.0096
    0.0079
    0.0103
    0.0148
    0.0104
    0.0095
    0.0082
    0.0077
    0.0077
    0.0090
    0.0117
    0.0093
    0.0087
    0.0087
    0.0077
    0.0077
    0.0089
    0.0097
    0.0097
    0.0097
    0.0103
    0.0103
    0.0116
    0.0106
    0.0094
    0.0104
    0.0130
    0.0142
    0.0118
    0.0099
    0.0099
    0.0099
    0.0099
    0.0107
    0.0114
    0.0114
    0.0144
    0.0159
    0.0159
    0.0131
    0.0112
    0.0112
    0.0112
    0.0147
    0.0172
    0.0154
    0.0187
    0.0168
    0.0144
    0.0123
    0.0096
    0.0104
    0.0135
    0.0123
    0.0123
    0.0161
    0.0337
    0.0149
    0.0102
    0.0102
    0.0086
    0.0102
    0.0135
    0.0148
    0.0202
    0.0147
    0.0115
    0.0115
    0.0128
    0.0108
    0.0090
    0.0096
    0.0117
    0.0126
    0.0151
    0.0121
    0.0107
    0.0107
    0.0107
    0.0088
    0.0088
    0.0109
    0.0138
    0.0163
    0.0134
    0.0102
    0.0094
    0.0094
    0.0115
    0.0165
    0.0165
    0.0133
    0.0104
    0.0088
    0.0088
    0.0099
    0.0129
    0.0142
    0.0142
    0.0142
    0.0122
    0.0122
    0.0131
    0.0131
    0.0125
    0.0133
    0.0114
    0.0103
    0.0115
    0.0131
    0.0168
    0.0124
    0.0095
    0.0095
    0.0108
    0.0096
    0.0105
    0.0119
    0.0098
    0.0098
    0.0092
    0.0091
    0.0091
    0.0105
    0.0120
    0.0142
    0.0142
    0.0106
    0.0078
    0.0078
    0.0094
    0.0108
    0.0097
    0.0097
    0.0103
    0.0087
    0.0073
    0.0084
    0.0084
    0.0095];

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     Have to reload and reset circuit each time
%     
    disp('Loading...');
    if RP.LoadCOF(circuitFile)
        disp(strcat('Loaded...',circuitFile));
    else
        disp('FAIL');
    end
    
    
    %%%Check for params from Labview
    checkforparams=1;
    disp(['Waiting for Labview...'])
    while checkforparams==1
        idn = fscanf(s2);
        disp([ ]);disp(idn);
       
        if length(idn)>=6
            try
                
                [frequency, stepsize, rate, overlap]=parse(idn);
                freq1=frequency;
         
                
                if isequal(stepsize, 1)
                    stepsize=.5;
                end
                
                freq2=freq1*(2^(1/12))^(stepsize);
            catch
                disp('Error: Labview output in incorrect format.  Resend data.')
                continue;
            end
            if param_check==1 
                checkforparams=0;
            end
        end
       
    end
  
   
 
   
   
    % Waits here until receives play or abort signal from Labview
    abort=determine_play_trial;
    if abort==1
        disp(['Abort...'])
        continue;
    end
   
    
    %%Run circuit
    if RP.Run
        disp('Running circuit!');
    else
        disp('Error running circuit!');
    end
    
    
    


    
    if all(bitget(RP.GetStatus,1:3))
        
        
        t=(0:(1/samplerate):(.04-(1/samplerate)));  %%%%.04= 40 ms per sound
        
        silence=(1/rate)-.04;
        
        t2=(0:(1/samplerate):(silence-(1/samplerate)));
        
% % %         calib=load ('calibration-20120713T180221.mat');
% % %         listfrequencies =[freq1, freq2];
        amps= interp1(frequencies, amplitudes, [freq1 freq2]);
        
        
        f1=sin(2*pi*t*freq1)*amps(1);
        f2=sin(2*pi*t*freq2)*amps(2);
        f4=f1+f2;
        
        f3=zeros(1, size(t2,2));
        
        if isequal(overlap, 0)
            doubleset=[f1 f3 f2 f3];
        elseif isequal(overlap, 1)
            doubleset=[f4 f3 f4 f3];
        else
            doubleset=[f1 f3 f2 f3];
        end
        
        doublesets=[];
        
        for q=1:20
            doublesets=[doublesets doubleset];
        end
        
       
        RP.SoftTrg(1);
        
        RP.WriteTagVEX('datain', 0, 'F32', doublesets);
    end
    
    
    
    
    
    
%    After circuit starts playing, waits for stop signal
    stop_circuit=determine_stop_trial;
    if stop_circuit==1
        
        RP.Halt;
        RP.ClearCOF;
        continue;
    end
   
   
end




    function [stop_circuit]=determine_stop_trial
        keepchecking2=1;
        while keepchecking2==1
            stop_circuit=[];
           
            % checking for abort signal
            idn = fscanf(s2);
            if strfind(idn,'9')
                stop_circuit=1;
                keepchecking2=0;
                disp(['Abort...'])
            end
           
        end
       
    end



    function abort=determine_play_trial
        keepchecking1=1;
        while keepchecking1==1
            abort=[];
            idn2 = fscanf(s2);
%             disp([ ]);disp(idn);disp([ ])
           
            if strfind(idn2,'9')
                abort=1;
                keepchecking1=0;
            end
            if strfind(idn2,'5')
                keepchecking1=0;
            end
           
        end
    end


    function [frequency, stepsize, rate, overlap]=parse(datasignal)
        % Have to initialize these values in case of evoking return
        % statement
        frequency='';
        stepsize='';
        rate='';
        overlap='';
        param_check=0;

        
        temp=datasignal;  
        hi_index=strfind(temp,'hi');
        bye_index=strfind(temp,'bye');
        
        % if data between contains 3 stim numbers
        if (bye_index-+hi_index)==12
            param_check=1;
        else
            return
        end
       
        frequency=str2double(temp(hi_index+2:hi_index+6));
        stepsize=str2double(temp(hi_index+7:hi_index+8));
        rate=str2double(temp(hi_index+9:hi_index+10));
        overlap=str2double(temp(hi_index+11));

             
    end

fclose(s2);
delete(s2);
clear s2
close all
clear -frequencies amplitudes

disp('DONE!')

end



