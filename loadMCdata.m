
function [] = loadMCdata(directory, matrix, totalchan, matrixnames, trial, fsamp, trialname)
% function [] = loadMCdata
% loadMCdata    load and plot multi-channel data; remove bad channels
%
% USAGE
%    loadMCdata(directory,filename,{[input 1 range] [input 2 range]... [input N range]},totalchan,{[input 1 label]...[input N label]},trial);
%    loadMCdata('Francesco','129frblnr130205142901',{[1:64],[65:128],[129:192],[193]},193,{'ADL' 'BIC' 'WF' 'Force'},1);
% %% this is a random line I am just adding to test stuff out%%
%
% DESCRIPTION
%
%   Loads OT Biolab .otb file and unzips .sig file. Outputs channels with
%   possible saturation (y-value > 2000).  Plots total number of inputs;
%   multi-channel matrices are plotted in groups of 16. Prompts user for
%   bad channels to remove.  Saves data with bad channels removed.
%
%   Implemented file structure for data processing; threshold analysis
%   moved to 'plotMC'; reduced memory overhead by eliminating 'fopen'
%   command; remove saturated channels automatically
%
%
% AUTHOR    Laura Miller
%           Christopher Thompson
%
% DATE CREATED      6-Feb-2013
%
% DATE MODIFIED   22-May-2013
%
%   To do: export more auxiliary channels; downsample force/auxiliary;
%   combine with 'SelectMCSeg' (or add this first);



% options.Resize='on';
% options.WindowStyle='normal';
% options.Interpreter='tex';

% prompt={'Directory, Trialname, Matrix, Totalchan, Matrixnames, Trial'}
% numlines=1;
% defaultanswer=('post methox_crossed ex with vibration','361111Ch130308154303',{[1:64],[65],[66],[67],[68]},68,{'soleus','Position','Torque','IM1','IM2'},1);
% defaultanswer=('post methox_crossed ex with vibration')
% x = inputdlg {prompt,'Input for loadMC data',numlines,defaultanswer}

% cd(strcat('C:\Users\cthompson\Desktop\Array recordings\',directory))
%%% cd(strcat('D:\Array recordings\',directory))

cd(strcat('C:\Array recordings\',directory))


mkdir (trialname) % create unzip directory

% copyfile(strcat(trialname,'.otb'),strcat('C:\Users\cthompson\Desktop\Array recordings\',directory,'\',trialname,'\'))
copyfile(strcat(trialname,'.otb'),strcat('C:\Array recordings\',directory,'\',trialname,'\'))

% cd(strcat('C:\Users\cthompson\Desktop\Array recordings\',directory,'\',trialname))
cd(strcat('C:\Array recordings\',directory,'\',trialname))
unzip(strcat(trialname,'.otb'))
delete(strcat(trialname,'.otb'))

% convert and read abstract file

% abstract=xml2struct('abstract.xml');
% a=abstract.root.signal_group.signal{1}.channels;
% b=abstract.root.signal_group.signal{1}.fsample;

% try
%     totalchan = str2num(cell2mat(struct2cell(abstract.root.signal_group.signal{1}.channels))); % total channel number
% catch
%     totalchan = str2num(cell2mat(struct2cell(abstract.root.signal_group{1}.signal{1}.channels))); % total channel number
% end
%
% try
%     fsamp = str2num(cell2mat(struct2cell(abstract.root.signal_group.signal{1}.fsample))); % sampling frequency
% catch
%     fsamp = str2num(cell2mat(struct2cell(abstract.root.signal_group{1}.signal{1}.fsample))); % sampling frequency
% end
% totalchan = str2num(cell2mat(struct2cell(abstract.root.signal_group{1}.signal{1}.channels))); % total channel number
% fsamp = str2num(cell2mat(struct2cell(abstract.root.signal_group{1}.signal{1}.fsample))); % sampling frequency

% movefile (strcat('\Users\Hornby\Desktop\Array recordings\',directory,'\abstract.xml'),...
%     strcat('\Users\Hornby\Desktop\Array recordings\',directory,'\',trialname,'\abstract.xml')) % move abstract
% movefile (strcat('\Users\Hornby\Desktop\Array recordings\',directory,'\',strcat(trialname,'.sig')),...
%     strcat('\Users\Hornby\Desktop\Array recordings\',directory,'\',trialname,'\',strcat(trialname,'.sig'))) % move .sig
% movefile (strcat('\Users\Hornby\Desktop\Array recordings\',directory,'\',strcat(trialname,'.xml')),...
%     strcat('\Users\Hornby\Desktop\Array recordings\',directory,'\',trialname,'\',strcat(trialname,'.xml'))) % move .xml

% cd(strcat('\Users\Hornby\Desktop\Array recordings\',directory,'\',trialname))

% f = fopen(strcat(trialname,'.sig'));
 %data = fread(f,[totalchan,204800],'short');

 %[a,b] = find(abs(data)>2000); % Moved to plotMC
% unique(a)
% pause

 %[x_coordstart,allbadchan, data, datafull] = plotMC(directory,trialname,matrix,totalchan,matrixnames);
[allbadchan, data] = plotMC(directory,trialname,matrix,totalchan,fsamp,matrixnames);


if allbadchan == 'n'
    return
else
    allbadchanauto = 1:totalchan;
    allbadchanauto(matrix{1}) = [];
    
    allbadchan = [allbadchan; allbadchanauto'];
    
   
    goodchan=[];
    for i = 1:(totalchan)
        if isempty(find(i==allbadchan))
            goodchan = [goodchan; i];
        end
    end
    

    % saveid = strcat(trialname,'EMG_goodch_',num2str(trial));
    % save(saveid,'EMGgood','Force','goodchan')
    
    %%% Sessanto Setup
    %EMG = data(goodchan,:);
    
    
    %%% Quattro Setup
    EMG = data(goodchan,:);
    Torque = data(totalchan-11,:);
    Velocity = data(totalchan-10,:);
    Position = data(totalchan-9,:);
    Trigger = data(totalchan-8,:);
    EMGTorque = data(totalchan-7,:);
    JR3Fx = data(totalchan-6,:);
    JR3Fy = data(totalchan-5,:);
    JR3Fz = data(totalchan-4,:);
    JR3Mx = data(totalchan-3,:);
    JR3My = data(totalchan-2,:);
    JR3Mz = data(totalchan-1,:);
    Feedback = data(totalchan,:);
          

%     respiration = dataall(totalchan-1,:);
%     BP = dataall(totalchan-2,:);
%     Force = dataall(65,:);
%     Start = x_coordstart;
%     End = x_coordend;
    
%     %EMG4=(data(totalchan-3,:));
%     %EMG3=figure; plot (data(totalchan-4,:));
%     %EMG2=figure; plot (data(totalchan-5,:));
%     %EMG1=figure; plot (data(totalchan-6,:));
%     %Thing2 = data(totalchan-7,:);
%     % Thing1 = data(totalchan-7,:);
%     Force2 = data(totalchan-7,:);
%     Force = data(totalchan-8,:);
%     Forcedown = downsample(Force,10);
%     Position = data(totalchan-9,:);
%     Positiondown = downsample(Position,10);

    % Forcedown = downsample(data(130,:),10);
    % Force = data(130,:);
    % Positiondown = downsample(data(129,:),10);
    % Position = data(129,:);
    % Trigger = data(totalchan,:);
    
    % Forcedown = downsample(data(130,:),10);
    % Positiondown = downsample(data(129,:),10);
    

    % Trigger = data(totalchan,:);
    
    % [b,a]=butter(2,[108/(2048/2),123/(2048/2)],'stop');
    % EMG=filtfilt(b,a,EMG')';
    
    % [b,a]=butter(2,[900/(fsamp/2)]);
    % EMG=filtfilt(b,a,EMG')';
    %
    % fsamp=fsamp/5
    
    % EMG = downsample(EMG',5)';
    
    % HighFreq = 1000;
    
    
    if max(matrix{1}) == 64 && min(matrix{1}) == 1
        %muscle = 'Sol'
        %muscle = 'R_ADel';
        muscle = 'TA'
        %muscle = 'IDL'
       % muscle = 'ECU'
        %muscle = 'Med_Bi'
        %muscle = 'Caud_TA'
        %muscle = 'R_Bic'
       % muscle = 'ECar_up'
    elseif max(matrix{1}) == 128 && min(matrix{1}) == 65
        muscle = 'ECar_low'
        %muscle = 'R_VL';
        % muscle = 'Rost_TA'
        %muscle = 'L_Bi'
       % muscle = 'DistTA'
    elseif max(matrix{1}) == 192 && min(matrix{1}) == 129
       %muscle = 'LG'
         %muscle = 'L_TA';
         %muscle = 'VM'
         muscle = 'TA3'
    elseif max(matrix{1}) == 256 && min(matrix{1}) == 193
       % muscle = 'TA'
        %muscle = 'L_VL'
        muscle = 'FF'
        %muscle = 'RF'
        elseif max(matrix{1}) == 320 && min(matrix{1}) == 257
        muscle = 'Tri'
        elseif max(matrix{1}) == 384 && min(matrix{1}) == 321
        muscle = 'FE'
    else
        disp('Check your channels')
    end

    
% muscle = 'Bicep'
    
    
    % mkdir(strcat('C:\Users\cthompson\Desktop\Array recordings\',directory,'\to_cluster'))
    mkdir(strcat('C:\Array recordings\',directory,'\to_cluster'))
    
    % cd(strcat('C:\Users\cthompson\Desktop\Array recordings\',directory,'\to_cluster'))
    cd(strcat('C:\Array recordings\',directory,'\to_cluster'))
    % saveid = strcat(directory,'_', trialname(16:20),'_soleus_',num2str(trial));
    % saveid = strcat(directory,'_', trialname(16:20),'_soleus_',num2str(trial));
    % saveid = strcat(directory,'_', trialname(16:20),'_MG_',num2str(trial));
    % saveid = strcat(directory,'_', trialname(16:20),'_soleus_',num2str(trial),'_',num2str(x_coord));
    
%    if isempty(x_coord) == 1

    %%% Sessanto Setup

       %  saveid = strcat(directory,'_', trialname(16:20),'_',muscle,'_',num2str(trial));
%     else

    %%% Quattro Setup
       saveid = strcat(directory,'_', trialname(17:21),'_',muscle,'_',num2str(trial));
%     end
    
    % LowFreq
    % Low Frequency of the bandpass filter (default = 20 Hz);
    % HighFreq
    % High Frequency of the bandpass filter (default = 500 Hz);
    % fsamp
    % Sample frequency of the data (default = 2048 Hz);
    % Finterf
    % Power line interference removal frequency (default = 60 Hz);
    % DecompRuns
    
    % HighFreq = 900;
    % DecompRuns = 100;
    
    
    %%% Quattro Setup
   
    save(saveid,'EMG','fsamp','allbadchan','Torque','Velocity','Position','Trigger','EMGTorque','JR3Fx','JR3Fy','JR3Fz','JR3Mx','JR3My','JR3Mz','Feedback')    
    %save(saveid,'EMG','fsamp','allbadchan', 'EMGTorque')

    
    %%% Sessanto Setup
 % save(saveid,'EMG','fsamp','allbadchan')

    
%    save(saveid,'EMG','fsamp','Trigger','allbadchan','DecompRuns')


end



