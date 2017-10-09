% function [varargout] = plotMC(subj,trialname, matrix, numchan, varargin) % why switch to 'subj'
% function [allbadchan, data] = plotMC(subj,trialname, matrix, numchan, varargin) % why switch to 'subj'
% function [x_coords, allbadchan, data] = plotMC(subj,trialname, matrix, numchan, varargin) % why switch to 'subj'

% function [x_coordstart,allbadchan, data, datafull] = plotMC(subj,trialname, matrix, numchan, varargin) % why switch to 'subj'
function [allbadchan, data] = plotMC(subj,trialname, matrix, numchan, fsamp, varargin) % why switch to 'subj'

% plotMC    plot multi-channel data; prompts for bad channels
%
% USAGE
%    [allbadchan] = plotMC('subj','38YaUtnr130205155003.sig',{[1:64]},64,{['ADL']})
%
% DESCRIPTION
%
%   Loads OT Biolab .sig file.  Plots total number of inputs;
%   multi-channel matrices are plotted in groups of 16. Prompts user for
%   bad channels to remove.  Outputs the bad channels.

%   Transposed bad channels for
%
% AUTHOR    Laura Miller
% BUTCHER   Christopher Thompson
% 
% DATE CREATED      6-Feb-2013
%
% DATE BUTCHERED   23-Mar-2013

if nargin > 4
    matrixnames = varargin{1};
else
    matrixnames =[];
end

% fsamp=10240;
% fsamp=5120;
% fsamp=2048;

% cd(strcat('\Users\cthompson\Desktop\Array recordings\',subj,'\',trialname))
cd(strcat('C:\Array recordings\',subj,'\',trialname))
f = fopen(strcat(trialname,'.sig'));


%%%Quattro Setup
data = fread(f,[numchan+8,inf],'short');
%%% Sessanto Setup
%data = fread(f,[numchan+2,inf],'short');


% data = fread(f,[numchan,70*fsamp],'short');
% data = data(:,1:35*fsamp);
% data = data(:,19*fsamp:end);

% clf
% plot (data(numchan,:));
% hold all; plot (data(65,:)-mean(data(65,end-100:end)));
% plot (data(size(data,1),1:30*fsamp));hold all
% plot (data(16,1:30*fsamp));

% 
% title('Enter zoom location')
% [x_coordzoom y_coord] = ginput(1);
% xlim([floor(x_coordzoom-fsamp*2) ceil(x_coordzoom+fsamp*2)]);
% 
% title('Enter start location')
% [x_coordstart y_coord] = ginput(1);
% x_coordstart = ceil(x_coordstart);
% % 
% % xlim([-inf inf])
% % title('Enter end location')
% % [x_coordend y_coord] = ginput(1);
% % x_coordend = floor(x_coordend);
% % 
% clf

% dataall = data;
    x_coordstart = 1;

% if x_coordstart<fsamp*3
%     x_coordstart = 1;
% end
% 
% if x_coordend > length(data)-5120*10
%     x_coordend = length(data);
%     data = data(:,x_coordstart:x_coordend);
%     x_coordend = floor(x_coordend);
% else
%     data = data(:,x_coordstart:x_coordend);
% end

% keyboard
% x_coord = []


% dataa = fread(f,[numchan,inf],'short');
% data = dataa(:,80*fsamp:length(dataa));
% clear dataa;

% data = fread(f,[numchan,107520],'short');
% keyboard
% datafull = data;
% data = diff(data);
data = double(data(:,x_coordstart:end));
% plot(data(15,:))


% [a,b] = find(abs(data)>2000); % Moved to plotMC
% unique(a')
% pause



% allbadchan=[ans'];
allbadchan = [];%% a = 'allbadchan has been changed!!!!'

% [b,a]=butter(2,[110/(2048/2),123/(2048/2)],'stop'); 
% [b,a]=butter(2,[900/(fsamp/2)]); 
% data=filtfilt(b,a,data')'; %ADDED




% for i=1:numchan
%     dataa (i,:)= downsample(data(i,:),5);
% end
%
% data = dataa;
%
% fsamp = fsamp/5

% Comment out to not plot the data
%
figure(101)

for k = 1:length(matrix)
    current  = matrix{k}; %current = which matrix/muscle [1:64] or [65:128], etc.
    if length(current) > 1
        % current16 = {[current(1:16)] [current(17:32)] [current(33:48)] [current(49:64)] [current(65:80)] [current(81:96)] [current(97:112)] [current(113:128)]}; %Current matrix broken into 16
        current16 = {[current(1:16)] [current(17:32)] [current(33:48)] [current(49:64)]}; %Current matrix broken into 16
        
        for m = 1:length(current16)
            if allbadchan == 'n'
                break
            else
                
                for i = 1:16 %current16{m}
                    plotdata(current16{m}(i),:) = data(current16{m}(i),:) - ones(1,length(data))*(i-1)*1000;
                end
                plot((1:length(plotdata))/fsamp,plotdata(current16{m},:)'); title(strcat(num2str(matrixnames{k}),'-',num2str(current16{m}(i-15)),':',num2str(current16{m}(i)))); %ylim([-70000 5000]);
                legend(num2str(current16{m}(1)),num2str(current16{m}(2)),num2str(current16{m}(3)),num2str(current16{m}(4)),num2str(current16{m}(5)),...
                    num2str(current16{m}(6)),num2str(current16{m}(7)),num2str(current16{m}(8)),num2str(current16{m}(9)),num2str(current16{m}(10)),num2str(current16{m}(11)),...
                    num2str(current16{m}(12)),num2str(current16{m}(13)),num2str(current16{m}(14)),num2str(current16{m}(15)),num2str(current16{m}(16)));
                badchan = input('Enter bad channels or empty matrix: e.g. [1,2]:');
                if badchan == 'n'
                    allbadchan = badchan
                else
                    allbadchan = [allbadchan; badchan'];
                end
                % allbadchan = [badchan'];
                % close
                clf
            end
        end
    else
        plotdata(current,:) = data(current,:);
        % figure; plot(plotdata(current,:)); title(strcat(num2str(matrixnames{k})));
        % badchan = input('Enter bad channels or empty matrix: e.g. [1,2]:');
        if badchan == 'n'
            allbadchan = badchan
        else
            allbadchan = [allbadchan; badchan'];
        end
        % allbadchan = [badchan'];
        % pause
        % close
        clf
    end
end

%for chained differential only - comment out for monopolar
% if length(matrixnames) == 1
%     allbadchan = [allbadchan; 64];
% elseif length(matrixnames) == 2
%     allbadchan = [allbadchan; 64; 128];
% elseif length(matrixnames) == 3
%     allbadchan = [allbadchan; 64; 128; 192];
% elseif length(matrixnames) == 4
%     allbadchan = [allbadchan; 64; 128; 192; 256];
% end

allbadchan = sort(allbadchan);
% varargout{1} = allbadchan;
% varargout{2} = data;

fclose all;
