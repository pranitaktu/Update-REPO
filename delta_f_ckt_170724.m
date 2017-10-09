% Determine delta f values from decomposed MU spike trains.
% Plot MU and delta f comparisons.
% Save results as PDF and Excel files.

function delta_f_ckt_170724(filename)
close all
% Open Decomposed MU File

load(filename);

if exist('fsamp') == 0
    fsamp = 2048;
end

% Plot Motor units across all contractions

for i = 1:length(MUFiring)
    plot (MUFiring{i}(2:end)./fsamp,1./diff(MUFiring{i})*fsamp-(10*(i-1)),'.'); hold all
end

% Display plot legend for different individual MUs
legend('show')

% TRUEMU is representative of the Total number of MUs present in contraciton
TRUEMU = [1:length(MUFiring)];
% Prompt user to exclude MU from analysis. Example: [3 7]
badMU = input('Enter bad motor units:');
% MUs that were input are removed from MUFiring array
MUFiring(badMU) = [];
% Clear Figure
clf
% Reassign TRUEMU varibale to represent user removed units
TRUEMU(badMU) = [];

% Plot MUs and prompt user for time frame for analysis
for i = 1:length(MUFiring)
    plot (MUFiring{i}(2:end)./fsamp,1./diff(MUFiring{i})*fsamp-(10*(i-1)),'.'); hold all
end

MUTime = input('Enter time for analysis (e.g. [0 30]):');

%

for i = 1:length(MUFiring)
    % Create cell array of MU datapoints in seconds
    loopMUFiring{i} = MUFiring{i}/fsamp;
    % Target MU firing instances within bounds of time frame selected above
    loopMUFiring{i} = loopMUFiring{i}(loopMUFiring{i}>MUTime(1) & loopMUFiring{i}<MUTime(2)).*fsamp;
    %%%Remove MU spikes with differential greater than 50ms
    loopMUFiring{i}(1./(diff(loopMUFiring{i})/fsamp)>50) = [];
    % Find Interspike Intervals greater than .4 seconds, assign to badISI variable
    badISIs = (find((diff(loopMUFiring{i})/fsamp)>.4));
    % Within the first 25% of data file, find max badISI value
    badISIsstart = max(badISIs(badISIs<.25*length(loopMUFiring{i})));
    %
    badISIs = (find((diff(loopMUFiring{i})/fsamp)>.8));
    %
    badISIsend = min(badISIs(badISIs>.75*length(loopMUFiring{i})));
    % If badISI varibale is blank, badISIstart defaults to 0, badISIend defaults to length of MUFiring
    if isempty(badISIsstart) == 1
        badISIsstart = 1;
    end
    if isempty(badISIsend) == 1
        badISIsend = length(loopMUFiring{i});
    end
    
    loopMUFiring{i} = loopMUFiring{i}(badISIsstart:badISIsend);
end
%
TRUEMU = [TRUEMU(~cellfun(@isempty, loopMUFiring)),0];
MUFiring = loopMUFiring(~cellfun(@isempty, loopMUFiring));

% Index MU

for MU = 1:size(MUFiring,2),
    index = MUFiring{MU};
    firing(MU,index) = 1;
end

firing(:,length(firing)+fsamp)=0;
AVERAGE = 1./sum(hanning(fsamp*2))*filtfilt(hanning(fsamp*2),1,firing')';
AVERAGE(size(AVERAGE,1)+1,:)=mean(AVERAGE);
figure ()
c = 1;
Tableaumap = Tableau;

for i = 1:length(MUFiring)
    windowDR{i} = AVERAGE(i,min(MUFiring{i}):max(MUFiring{i}));
    
    onset(i) = min(MUFiring{i})/fsamp;
    windowST{i} = min(MUFiring{i})/fsamp:1/fsamp:max(MUFiring{i})/fsamp;
    
    try
        plot(windowST{i},windowDR{i},'LineWidth',1.75,'Color',Tableaumap(c,:)); hold all
        c=c+1;
    catch
        c=1;
        plot(windowST{i},windowDR{i},'LineWidth',1.75,'Color',Tableaumap(c,:)); hold all
    end
end
onset';

windowDR{size(AVERAGE,1)} = AVERAGE(size(AVERAGE,1),:);
onset(size(AVERAGE,1)) = 1/fsamp;
windowST{size(AVERAGE,1)} = 1/fsamp:1/fsamp:(size(AVERAGE,2))/fsamp;

hold all;
title({[filename(1:end),'_',num2str(MUTime(1)),'_',num2str(MUTime(2))],['All ',num2str(length(MUFiring)),' units']},'interpreter','none')
c = 1;

% Plot low vs high threshold motor units, find delta f values

for i = 1:length(windowST)
    for j = 1:length(windowST)
        
        if i == j
            continue
        end
        
        Tdiff = min(windowST{j}) - min(windowST{i});
        if Tdiff>0 % && min(windowST{j})<max(windowST{i})
            try
                figure (4);clf;
                plot(windowST{i},windowDR{i},'b','linewidth',2);hold all
                plot(MUFiring{i}(2:end)/fsamp,1./(diff(MUFiring{i})/fsamp),'b.')
                plot(windowST{j},windowDR{j}+10,'r','linewidth',2)
                plot(MUFiring{j}(2:end)/fsamp,1./(diff(MUFiring{j})/fsamp)+10,'r.')
                ylim([0 30])
            end
            
            Fstart = windowDR{i}(round(Tdiff*fsamp));
            Tstart = windowST{i}(windowDR{i}==Fstart);
            
            Fmax = max(windowDR{i}(round(Tdiff*fsamp):end));
            Tmax = windowST{i}(windowDR{i}==Fmax);
            
            ratemod = Fmax - Fstart;
            
            if max(windowST{j})<max(windowST{i})
                Tend = max(windowST{j});
            else
                Tend = max(windowST{i});
            end
            
            
            Fend = windowDR{i}(round(windowST{i}*fsamp)==round(Tend*fsamp));
            deltaf =  Fstart - Fend;
            
            try
                plot([Tstart windowST{j}(1)],[Fstart, windowDR{j}(1)+10],'k--')
                plot([Tend windowST{j}(end)],[Fend, windowDR{j}(end)+10],'k--')
                plot([Tstart windowST{i}(end)+1],[Fstart Fstart],'k:')
                plot([Tstart windowST{i}(end)+1],[Fend Fend],'k:')
                plot([Tmax Tmax],[Fmax Fstart],'k','linewidth',2)
                
                ratemodtxt = sprintf('%0.1f pps',ratemod);
                text(Tmax+1,Fmax+1,ratemodtxt)
                plot([Tend+1 windowST{j}(end)+1],[Fend, Fstart],'k','linewidth',2)
                deltaftxt = sprintf('%0.1f pps',deltaf);
                text(Tend+1.2,(Fend+Fstart)/2,deltaftxt)
                plot([windowST{i}(1) windowST{j}(1)],[windowDR{j}(1)+10 windowDR{j}(1)+10],'k','linewidth',2)
                Tdifftxt = sprintf('%0.1f s',Tdiff);
                text(windowST{i}(1),windowDR{j}(1)+9,Tdifftxt)
            end
            
            if length(deltaf)==1
            else
                deltaf=0;
            end
            
            initiali = windowDR{i}(1);
            maxi = max(windowDR{i});
            finali = windowDR{i}(end);
            initialj = windowDR{j}(1);
            maxj = max(windowDR{j});
            finalj = windowDR{j}(end);
            
            try
                
                %%% find peak on CST
                
                [M,I] = max(AVERAGE(end,:));
                
                %%% Segment windowDRicut
                
                thresh = I/fsamp;
                
                windowDRicut = windowDR{i}((windowST{i}>=thresh));
                windowDRjcut = windowDR{j}((windowST{j}>=thresh));
                
                windowSTicut = windowST{i}((windowST{i}>=thresh));
                windowSTjcut = windowST{j}((windowST{j}>=thresh));
                
                if windowSTjcut(1)<thresh
                    windowDRicut = windowDRicut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    windowDRjcut = windowDRjcut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    
                    windowSTicut = windowSTicut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    windowSTjcut = windowSTjcut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    
                else
                    windowDRicut = windowDRicut(round(1+((windowSTjcut(1)-thresh)*fsamp):min([length(windowDRicut) length(windowDRjcut)])+((windowSTjcut(1)-thresh)*fsamp)));
                    windowDRjcut = windowDRjcut(1:min([length(windowDRicut) length(windowDRjcut)]));
                    
                    windowSTicut = windowSTicut(round(1+((windowSTjcut(1)-thresh)*fsamp):min([length(windowDRicut) length(windowDRjcut)])+((windowSTjcut(1)-thresh)*fsamp)));
                    windowSTjcut = windowSTjcut(1:min([length(windowDRicut) length(windowDRjcut)]));
                end
                
                plot(windowSTicut,windowDRicut,'b','linewidth',4)
                plot(windowSTjcut,windowDRjcut+10,'r','linewidth',4)
                
                %windowDRicut = windowDR{i}(round(Tdiff*fsamp):round((Tend*fsamp)-(min(windowST{i})*fsamp)));
                [rvalue slope inter] = regression(windowDRicut,windowDRjcut);
            catch
                rvalue = 0;
            end
            
            titletxt1 = sprintf('Units %0.0f and %0.0f - deltaf = %0.1f, rate-rate slope = %0.2f, rate-rate r-value = %0.2f',TRUEMU(i),TRUEMU(j),deltaf,slope,rvalue);
            titletxt2 = sprintf('Unit %0.0f - Tstart = %0.1f, Tend = %0.1f, Fstart = %0.1f, Fmax = %0.1f, Fend = %0.1f',TRUEMU(i),windowST{i}(1),windowST{i}(end),initiali,maxi,finali);
            titletxt3 = sprintf('Unit %0.0f - Tstart = %0.1f, Tend = %0.1f, Fstart = %0.1f, Fmax = %0.1f, Fend = %0.1f',TRUEMU(j),windowST{j}(1),windowST{j}(end),initialj,maxj,finalj);
            title({titletxt1,titletxt2,titletxt3})
            
            if or(Tdiff > 0.5 && rvalue > 0.7 && ratemod > 0.5,i==length(windowST));
                output(c,:) = [TRUEMU(i),TRUEMU(j),deltaf,Tdiff,rvalue,slope,ratemod,initiali,maxi,finali,initialj,maxj,finalj];
                c=c+1;
                
                
                 if i~=length(windowST) %%%% CHECK WHY THIS IS HERE
                 end
                 
                if i~=length(windowST)
                    %%%Folder Loop%%%
                    currentdir = pwd;
                    pdfdir = [currentdir '\deltaf_pdf'];
                    if ~exist('deltaf_pdf', 'dir')
                        mkdir('deltaf_pdf');
                    end
                    cd (pdfdir)                    
                    print([filename(1:end),'_',num2str(MUTime(1)),'-',num2str(MUTime(2)),'_Units_', num2str(TRUEMU(i), '%02i'),'_vs_',num2str(TRUEMU(j), '%02i')],'-dpdf')
                    cd (currentdir)
                end
            else
            end
        else
        end
    end
end

plot_axis = [MUTime(1) MUTime(2) 0 20];
axis(plot_axis)
savename = [filename(1:end),'_all units_win'];
c=1;


                    currentdir = pwd;
                    pdfdir = [currentdir '\deltaf_pdf'];
                    if ~exist('deltaf_pdf', 'dir')
                        mkdir('deltaf_pdf');
                    end

for i = 1:length(MUFiring)
    cd (pdfdir)
    %hold all; %plot(1/fsamp:1/fsamp:length(AVERAGE)/fsamp,AVERAGE(1,:))
    starttime = MUFiring{i}(1)/fsamp;
    endtime = MUFiring{i}(end)/fsamp;
    if c<21
        figure;plot(MUFiring{i}(2:end)/fsamp,1./(diff(MUFiring{i})/fsamp),'o','Color',Tableaumap(c,:))
        hold on; plot(starttime:1/fsamp:endtime,AVERAGE(i,starttime*fsamp:endtime*fsamp),'Color',Tableaumap(c,:),'LineWidth', 2)
    else
        c=1;
        figure;plot(MUFiring{i}(2:end)/fsamp,1./(diff(MUFiring{i})/fsamp),'o','Color',Tableaumap(c,:))
        hold on; plot(starttime:1/fsamp:endtime,AVERAGE(i,starttime*fsamp:endtime*fsamp),'Color',Tableaumap(c,:),'LineWidth', 2)
    end
    title({filename(1:end), ['Unit ', num2str(TRUEMU(i), '%02i')]}, 'Interpreter', 'none');
    axis(plot_axis)
    print([filename(1:end),'_Unit_', num2str(TRUEMU(i), '%02i'),'_win'],'-dpdf')
    c=c+1;
end
cd (currentdir)
% Save the data to excel file in folder "results_excel" in current directory

try
    currentdir = pwd;
    exdir = [currentdir '\results_excel'];
    if exist('results_excel')== 7
        cd (exdir)
        dlmwrite([filename(1:end),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_deltaf.xls'],output,'delimiter','\t')
        cd (currentdir)
    else
        mkdir('results_excel')
        cd (exdir)
        dlmwrite([filename(1:end),'_',num2str(MUTime(1)),'_',num2str(MUTime(2)),'_deltaf.xls'],output,'delimiter','\t')
        cd (currentdir)
    end
end
end

% Load colormap

function [tableau] = Tableau
tableau =     [0.1216    0.4667    0.7059
    0.6824    0.7804    0.9098
    1.0000    0.4980    0.0549
    1.0000    0.7333    0.4706
    0.1725    0.6275    0.1725
    0.5961    0.8745    0.5412
    0.8392    0.1529    0.1569
    1.0000    0.5961    0.5882
    0.5804    0.4039    0.7412
    0.7725    0.6902    0.8353
    0.5490    0.3373    0.2941
    0.7686    0.6118    0.5804
    0.8902    0.4667    0.7608
    0.9686    0.7137    0.8235
    0.4980    0.4980    0.4980
    0.7804    0.7804    0.7804
    0.7373    0.7412    0.1333
    0.8588    0.8588    0.5529
    0.0902    0.7451    0.8118
    0.6196    0.8549    0.8980];
end