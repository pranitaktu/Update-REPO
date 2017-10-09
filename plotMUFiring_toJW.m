function [] = plotMUFiring_toJW(filename)
% plotMUFiring    loads motor units decomposed from cluster; plots units
%
% USAGE
%    plotMUFiring('071816_52405_Sol_4_JNEdecomposed.mat')
%    
%   

[tableau] = Tableau;
load(filename)

filename_index = strfind(filename,'_');

load([filename(1:filename_index(end)-1) '.mat'])

try
    Torque_lowess = smooth(Torque,fsamp/10,'loess')';
catch
    Torque_lowess = zeros(1,length(EMG));
    Torque = zeros(1,length(EMG));
end


[B,A] = butter(3,[10]*2/fsamp,'low');
Torque_butter = filtfilt(B,A,Torque')';
 
Torque_time = 1/fsamp:1/fsamp:length(Torque)/fsamp;

h=figure;
t=1;

for j = 1:size(MUFiring,2)
    MUFiring_loop = [];
    MUFiring_loop = MUFiring{j};

    ISI = diff(MUFiring_loop)/fsamp; % Interspike interval in sec
    IDR = 1./ISI; % instantaneous discharge rate
    
    MUTime = MUFiring_loop(2:end)/fsamp; % Time of discharge in sec
    
    %%
    
  
      display(min(MUFiring{j}))
    
    
    %%
    
    subplot(2,1,2);
    plot(MUTime,IDR-(j-1)*20,'.','MarkerSize',8,'Color',tableau(t,:));hold all    % offset axis
    xlim([0 max(Torque_time)]);
    t=t+2;
    if t>20
        t=1;
    end
end

%%%


%%%
subplot(2,1,1)
plot(Torque_time,Torque);hold all
plot(Torque_time,Torque_lowess,'--','LineWidth',2);
plot(Torque_time,Torque_butter,'LineWidth',2);
%

%
try
axis([0 max(Torque_time) min(Torque) max(Torque)]);
catch
    axis([0 max(Torque_time) 0 1])
end

p = get(subplot(2,1,1),'position');
p(2) = p(2)*1.4;  % bottom
p(4) = p(4)*0.4;  % height
set(subplot(2,1,1), 'position', p);
xlim([0 max(Torque_time)]);
box off

title(filename(1:end-4),'interpreter','none')

subplot(2,1,2);

y = get(subplot(2,1,2),'ylim');

if y(2)>100
    y(2)=100;
end

if ~isempty(MUFiring)

y(1) = -(j-1).*20;

set(subplot(2,1,2),'ylim',y)

p = get(subplot(2,1,2),'position');
p(2) = p(2)*.1; % Subtracts 90 percent from bottom
p(4) = p(4)*2.3;  % Add 200 percent to height
set(subplot(2,1,2), 'position', p);
set(gca,'YGrid','on');
set(gca,'GridLineStyle','-');
set(gca, 'YColor', [.4, .4, .4]);

Ticklines = [(length(MUFiring)*20*-1):20:20];


set(gca,'YTick',Ticklines);
end

box off

% set(gca,'XTickLabel',Ticklinevalues)
%
warning('off','MATLAB:print:FileName')
savefilename = strcat(filename(1:end-4), '.pdf');
print (h,'-dpdf',savefilename);

pause(1)

%close(h)

end

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
    0.6196    0.8549    0.8980]; % load colormap
end


