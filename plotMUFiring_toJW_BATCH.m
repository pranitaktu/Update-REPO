function plotMUFiring_toJW_BATCH

filenames = dir('*JNEdecomposed.mat');
for i = 1:length(filenames)
    plotMUFiring_toJW(filenames(i).name)
end