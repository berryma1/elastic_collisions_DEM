close all
clear all
clc

% change directory to location of the csv file 
data = load ('elastic_collisions_res.csv');
n = size(data,2)/2;

figure;
xlim([0 n+1]);
ylim([0 n+1]);
hold on
video1 = VideoWriter('elastic_collisions.avi');
open(video1)
i = 1;
while i < size(data, 1)
    for j = 1:n
        scatter(data(i,2*j-1),data(i,2*j), 200, 'filled')%'or','MarkerSize',35);
    end
    pause(0.01) % prevents getframe plots from missing data
    M = getframe(1);
    writeVideo(video1, M)
    cla
    i = i+3; % plot and record every third frame
end
close(video1)