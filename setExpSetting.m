
clear all

subNum = 1; 

rng(datenum(date)*subNum) %generate rand label
unSS = [1 3]; %unique set size 1-3
unDel = [1 7]; %unique delay 1-7s

nTrials = 700;

ss = repelem([unSS(1) unSS(2)], 1, nTrials/2); %each row is 1 trial, with cols: # of oriented items, s of delay
delays = repelem([unDel(1) unDel(2) unDel(1) unDel(2)], 1, nTrials/4);


nBins = 25;
binEdges = linspace(-pi/2, pi/2, nBins+1);
% target distribution (uniform for now)
for ij = 1:nBins
    ori(ij,:) =  binEdges(ij) + (binEdges(ij+1)-binEdges(ij)).*rand((length(unSS) * length(unDel))*7,1); 
end

%%inside rand: calculates the number of orientations to sample from each
%%bin, based on # trials wanted

expMat = [ss' delays' ori(:)];
indi = randperm(length(expMat));
expMat = expMat(indi,:);

%%
block = 1;
bayesDriftExp(expMat,block,subNum) 




