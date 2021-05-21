close all;
clear all;
clc;

n_tests = 50;

prob_x = rand(n_tests, 15);
mean_T = zeros(n_tests,size(prob_x,2));

cumReward = zeros(1, 10001);
idealReward = zeros(1, 10001);

horizon = 10000;
index = 1:1:horizon;

num_bins = 30;

regret_ts = zeros(n_tests, horizon+1);
regret_bts = zeros(n_tests, horizon+1);
reward_ts = zeros(n_tests, horizon+1);
reward_bts = zeros(n_tests, horizon+1);

for x=1:n_tests
    
    disp(x);
    prob = prob_x(x,:);
            
    X = ones(1,length(prob));
    T = ones(1,length(prob));
    N = length(prob);
    Q = zeros(1,length(prob));

    for i=1:horizon

        for j=1:length(prob)
            Q(j) = machine(X(j), T(j), num_bins);
        end

        m = max(Q);
        I = find(Q == m);
        r = rand(1, length(prob));

        if r(I(1)) < prob(I(1))
            X(I(1)) = X(I(1)) + 1;
            reward_ts(x, i+1) = reward_ts(x, i) + 1;
        else
            reward_ts(x, i+1) = reward_ts(x, i);
        end

        regret_ts(x, i+1) = regret_ts(x, i) + max(prob) - prob(I(1));

        T(I(1)) = T(I(1)) + 1;
        N = N + 1;
    end
end

function Q = machine(X, T, num_bins) 

    delta = 1/num_bins;
    bin_count = zeros(1,num_bins);
    val = rand(1,T);
    
    for i=1:T
        a = val(i);
        for j=1:num_bins
            if (a>(j-1)*delta) && (a<j*delta)
                bin_count(j) = bin_count(j) + 1;
                break;
            end
        end
    end
    
    count = 0;
    index = 0;
    while count < X
        index = index + 1;
        count = count + bin_count(index);
    end
    Q = ((delta*index)+(delta*(index-1)))/2;
end

function [lineOut, fillOut] = stdshade(amatrix,alpha,acolor,F,smth)
% usage: stdshading(amatrix,alpha,acolor,F,smth)
% plot mean and sem/std coming from a matrix of data, at which each row is an
% observation. sem/std is shown as shading.
% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)
% - smth defines the smoothing factor (default is no smooth)
% smusall 2010/4/23
if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end
if exist('F','var')==0 || isempty(F)
    F=1:size(amatrix,2);
end
if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1; %no smoothing by default
end  
if ne(size(F,1),1)
    F=F';
end
amean = nanmean(amatrix,1); %get man over first dimension
if smth > 1
    amean = boxFilter(nanmean(amatrix,1),smth); %use boxfilter to smooth data
end
% astd = nanstd(amatrix,[],1); % to get std shading
astd = nanstd(amatrix,[],1)/sqrt(size(amatrix,1)); % to get sem shading
if exist('alpha','var')==0 || isempty(alpha) 
    fillOut = fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none');
    acolor='k';
else
    fillOut = fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none');
end
if ishold==0
    check=true; else check=false;
end
hold on;
lineOut = plot(F,amean, 'color', acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line
if check
    hold off;
end
end

function dataOut = boxFilter(dataIn, fWidth)
% apply 1-D boxcar filter for smoothing
fWidth = fWidth - 1 + mod(fWidth,2); %make sure filter length is odd
dataStart = cumsum(dataIn(1:fWidth-2),2);
dataStart = dataStart(1:2:end) ./ (1:2:(fWidth-2));
dataEnd = cumsum(dataIn(length(dataIn):-1:length(dataIn)-fWidth+3),2);
dataEnd = dataEnd(end:-2:1) ./ (fWidth-2:-2:1);
dataOut = conv(dataIn,ones(fWidth,1)/fWidth,'full');
dataOut = [dataStart,dataOut(fWidth:end-fWidth+1),dataEnd];
end