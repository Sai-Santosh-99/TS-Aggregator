close all;
clear all;
clc;

probo=0:1/50:1;
size_probo=51;
horizon = 10000;
index = 1:1:horizon;

n_tests = 50;

regret_ucb = zeros(n_tests, 10001);
regret_bts_120 = zeros(n_tests, 10001);
regret_bts_2 = zeros(n_tests, 10001);

mean_T = zeros(n_tests, 4);

for x=1:n_tests

    index_1 = randperm(size_probo,6);
    
    prob(x,:) = [0.54, 0.53, 0.52, 0.51];
    
    disp(x);
    [best_arm_mu best_arm] = max(prob(x,:));
    
    X = ones(1,length(prob(x,:)));
    T = ones(1,length(prob(x,:)));
    Q = zeros(1, length(prob(x,:)));
    N = length(prob(x,:));
    
    bin_count=100;
    B= zeros(bin_count, length(prob(x,:)));
    
    for i=1:horizon      
        if(i==1)
             for j=1:length(prob(x,:))
                 B(:,j) = B_update(B(:,j),bin_count);
             end
        end
        for j=1:length(prob(x,:))
            for ka=1:1:1
                val = find(B(:,j)>0);   

                if (find (val > 0))
                    val_ind = randi(size(val,1));
                    B(val(val_ind),j) = B(val(val_ind),j)-1;
                    B(:,j) = B_update(B(:,j),bin_count);
                end
            end
            Q(j) = machine_bts_2(X(j), T(j), bin_count,B(:,j));
        end

        m = max(Q);
        I = find(Q == m);
        r = rand(1, length(prob(x,:)));

        if r(I(1)) < prob(x,I(1))
            X(I(1)) = X(I(1)) + 1;
        end

        regret_bts_2(x, i+1) = regret_bts_2(x, i) + best_arm_mu - prob(x,I(1));

        T(I(1)) = T(I(1)) + 1;
        B(:,I(1)) = B_update(B(:,I(1)),bin_count);
        N = N + 1;
    end
    mean_T(x,:) = T;
    best_arm_sel_bts_2(x) = T(best_arm);
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

function Q = machine(X, T) 
    val = rand(1,T);
    val = sort(val);
    Q = val(X);
end

function Q = machine_bts(X, T)
    Q = betarnd(X, T-X);
end

function Q = machine_ucb(X, T,N)
    Q =X/T+sqrt((2*log(N))/T);
end

function Q = machine_bts_1(X, T, num_bins) 

    delta = 1/num_bins;
    bin_count = zeros(1,num_bins);
   % val = rand(1,T);
    
    for i=1:T
      %  a = val(i);
        a = rand(1,1);
       index = round(a*num_bins);
       if (index==0)
           index=1;
       end
       bin_count(index) = bin_count(index) + 1;
    end
    
    count = 0;
    index = 0;
    while count < X
        index = index + 1;
        count = count + bin_count(index);
    end
    Q = ((delta*index)+(delta*(index-1)))/2;
end


function bin_count = B_update(bin_count,num_bins) 

    a = rand(1,1);
   index = round(a*num_bins);
   if (index==0)
       index=1;
   end
   bin_count(index) = bin_count(index) + 1;
end


function Q = machine_bts_2(X, T, num_bins,bin_count) 

    delta = 1/num_bins;
       
    count = 0;
    index = 0;
    while count < X
        index = index + 1;
        count = count + bin_count(index);
    end
    Q = ((delta*index)+(delta*(index-1)))/2;
end