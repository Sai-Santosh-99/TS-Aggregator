clear all;
close all;
clc;

experiment = 200;
arms = 15;
probs = rand(experiment, arms); % probability distribution of all arms

for expt= 1:experiment
    
    disp(expt);

    close all;
    clearvars -except alg_sel regret_exp_agg regret_exp_ucb regret_exp_ts expt probs regret_exp_aggr_orig sel_exp_agg sel_exp_ucb sel_exp_ts sel_exp_aggr_orig
    
    prob = probs(expt, :);
    
    m = max(prob); 
    arm_best = find(prob == m);
    arm_best = arm_best(1+floor(length(arm_best)*rand));
    
    horizon = 10000; % horizon value
    best_reward =  cumsum(ones(1,horizon).*m);
    nb_arms = length(prob); % number of arms in the game
    nb_children = 2; % number of algorithms in the game

    regret_aggr = zeros(nb_children+2,horizon);

    T = ones(1, nb_arms); % the number of selections of each arm by the aggregator
    TAggr = ones(1, nb_arms); % the number of selections of each arm by the aggregator
    TUCB = ones(1, nb_arms); % the cummulative reward of each arm selected by the aggregator
    TTS = ones(1, nb_arms); % the number of selections of each arm by the aggregator

    % the trust value for each algorithm in the aggregator
    % start with uniform trust values
    trust_values = ones(1, nb_children)./nb_children;
    trust_values_aggr = ones(1, nb_children)./nb_children;

    % to keep track of the choices of each algorithm in each time slot
    choices = -100.*ones(1, nb_children);
    choices_aggr = -100.*ones(1, nb_children);

    % update the choices for each algorithm for the first slot.
    choices(1) =  policyUCB(0, 0, nb_arms, true); % arm selected by UCB
    choices(2) =  policyThompson(0, 0, nb_arms, true); % arm selected by Thompson
    choices(3) = policyUCBO(0,0,nb_arms, true);
    choices(4) = policyThompsonO(0,0,nb_arms, true);
    
    choices_aggr(1) = policyUCBAggr(0, 0, nb_arms, true);
    choices_aggr(2) = policyThompsonAggr(0, 0, nb_arms, true);
       
    i = 0;
    epoch=1;
    
    alg_sel_proposed = [0 0];
    alg_sel_existing = [0 0];
    
    for t = 1:horizon
       
        % learning rate of the aggregator. If horizon is fixed, then replace 't'
        % with the 'horizon' value, else rate gradually decreases with time
        learning_rate = sqrt(1*log(nb_children)/(t*nb_arms));

        % weighted random sampling of the arm to be played in the current slot from the choice of all
        % algorithms with trust values as weights. In hardware, you can fetch the trust value of one
        % algorithm and then spawn a random number. If number below the 
        % algorithm's trust value, choose that algorithm, else repeat the step for others.
        
        switch_till = 500; % no. of time slots till which switching should take place.
        switch_every = 2^epoch; % no. of time slots after which algorithms should switch.

        if t <= switch_till
            if i == 0
                chose = 1;
                choice_t = choices(1);
                b=1;
            else
                choice_t = choices(2);
                b=2;
            end
        else
            if t == switch_till+1
               disp(trust_values);
            end
            a = max(trust_values);
            b = find(trust_values == a);
            choice_x = choices(1:2);
            choice_t = choice_x(b(1+floor(length(b)*rand)));
        end
        alg_sel_proposed(b) = alg_sel_proposed(b)+1;
        
        choice_t_aggr = randsample(choices_aggr(1:2), 1, true, trust_values_aggr);

        % non-bernoulli rewards
        reward_ct = 0.1.*randn + prob(choice_t);
        reward_ct_aggr = 0.1.*randn + prob(choice_t_aggr);
        reward_ctucbo = 0.1.*randn + prob(choices(3));
        reward_cttso = 0.1.*randn + prob(choices(4));
        
        % update the reward of the aggregator in each time slot
        T(choice_t) = T(choice_t) + 1;
        TAggr(choice_t_aggr) = TAggr(choice_t_aggr) + 1;
        TUCB(choices(3)) = TUCB(choices(3)) + 1;
        TTS(choices(4)) = TTS(choices(4)) + 1;
                
        if(t > 1)
            regret_aggr(1,t) =  m*(t+nb_arms) - max(cumsum(T.*prob));
            regret_aggr(2,t) =  m*(t+nb_arms) - max(cumsum(TUCB.*prob)); %regret_ucbo
            regret_aggr(3,t) =  m*(t+nb_arms) - max(cumsum(TTS.*prob)); %regret_tso
            regret_aggr(4,t) =  m*(t+nb_arms) - max(cumsum(TAggr.*prob));
        end
        
        % unbiased reward value updation. 
        prob_of_observing_arm = sum(trust_values(i+1));
        reward = reward_ct./prob_of_observing_arm;
        
        prob_of_observing_arm_aggr = sum(trust_values_aggr(choices_aggr(1:2) == choice_t_aggr));
        reward_aggr = reward_ct_aggr./prob_of_observing_arm_aggr;
        
        % update the trust values for only the correct algorithms and then normalize the values to get a
        % probability of trusting a particular algorithm
        % Update trust values only till the switching is taking place.
        if t <= switch_till
            trust_values(i+1) = trust_values(i+1).*exp(learning_rate*reward);
            if i == 0 && rem(t, switch_every) == 0
                    i = 1;
                    epoch = epoch+1;
            elseif i== 1 && rem(t, switch_every) == 0
                    i = 0;
            end
            trust_values = trust_values./sum(trust_values);
        end
        
        trust_values_aggr(choices_aggr(1:2) == choice_t_aggr) = trust_values_aggr(choices_aggr(1:2) == choice_t_aggr).*exp(learning_rate*reward_aggr);
        trust_values_aggr = trust_values_aggr./sum(trust_values_aggr);

        % give the reward of the current time slot to all algorithms
        % irrespective of if they were correct or not
        choices(1) = policyUCB(reward_ct, choice_t, nb_arms, false);
        choices(2) = policyThompson(reward_ct, choice_t, nb_arms, false);
        choices(3) = policyUCBO(reward_ctucbo, choices(3), nb_arms, false);
        choices(4) = policyThompsonO(reward_cttso, choices(4), nb_arms, false);
        choices_aggr(1) = policyUCBAggr(reward_ct_aggr, choice_t_aggr, nb_arms, false);
        choices_aggr(2) = policyThompsonAggr(reward_ct_aggr, choice_t_aggr, nb_arms, false);
        
        [sda b]= max(trust_values_aggr);
        alg_sel_existing(b) = alg_sel_existing(b)+1;
        
    end
            
    regret_exp_agg(expt,:) =  regret_aggr(1,1:100:horizon);
    regret_exp_ucb(expt,:) =  regret_aggr(2,1:100:horizon);
    regret_exp_ts(expt,:) =  regret_aggr(3,1:100:horizon);
    regret_exp_aggr_orig(expt,:) =  regret_aggr(4,1:100:horizon);
    sel_exp_agg(expt,:) =  T;
    sel_exp_ucb(expt,:) = TUCB ;
    sel_exp_ts(expt,:) =  TTS;
    sel_exp_aggr_orig(expt,:) = TAggr;
    alg_sel(expt,1:2)= alg_sel_existing;
    alg_sel(expt,3:4)= alg_sel_proposed;
end

disp([mean(regret_exp_agg(:, end)) mean(regret_exp_ucb(:, end)) mean(regret_exp_ts(:, end)) ...
    mean(regret_exp_aggr_orig(:, end))]);

figure;
set(gcf, 'color', 'w');
hold on;
[a_ours, b_ours] = stdshade(regret_exp_agg, 0.5, 'r');
[a_orig, b_orig] = stdshade(regret_exp_aggr_orig, 0.5, 'g');
[a_ucb, b_ucb] = stdshade(regret_exp_ucb, 0.5, 'b');
[a_ts, b_ts] = stdshade(regret_exp_ts, 0.5, 'm');
hold off;
legend([a_ours, a_orig, a_ucb, a_ts], 'Our Aggregator','Original Aggregator','UCB','TS');
 
reg = [regret_exp_agg(:, end) regret_exp_aggr_orig(:, end) regret_exp_ucb(:, end) regret_exp_ts(:, end)];

figure;
set(gcf, 'color', 'w');
bar(reg);
legend('Our Aggregator','Original Aggregator','UCB','TS');

function choiceUCB = policyUCB(reward, arm, nb_arms, start_algo)
    persistent X;
    persistent T;
    persistent N;
    persistent choiceUCB_prev;
    if start_algo == true
        X = ones(1, nb_arms);
        T = ones(1, nb_arms);
        N = nb_arms;
    else
        X(arm) = X(arm) + reward;
        T(arm) = T(arm) + 1;
        N = N + 1;
    end
        ucb = X./T + sqrt(2*log(N)./T);
        m = max(ucb); 
        arm = find(ucb == m);
        arm = arm(1+floor(length(arm)*rand));
    

    choiceUCB = arm;
    choiceUCB_prev = arm;
end

function choiceUCBAggr = policyUCBAggr(reward, arm, nb_arms, start_algo)
    persistent XAggr;
    persistent TAggr;
    persistent NAggr;
    if start_algo == true
        XAggr = ones(1, nb_arms);
        TAggr = ones(1, nb_arms);
        NAggr = nb_arms;
    else
        XAggr(arm) = XAggr(arm) + reward;
        TAggr(arm) = TAggr(arm) + 1;
        NAggr = NAggr + 1;
    end
        ucb = XAggr./TAggr + sqrt(2*log(NAggr)./TAggr);
        m = max(ucb); 
        arm = find(ucb == m);
        arm = arm(1+floor(length(arm)*rand));
    

    choiceUCBAggr = arm;
end

function choiceUCBO = policyUCBO(reward, arm, nb_arms, start_algo)
    persistent XO;
    persistent TO;
    persistent NO;
    persistent choiceUCB_prevO;
    if start_algo == true
        XO = ones(1, nb_arms);
        TO = ones(1, nb_arms);
        NO = nb_arms;
    else
        XO(arm) = XO(arm) + reward;
        TO(arm) = TO(arm) + 1;
        NO = NO + 1;
    end
        ucb = XO./TO + sqrt(2*log(NO)./TO);
        m = max(ucb); 
        arm = find(ucb == m);
        arm = arm(1+floor(length(arm)*rand));
    

    choiceUCBO = arm;
    choiceUCB_prevO = arm;
end

function choiceThompson = policyThompson(reward, arm, nb_arms, start_algo)

    persistent XTs;
    persistent TTs;
    persistent NTs;
    persistent persistentPrev;
    if start_algo == true
        XTs = ones(1, nb_arms);
        TTs = ones(1, nb_arms);
        NTs = nb_arms;
    else
        XTs(arm) = XTs(arm) + reward;
        TTs(arm) = TTs(arm) + 1;
        NTs = NTs + 1;
    end    
        ucb = betarnd(XTs, TTs-XTs);
        m = max(ucb); 
        arm = find(ucb == m);
        arm = arm(1+floor(length(arm)*rand));
    
    
    choiceThompson = arm;
    persistentPrev = arm;
end

function choiceThompsonO = policyThompsonO(reward, arm, nb_arms, start_algo)

    persistent XTsO;
    persistent TTsO;
    persistent NTsO;
    persistent persistentPrevO;
    if start_algo == true
        XTsO = ones(1, nb_arms);
        TTsO = ones(1, nb_arms);
        NTsO = nb_arms;
    else
        XTsO(arm) = XTsO(arm) + reward;
        TTsO(arm) = TTsO(arm) + 1;
        NTsO = NTsO + 1;
    end    
        ucb = betarnd(XTsO, TTsO-XTsO);
        m = max(ucb); 
        arm = find(ucb == m);
        arm = arm(1+floor(length(arm)*rand));
    
    
    choiceThompsonO = arm;
    persistentPrevO = arm;
end

function choiceThompsonAggr = policyThompsonAggr(reward, arm, nb_arms, start_algo)
    persistent XTsAggr;
    persistent TTsAggr;
    persistent NTsAggr;
    if start_algo == true
        XTsAggr = ones(1, nb_arms);
        TTsAggr = ones(1, nb_arms);
        NTsAggr = nb_arms;
    else
        XTsAggr(arm) = XTsAggr(arm) + reward;
        TTsAggr(arm) = TTsAggr(arm) + 1;
        NTsAggr = NTsAggr + 1;
    end    
        ucb = betarnd(XTsAggr, TTsAggr-XTsAggr);
        m = max(ucb); 
        arm = find(ucb == m);
        arm = arm(1+floor(length(arm)*rand));
    
    
    choiceThompsonAggr = arm;
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
F = 1:100:10000;
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