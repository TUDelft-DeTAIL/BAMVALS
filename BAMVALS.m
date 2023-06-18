%% SISO Watertank Benchmark

%configuration is set as it is used in the experiments in the paper.
%for the code to run you must download the following files:

%https://github.com/kbatseli/MVMALS [1]
%https://data.4tu.nl/articles/_/12960104 [2]


%[1] Batselier, Kim, Zhongming Chen, and Ngai Wong. "Tensor Network alternating linear scheme for MIMO Volterra system identification." Automatica 84 (2017): 26-35.
%[2] Schoukens, Maarten, et al. "Cascaded tanks benchmark combining soft and hard nonlinearities." Workshop on nonlinear system identification benchmarks. 2016.


clear all
close all
load dataBenchmark %see [2]
%%

M = 95;                     %memory of Volterra system
yValOr = yVal;              %storinng true outputs
order = 3;                  %actual order plus one to account for expected format of rank array in mvlas algorithm
normcore = 2;               %location of norm - prior knowledge will be imposed in this core
sigma_sq = 0.33/2*10^-1;    %value for sigma_square

ITR = 1;                    %for multiple runs, change it to desired number of iterations

% RKS = [3,6,9];            %ranks used in the paper. Split in severala
% RKS = [72,96];            %arrays for quicker runtime
RKS = [48];                 %tank used for Fig5

LMBD = [20000];             %array of diffrerent hyperparameters 1/lambda

%%
for ranks = 1:length(RKS)
    % allocate arrays
    r = RKS(ranks).*ones(1,order);  %rank array as passed to mvals

    time_mvlas = zeros(1,ITR);      %store tic/toc mvals
    rel_err_mvals = zeros(1,ITR);   %store rel error mvlas
    RMSE_mvals = zeros(1,ITR);      %store RMSE mals

    time_tikh = zeros(length(LMBD),ITR);    %store tic/toc update last core
    rel_err_tikh = zeros(length(LMBD),ITR); %store rel error bamvals
    RMSE_tikh = zeros(length(LMBD),ITR);    %store RMSE bamvals
    for i = 1:ITR
        %Use MVALS, see [1]
        tic;
        [TT,e2, Vm,Vp]=mvals(yEst,uEst,M,r(2:end-1),normcore);
        time_mvals(i) = toc;
        yhat=sim_volterraTN(uVal,TT);

        %Compute relative error and RMSE for MVALS
        rel_err_mv = norm(yhat(M:end) - yValOr(M:end))/norm(yValOr(M:end));
        RMSE_mv = sqrt(mean((yhat(M:end)-yValOr(M:end)).^2));
        disp(['relative error mvals: ', num2str(rel_err_mv)])
        disp(['RMSE mvals: ', num2str(RMSE_mv)])
        rel_err_mvals(i) = rel_err_mv;
        RMSE_mvals(i) = RMSE_mv;

        %iterate over different 1/lambda
        for l = 1:length(LMBD)
            lambda = LMBD(l); 
                   
            tic;
            %see equation (19)
            PhiW = compute_PhiW(uEst,TT,Vm,Vp,r(normcore-1:normcore));
            I = eye(size(PhiW,2));
            Pplus = PhiW'*PhiW/sigma_sq + lambda*I;
            %see equation (20)
            wd_tikh = Pplus\PhiW'*yEst/sigma_sq;

            time_tikh(l,i) = toc;

            %update core with results from equation (19) and (20)    
            TT_tikh = TT;
            TT_tikh.core{normcore} = reshape(wd_tikh,size(TT.core{normcore}));
            yhat_tikh = sim_volterraTN(uVal,TT_tikh);

            rel_err_ti = norm(yhat_tikh(M:end) - yValOr(M:end))/norm(yValOr(M:end));
            RMSE_ti = sqrt(mean((yhat_tikh(M:end)-yValOr(M:end)).^2));
        %     disp(['relative error tikh: ', num2str(rel_err_ti)])
        %     disp(['RMSE tikh: ', num2str(RMSE_ti)])
            rel_err_tikh(l,i) = rel_err_ti;
            RMSE_tikh(l,i) = RMSE_ti;
        end             
    end
    disp('=========================================')
        disp(['R = ', num2str(RKS(ranks))])
        disp(['runtime mvals = ', num2str(mean(time_mvals)), ' +- ', num2str(std(time_mvals))])
        disp(['rel err mvals = ', num2str(mean(rel_err_mvals)), ' +- ', num2str(std(rel_err_mvals))])
        disp(['RMSE mvals = ', num2str(mean(RMSE_mvals)), ' +- ', num2str(std(RMSE_mvals))])
    for l = 1:length(LMBD)
        disp('-----------------------------------------')
        disp(['lambda = ', num2str(LMBD(l))])        
        disp(['runtime tikh = ', num2str(mean(time_tikh(l,:))), ' +- ', num2str(std(time_tikh(l,:)))])
        disp(['runtime mvals + tikh = ', num2str(mean(time_mvals + time_tikh(l,:))), ' +- ', num2str(std(time_mvals + time_tikh(l,:)))])           
        disp(['rel err tikh = ', num2str(mean(rel_err_tikh(l,:))), ' +- ', num2str(std(rel_err_tikh(l,:)))])
        disp(['RMSE tikh = ', num2str(mean(RMSE_tikh(l,:))), ' +- ', num2str(std(RMSE_tikh(l,:)))])
    end
end



%% create plot for Fig. 5
grid on
hold on
plot(yValOr(M:end));

plot(yhat_tikh(M:end));
legend('original','tikh')

% compute confidence bounds
PhiStarW = compute_PhiW(uVal,TT,Vm,Vp,r(normcore-1:normcore));
confidenceMatrix = PhiStarW*inv(Pplus)*PhiStarW';
confidenceBounds = diag(confidenceMatrix(M:end,M:end));


x = [1:length(yhat_tikh(M:end))]';
yTop = yhat_tikh(M:end)+ confidenceBounds;
yBottom = yhat_tikh(M:end)- confidenceBounds;

fill([x;flipud(x)], [yTop; flipud(yBottom)], 'g', 'FaceAlpha',0.5);


