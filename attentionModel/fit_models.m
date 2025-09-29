function fit_result = fit_hybrid_bias_models(fileName, savefilename)

% read the data prepared before, fit the hybrid model on 
%fileName = 'Z:\HongliWang\Juvi_ASD Deterministic\TSC2\Summary\Results\data4modelAB-AB1.mat';
load(fileName);
subjects = unique(X.Animal_ID);
% subjects = cellfun(@(x) str2double(x(4:end)), subjects);

beta_mu = 5;
beta_sigma = 7;

%% Define priors of model parameters
% Define functions to sample model parameters from uniform distributions
beta_sample = @(x) betarnd(5, 10); % Sampling function for beta parameter
alpha_sample = @(x) unifrnd(0, 1); % Sampling function for alpha parameter
forget_sample = @(x) unifrnd(0, 1); % Sampling function for forget parameter
stick_sample = @(x) unifrnd(-1, 1); % Sampling function for stickiness parameter
bias_sample = @(x) unifrnd(0, 1); % Sampling function for bias parameter
w_sample = @(x) unifrnd(0, 1); % Sampling function for w parameter
eps_sample = @(x) unifrnd(0, 1); % Sampling function for epsilon parameter
lapse_sample = @(x) unifrnd(0, 1); % Sampling function for lapse parameter (T_1^0)
rec_sample = @(x) unifrnd(0, 1); % Sampling function for recover parameter (T_0^1)

%% set params
Ms = [];

% RL_eps 
% curr_model = [];
% curr_model.name = 'RL_epsilon';
% curr_model.pMin = [1e-6 1e-6 1e-6 -1 1e-6];
% curr_model.pMax = [20 1 1 1 1];
% curr_model.pdfs = {beta_sample, alpha_sample, alpha_sample, stick_sample, eps_sample}; % Sampling functions for model parameters
% curr_model.pnames = {'beta','alpha-','alpha+','stick','epsilon'};
% 
% Ms{1}=curr_model;

% % RL1s_lapse
% curr_model = [];
% curr_model.name = 'RL_lapsestick';
% curr_model.pMin = [1e-6 1e-6 1e-6 1e-6 1e-6 -1];
% curr_model.pMax = [20 1 1 1 1 1];
% curr_model.pdfs = {beta_sample, lapse_sample, rec_sample, alpha_sample, alpha_sample, stick_sample}; % Sampling functions for model parameters
% curr_model.pnames = {'beta','lapse','rec','alpha-','alpha+','stick'};
% 
% Ms{2}=curr_model;

% RL4s
% curr_model = [];
% curr_model.name = 'a0b3s';
% curr_model.pMin = [1e-6 1e-6 -1 -1 -1];
% curr_model.pMax = [20 1 1 1 1];
% curr_model.pdfs = {beta_sample, alpha_sample, stick_sample, stick_sample, stick_sample}; % Sampling functions for model parameters
% curr_model.pnames = {'beta','alpha+','s1','s2','s3'};
% 
% Ms{2}=curr_model;

% RL4s_hybrid
% curr_model = [];
% curr_model.name = 'a0b3s_hybrid';
% curr_model.pMin = [1e-6 1e-6 -1 -1 -1 1e-6 1e-6 1e-6];
% curr_model.pMax = [20 1 1 1 1 1 1 1];
% curr_model.pdfs = {beta_sample, alpha_sample, stick_sample, stick_sample, stick_sample, lapse_sample, rec_sample, bias_sample}; % Sampling functions for model parameters
% curr_model.pnames = {'beta','alpha+','s1','s2','s3','lapse','rec','bias'};
% 
% Ms{3}=curr_model;

% RL4s_hybrid
% curr_model = [];
% curr_model.name = 'a0b2s_hybrid';
% curr_model.pMin = [1e-6 1e-6 -1 -1 1e-6 1e-6 1e-6];
% curr_model.pMax = [20 1 1 1 1 1 1];
% curr_model.pdfs = {beta_sample, alpha_sample, stick_sample, stick_sample, lapse_sample, rec_sample, bias_sample}; % Sampling functions for model parameters
% curr_model.pnames = {'beta','alpha+','s1','s2','lapse','rec','bias'};
% 
% Ms{4}=curr_model;

% RL4s_hybrid
curr_model = [];
curr_model.name = 'a0b1s_hybrid';
curr_model.pMin = [1e-6 1e-6 -1 1e-6 1e-6 1e-6];
curr_model.pMax = [20 1 1 1 1 1];
curr_model.pdfs = {beta_sample, alpha_sample, stick_sample, lapse_sample, rec_sample, bias_sample}; % Sampling functions for model parameters
curr_model.pnames = {'beta','alpha+','s1','lapse','rec','bias'};

Ms{1}=curr_model;

% RL4s_hybrid
% curr_model = [];
% curr_model.name = 'a01s_hybrid';
% curr_model.pMin = [1e-6 -1 1e-6 1e-6 1e-6];
% curr_model.pMax = [1 1 1 1 1];
% curr_model.pdfs = {alpha_sample, stick_sample, lapse_sample, rec_sample, bias_sample}; % Sampling functions for model parameters
% curr_model.pnames = {'alpha+','s1','lapse','rec','bias'};
% 
% Ms{6}=curr_model;

%% Fit models
for sess = 1
    All_Params = cell(length(Ms), 1);
    All_fits = cell(length(Ms), 1);
    for m = 1:length(Ms)
        fit_model = Ms{m};
        pmin = fit_model.pMin;
        pmax = fit_model.pMax;
        pdfs = fit_model.pdfs;
    
        fitmeasures = cell(length(subjects), 1);
        fitparams = cell(length(subjects), 1);
    
        for k = 1:length(subjects) % no parallel processing
        %parfor k = 1:length(subjects) % parallel processing
            s = subjects(k);
            T = find((X.Animal_ID==s) & X.schedule<3);
            this_data = [X.schedule(T) X.action(T) X.reward1(T)>0];
            this_data = this_data(~isnan(this_data(:,2)),:);
    
            % Sample parameter starting values
            par = zeros(length(pmin), 1);
            for p_ind = 1:length(pmin)
                par(p_ind) = pdfs{p_ind}(0); % 0 is the random seed
            end
    
            % Set starting values of dynamic model parameters to the best fit
            % static model parameters
            if contains(fit_model.name, 'a0b3s_hybrid')
                model_dynamic = fit_model;
                model_static_ind = find(strcmp(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'a0b3s'));
                model_static = Ms{model_static_ind};
                for z = 1:length(model_dynamic.pnames)
                    this_p = model_dynamic.pnames{z};
                    if sum(strcmp(this_p, model_static.pnames)) > 0
                        par(strcmp(this_p, model_dynamic.pnames)) = All_Params{model_static_ind}(k, strcmp(this_p, Ms{model_static_ind}.pnames));
                    end
                end
%                 par(strcmp('lapse', model_dynamic.pnames)) = All_Params{model_static_ind}(k, strcmp('epsilon', model_static.pnames));
%                 par(strcmp('recover', model_dynamic.pnames)) = 1 - All_Params{model_static_ind}(k, strcmp('epsilon', model_static.pnames));
            end
    
            % Define the objective function for optimization
            llhfun = @(p) feval([fit_model.name, '_llh'], p, this_data);
            if sum(strcmp(fit_model.pnames, 'beta')) > 0
                beta_idx = strcmp(fit_model.pnames, 'beta');
                myfitfun = @(p) llhfun(p) + sum((p(beta_idx) - beta_mu).^2 ./ (2*beta_sigma.^2));
            else
                myfitfun = @(p) llhfun(p);
            end
            rng default % For reproducibility
            fmincon_opts = optimoptions(@fmincon, 'Algorithm', 'sqp');
            problem = createOptimProblem('fmincon', 'objective', myfitfun, 'x0', par, 'lb', pmin, 'ub', pmax, 'options', fmincon_opts);
            gs = GlobalSearch;
            [param, llh] = run(gs, problem);
        
            % Calculate fit measures (AIC, BIC, etc.)
            ntrials = size(this_data,1);
            AIC = 2 * llh + 2 * length(param);
            BIC = 2 * llh + log(ntrials) * length(param);
            AIC0 = -2 * log(1/3) * ntrials;
            psr2 = (AIC0 - AIC) / AIC0;
        
            % Store fit measures and parameters for each subject
            fitmeasures{k} = [k llh AIC BIC psr2 AIC0];
            fitparams{k} = param';
        end
    
        % Store fit measures and parameters for each model
        All_Params{m} = cell2mat(fitparams);
        All_fits{m} = cell2mat(fitmeasures);
    end
    
    % Reformat All_fits matrix
    temp = All_fits;
    All_fits = zeros(length(subjects), size(temp{1}, 2), length(Ms));
    for i = 1:length(Ms)
        All_fits(:, :, i) = temp{i};
    end

    save(savedatapath, 'All_fits')
end

%% Plot fit (AIC)
% clear all
% sess = 1;
% load(['FitModels_TSC_sess',num2str(sess)])
% 
% AICs = squeeze(All_fits(:, 3, :));
% mAICs = AICs - repmat(mean(AICs, 2), 1, size(AICs, 2));
% 
% figure('Position', [300 300 900 400])
% subplot(1, 2, 1)
% hold on
% bar(mean(mAICs))
% errorbar(mean(mAICs), std(mAICs) / sqrt(size(mAICs, 1)))
% xticks(1:length(Ms));
% xticklabels(cellfun(@(x) x.name, Ms, 'UniformOutput', false));
% set(gca, 'TickLabelInterpreter', 'none')
% ylabel('\Delta AIC')
% title('Models')
% 
% % static_model_ind = find(contains(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'static'));
% % dynamic_model_ind = find(contains(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'dynamic'));
% static_model_ind = find(strcmp(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'a0b3s'));
% dynamic_model_ind = find(contains(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'a0b3s_hybrid'));
% p = signrank(AICs(:, dynamic_model_ind), AICs(:, static_model_ind), 'tail', 'left');
% 
% subplot(1, 2, 2)
% yline(0, '--')
% hold on
% plot(sort(AICs(:, dynamic_model_ind) - AICs(:, static_model_ind), 'descend'), '.', 'MarkerSize', 15)
% title('a0b3s vs. a0b3s hybrid')
% ylabel('\Delta AIC')
% xlabel('sorted participant')
% set(gca, 'fontsize', 14)
% sgtitle(['AB Sess ' num2str(sess) ' (p=' num2str(p) ')'])
% 
% % Save figures (Fig 4)
% saveas(gcf, ['plots/TSC_fit_sess', num2str(sess), '.png'])
% saveas(gcf, ['plots/TSC_fit_sess', num2str(sess), '.svg'])
% 
% %% Plot fit (BIC)
% clear all
% sess = 1;
% load(['FitModels_TSC_sess',num2str(sess)])
% 
% BICs = squeeze(All_fits(:, 4, :));
% mBICs = BICs - repmat(mean(BICs, 2), 1, size(BICs, 2));
% 
% figure('Position', [300 300 900 400])
% subplot(1, 2, 1)
% hold on
% bar(mean(mBICs))
% errorbar(mean(mBICs), std(mBICs) / sqrt(size(mBICs, 1)))
% xticks(1:length(Ms));
% xticklabels(cellfun(@(x) x.name, Ms, 'UniformOutput', false));
% set(gca, 'TickLabelInterpreter', 'none')
% ylabel('\Delta BIC')
% title('Models')
% 
% % static_model_ind = find(contains(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'static'));
% % dynamic_model_ind = find(contains(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'dynamic'));
% static_model_ind = find(contains(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'a0b1s_hybrid'));
% dynamic_model_ind = find(contains(cellfun(@(x) x.name, Ms, 'UniformOutput', false), 'a0b3s_hybrid'));
% % dynamic_model_ind = 2;
% [p, h, z] = signrank(BICs(:, dynamic_model_ind), BICs(:, static_model_ind));
% 
% subplot(1, 2, 2)
% yline(0, '--')
% hold on
% plot(sort(BICs(:, dynamic_model_ind) - BICs(:, static_model_ind), 'descend'), '.', 'MarkerSize', 15)
% title('a0b1s vs. a0b2s hybrid')
% ylabel('\Delta BIC')
% xlabel('sorted participant')
% set(gca, 'fontsize', 14)
% sgtitle(['AB Sess ' num2str(sess) ' (p=' num2str(p) ', z=' num2str(z.zval) ')'])
% 
% % Save figures (Fig 4)
% saveas(gcf, ['plots/TSC_fit_sess', num2str(sess), '_BIC.png'])
% saveas(gcf, ['plots/TSC_fit_sess', num2str(sess), '_BIC.svg'])