%------------------------------------------
% Date created: 03-10-2023
% @Northwestern Polytechnical University 
% Please contact Lei Du and Jin Zhang(jinzhang@mail.nwpu.edu.cn) for any comments or questions.
% -----------------------------------------
close all; 
clc;
% Load data
% load('Data_SNP_Protein_Imaging.mat');
load('data_using.mat');
sub = snp_position;
snp_name = SNP_ID(sub);
protein_name = plasma_ID;
img_name1 =img_name;
Y{1} = snp(:, sub);
Y{2} = plasma_adj;
Y{3} = img_vbm_adj;
% BL_DX = BL_DX;
% Normalization
X{1}  = getNormalization(Y{1}, 'normalize');
X{2} = getNormalization(Y{2}, 'normalize');
X{3} = getNormalization(Y{3}, 'normalize');
% time_start = cputime;
% Tuned parameters
% mtgpic 
opts.mtgpic.lambda_1 = 0.1; % L1-norm for SNP
% opts.mtgpic.lambda_11 = 1; % FGL 
opts.mtgpic.lambda_21u = 0.1; % L21-norm for SNP
opts.mtgpic.lambda_FGL21 = 1; % FGL21-norm for SNP
opts.mtgpic.lambda_2 = 0.1;  % L1-norm for Protein
opts.mtgpic.lambda_21v = 1; % L21-norm for Protein
opts.mtgpic.lambda_3 = 1;  % L1-norm for Imaging
opts.mtgpic.lambda_GGL = 1;  % GGL-norm for Imaging
opts.mtgpic.lambda = 0.1;% L1-norm for Q
opts.mtgpic.alpha = 1;% alpha
% Kfold Cross validation
n = size(X{1}, 1);
k_fold = 5;
indices = crossvalind('Kfold', n, k_fold);
fprintf('===================================\n');
for k = 1 : k_fold
    fprintf('Current fold: %d\n', k);
    % Split training data and test data
    test = (indices == k);
    train = ~test;
    for i = 1 : numel(X)
        trainData.X{i} = getNormalization(X{i}(train, :), 'normalize');
        testData.X{i}  = getNormalization(X{i}(test, :),'normalize');       
    end
% Training step
    tic;
    [W.mtgpic{k}, u.mtgpic{k}, v.mtgpic{k}, w.mtgpic{k}, Q_mtgpic{k}] = MTGPIC(trainData, opts.mtgpic);
    fprintf('mtgpic: %.3fs\n', toc);
end
% Weights
% mtgpic
u.mtgpic_mean = 0;
v.mtgpic_mean = 0;
w.mtgpic_mean = 0;
for k = 1 : k_fold
    u.mtgpic_mean = u.mtgpic_mean + u.mtgpic{k};
    v.mtgpic_mean = v.mtgpic_mean + v.mtgpic{k};
    w.mtgpic_mean = w.mtgpic_mean + w.mtgpic{k};
end
u.mtgpic_mean = u.mtgpic_mean / k_fold;
v.mtgpic_mean = v.mtgpic_mean / k_fold;
w.mtgpic_mean = w.mtgpic_mean / k_fold;

Q_mtgpic_ave = 0;
for k = 1 : k_fold
    Q_mtgpic_ave = Q_mtgpic_ave + Q_mtgpic{k};
end
Q_mtgpic_ave = Q_mtgpic_ave / k_fold;

% Figures
% ---------------draw snp heatmap-----------------------
figure(1)
ifontsize = 15;
caxis_range = 0.1;
colormap jet;
imagesc(u.mtgpic_mean');
caxis([-1 * caxis_range caxis_range]);

% ---------------draw protein heatmap-----------------------
figure(2)
ifontsize = 15;
caxis_range = 0.1;
colormap jet;
imagesc(v.mtgpic_mean');
caxis([-1 * caxis_range caxis_range]);

% ---------------draw imaging QT heatmap-----------------------
figure(3)
ifontsize = 10;
caxis_range = 0.1;
colormap jet;
imagesc(w.mtgpic_mean');
caxis([-1 * caxis_range caxis_range]);

% ---------------draw snp-protein heatmap-----------------------
figure(4)
ifontsize = 15;
caxis_range = 0.1;
colormap jet;
imagesc(Q_mtgpic_ave);
% set(gca, 'XTick', [], 'YTick', [], 'TickLength', [0 0]);
% set(gca,'Xtick',1:length(Xname), 'xticklabel',Xname, 'XTickLabelRotation',90, 'FontSize',12);
% set(gca,'Ytick',1:length(Yname), 'yticklabel',Yname, 'YTickLabelRotation',0, 'FontSize',6);
caxis([-1 * caxis_range caxis_range]);
