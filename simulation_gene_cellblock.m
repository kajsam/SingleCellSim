function [X, p] = simulation_gene_cellblock(fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no), October 9th 2018

% Requires: structure_matrix.m, cell_gene_effect.m, observation_matrix.m

% Creates a structure matrix, the cell and/or gene effects, and finally the
% observed binary gene expression matrix.

%% For the figures
% The size
n = 1000; 
d = 5000;

rng('default')
clf(figure(fig_nr))

% IM = ones(n,d); This is just confusing
xlab = 'Genes';
ylab = 'Cells';

%% The structure matrix

C = cell(1,4);
C{1} = 1: n;
C{2} = 1: ceil(n/3);
C{3} = ceil(n/3)+1:ceil(n/3)+ceil(n/4);
C{4} = setdiff(C{1},union(C{2},C{3}));

G = cell(1,4);
G{1} = 1: floor(0.2*d);
G{2} = floor(0.2*d)+1: floor(0.3*d);
G{3} = floor(0.3*d)+1: ceil(0.45*d);
G{4} = floor(0.45*d)+1: ceil(0.55*d);

S = structure_matrix(n,d,C,G);

if fig_nr
  figure(fig_nr), subplot(1,4,1), colormap(gray)
  imagesc(S)
  title('S')
  set(gca,'xaxisLocation','top')
  xlabel(xlab)
  ylabel(ylab)
  drawnow
end

%% The Bernoulli parameter matrix. 
% Each entry gives defines the Bernoulli distribution from which x_{ij} is drawn. 

p_const = 0.65; % An overall probability of X_{ij} = s_{ij}

block = 5;
distr = 'Normal';
param = [0 0.1];

p = cell_gene_effect(distr, param, block, p_const, n,d,C(2:4), [0.1 -0.1 0.05], [401:3500]);


figure(fig_nr+1), imagesc(p), colormap(gray), title('True gene effect')
% figure(fig_nr), subplot(1,4,2), 

figure(fig_nr), subplot(1,4,2) , histogram(p(:),100), title(strcat('Truncated ',{' '}, distr))

%% The simulated observation matrix
[X, Pi] = observation_matrix(S, p_const, p);

figure(fig_nr), subplot(1,4,3), imagesc(Pi, [0 1]), colormap(gray)
title(strcat('Bernouilli \pi = ',num2str(p_const)))
xlabel(xlab)
ylabel(ylab)
drawnow

figure(fig_nr), subplot(1,4,4), colormap(gray)
imagesc(X, [0 1])
title('X')
xlabel(xlab)
ylabel(ylab)




