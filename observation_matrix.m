function X = observation_matrix(S, p_const, p, fig_nr, xlab, ylab)

% Kajsa Mollersen (kajsa.mollersen@uit.no) October 9th 2018

% This function simulates a two-class binary single cell RNA expression 
% matrix.

% Input:    
% S:        The structure matrix 
% p_const:  The overall uncertainty
% p:        The cell/gene effect

rng('default') % for reproducibility


%% Generate the probability matrix P of size (n,d), where each entry is the 
% probabilty of success (x_ij = 1)
Pi = (1 - S) + (2.*S - 1).*p_const;
[n, d] = size(Pi);
IM = ones(n,d);

% No cell or gene effects: the probability is a scalar
if all(size(p) == [1 1])
  p = p.*zeros(n,d);
  effect = 'random effect';
elseif any(size(p) == n) % Cell effect: the probability is a nx1 vector
  if size(p,2) == n
    p = p';
  end
  effect = ' + cell effect';
  p = p*ones(1,d);
elseif  any(size(p) == d) % Gene effect: the probability is a 1xd vector
  if size(p,1) == d
    p = p';
  end
  p = ones(n,1)*p;
  effect = ' + gene effect';
elseif all(size(p) == [n d]) % Cell and gene effect: the probability is a nxd matrix
  effect = '+ gene and cell effect';
else
  'Wrong dimension of p'
  return
end

Pi = Pi + p;
figure(fig_nr), subplot(1,3,1), colormap(gray)
imagesc((IM - Pi), [0 1])
title(strcat('Bernouilli \pi = ',num2str(p_const), effect))
xlabel(xlab)
ylabel(ylab)
drawnow

uPi = unique(Pi);
luPi = length(uPi);
if luPi < n*d
  idx = cell(1,luPi);
  lidx = zeros(1,luPi);
  for i = 1: luPi
    idx{i} = find(Pi==uPi(i));
    lidx(i) = length(idx{i});
  end
end
    
X = false(n*d,1);
for i = 1: luPi
  X(idx{i}) = random('Binomial', 1, uPi(i)*ones(1,lidx(i))); % cell i, gene j
end
X = vec2mat(X,n);
X = X';

subplot(1,3,2), colormap(gray)
imagesc((IM - X), [0 1])
title('X')
xlabel(xlab)
ylabel(ylab)
