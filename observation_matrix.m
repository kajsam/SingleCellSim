function [X, Pi] = observation_matrix(S, p_const, p)

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
% IM = ones(n,d); % Confusing

% No cell or gene effects: the probability is a scalar
if all(size(p) == [1 1])
  p = p.*zeros(n,d);
  effect = 'random effect';
elseif all(size(p) == [n d]) % Cell and gene effect: the probability is a nxd matrix
  effect = '+ gene and cell effect';
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
else
  disp('Wrong dimension of p')
  return
end

Pi = Pi + p;
find(Pi>1)
find(Pi<0)
figure, imagesc(Pi), colormap(gray)

% uPi = unique(Pi);
% luPi = length(uPi);
% if luPi < n*d
%   idx = cell(1,luPi);
%   lidx = zeros(1,luPi);
%   for i = 1: luPi
%     idx{i} = find(Pi==uPi(i));
%     lidx(i) = length(idx{i});
%   end
% end
    
X = false(n*d,1);
% for i = 1: luPi
%   X(idx{i}) = random('Binomial', 1, uPi(i)*ones(1,lidx(i))); % cell i, gene j
% end
for i = 1: n*d
  X(i) = random('Binomial', 1, Pi(i)); % cell i, gene j
end
X = vec2mat(X,n);
X = X';

