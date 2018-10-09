function p = cell_gene_effect(distr, param, block, p_const,m)

% Kajsa Mollersen (kajsa.mollersen@uit.no) October 9th 2018

% Input:    distr - which distr, string
%           param - vector of parameters
%           block - these rows/columns will have the same distribution
%           p_constant - (0.5,1) noise
%           m - length of cell/gene effect vector

if strcmp(distr,'Normal')
  mu = param(1);
  sigma = param(2);
  
  p_distr = random(distr,mu,sigma,ceil(m/block),1);
  
  p_distr(p_distr < p_const - 1) = p_const - 1;
  p_distr(p_distr > 1 - p_const) = 1 - p_const;
  
  p = zeros(m,1);
  for i = 1: m
    p(i) = p_distr(ceil(i/block));
  end
end

if strcmp(distr,'Beta')
  % Not ready yet
  
end

if strcmp(distr,'Uniform')
  % Not ready yet
end

