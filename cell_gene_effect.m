function p = cell_gene_effect(distr, param, block, p_const,n,d,C, cell_effect, cell_genes)

% Kajsa Mollersen (kajsa.mollersen@uit.no) October 9th 2018

% Input:    distr - which distr, string
%           param - vector of parameters
%           block - these rows/columns will have the same distribution
%           p_constant - (0.5,1) noise
%           m - length of cell/gene effect vector

if strcmp(distr,'Normal')
  sigma = param(2); 
  p = zeros(n,d);
  dg = cell_genes;
  for b = 1: length(cell_effect)
    rng('default')
    mu = param(1)+cell_effect(b);
  
    p_distr = random(distr,mu,sigma,ceil(length(dg)/block),1);
  
    p_distr(p_distr < p_const - 1) = p_const - 1;
    p_distr(p_distr > 1 - p_const) = 1 - p_const;
  
  
    for j = 1: length(dg)
      p(C{b},dg(j)) = p_distr(ceil(j/block));
    end
  end
  
  rng('default')
  mu = param(1);
  dng = setdiff(1:d, dg);
  whos dng
  p_distr = random(distr,mu,sigma,ceil(length(dng)/block),1);
  whos p_distr
  
  p_distr(p_distr < p_const - 1) = p_const - 1;
  p_distr(p_distr > 1 - p_const) = 1 - p_const;
  
  for j = 1: length(dng)
    p(:,dng(j)) = p_distr(ceil(j/block));
  end
  
end

if strcmp(distr,'Beta')
  % Not ready yet
  
end

if strcmp(distr,'Uniform')
  % Not ready yet
end

