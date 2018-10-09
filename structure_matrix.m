function S = structure_matrix(n,d,C,G,noise)

% Kajsa Mollersen, October 8th 2018

% Creates a structure matrix S with entries corresponding to 
% s_{ij} = 1 if P(X_{ij} = 1) > 0.5
% s_{ij} = 0 if P(X_{ij} = 1) < 0.5
% s_{ij} ~ Bernoulli(0.5) otherwise

% Input:    n : number of cells
%           d : number of genes
%           C : 1 x B cellstructure
%           G : 1 x B cellstructure
%           noise : 1 x 2 cellstructure 

% C(b)  : cell entries for block b of s_{ij} = 1
% G(b)  : gene entries for block b of s_{ij} = 1
% noise : cell and gene entries for P(X_{ij} = 1) = 0.5

B = length(C); % The number of blocks

S = false(1,d);

for b = 1 : B
  c = false(n,1);
  c(C{b}) = 1;
  g = false(1,d);
  g(G{b}) = 1;
  S = S | c*g;
end

% I haven't incorporated the noise yet

  
    
    
  



