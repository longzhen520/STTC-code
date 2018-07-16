function [htt,alpha,k] = update_Z(x,tau,ht)
U = cell(1, ht.nr_nodes);
B = cell(1, ht.nr_nodes);
x_ = full(x);
ht_is_leaf = ht.is_leaf;
ht_dims = ht.dims;
k = zeros(1, ht.nr_nodes);
for ii=find(ht_is_leaf)
  Z = matricize(x, ht_dims{ii}, [1:ht_dims{ii}-1, ht_dims{ii}+1:ndims(x)], false);
  [U_,n,S] = Pro2TraceNorm_X(Z, tau(ii));
   alpha(ii)=sum(S(1:n));
    k(ii)=n;
   U{ii} = U_(:,1:k(ii));
  x_ = ttm(x_, U{ii}, ht_dims{ii}, 'h');
end
% Set x to be the reduced tensor x_
x = x_;
ht_lvl = ht.lvl;
for lvl_iter = max(ht_lvl):-1:0
  % Go through all nodes at given level
  for ii=find(ht_lvl == lvl_iter)
    
      
    % Leaves have already been treated
    if(ht_is_leaf(ii))
      continue;
    end
    
    % Matricization of x corresponding to node ii
    x_mat = matricize(x, ht_dims{ii});
    
    % special case root node: matricization is a vector
    if(ii == 1)   
      B_ = x_mat;
      k(ii) = 1;
      alpha(ii)=svd(B_(:));
    else
          [B_,n,S] = Pro2TraceNorm_X( x_mat, tau(ii));
          k(ii)=n;
         alpha(ii)=sum(S(1:n));
              end
   B_ = B_(:, 1:k(ii));
    % Child nodes' indices
    ii_left  = ht.children(ii, 1);
    ii_right = ht.children(ii, 2);
    
    % reshape B{ii} from matrix U_ to a 
    % k(ii) x k(ii_left) x k(ii_right) tensor, 
    B{ii} = dematricize(B_, [k(ii_left), k(ii_right), k(ii)], ...
                        [1 2], 3, false);
    x_mat_ = matricize(x_, ht_dims{ii});
    
    % calculate B{ii}'*x_mat_
    B_x_mat = B_'*x_mat_;
    tsize_red = size(x_); tsize_red(end+1:ndims(ht)) = 1;
    tsize_red(ht_dims{ii_left }(1)) = k(ii);
    tsize_red(ht_dims{ii_right}(1)) = 1;
    
    % Reshape x_mat_ to tensor x_
    x_ = dematricize(B_x_mat, tsize_red, ht_dims{ii});
  end
  
  % Set x to be the reduced tensor x_  
  x = x_;
end

% Call htensor constructor
htt = htensor(ht.children, ht.dim2ind, U, B, true);

end
