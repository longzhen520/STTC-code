function [X,iter,relChgXPath] = Smoothlowrank_TVas( x0, known, rho,T,weigt)
%-----------------------------------------------------------------------------------------------------
% Model: min  sum alpha(t)rank(Z(t))+ rho*(||L||_{TV})
%        s.t. Z=X, X(Omega)=T(Omega), L=Dy,Y=Z
%
% 
%--------------------------------INPUTs----------------------------------------------------------



beta(1)=1/norm(x0(:),'fro');
beta(2)=1/norm(x0(:),'fro');
beta(3)=0.010;

% init. variables 
sizeD=size(x0);
N    = prod(sizeD);

%%  check and initilization
if ~exist('opts','var'); opts=[]; end
if isfield(opts,'maxIter'); maxIter = opts.maxIter; else maxIter =500; end
if isfield(opts,'tol');  tol = opts.tol; else tol = 1e-6; end
% if isfield(opts,'beta'); beta = opts.beta; else beta(1)=0.01; beta(2)=0.001; end
if isfield(opts,'isCont')
    isCont = opts.isCont; contpar = opts.contpar; contfactor = opts.contfactor;
 else
    isCont = true; contpar = 0.95; contfactor = 1+2*(length(known)/N);
end
h    = sizeD(1);
w    = sizeD(2);
d    = sizeD(3);

Z    = zeros(sizeD);
Y    = Z;
X    = x0;
nois  = zeros(sizeD);
L    = zeros(3*N, 1);

v1    = zeros(sizeD);
v2    = zeros(sizeD);
v3    = zeros(3*N, 1);
I=ones(sizeD);
Eny_x   = ( abs(psf2otf([-1; +1], [h,w,d])) ).^2  ;
Eny_y   = ( abs(psf2otf([-1, +1], [h,w,d])) ).^2  ;
Eny_z   = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
Eny_z   =  permute(Eny_z, [3, 1 2]);
% denom1  =  Eny_x + Eny_y ;
denom1  =  (weigt(1)^2)*Eny_x +  (weigt(2)^2)*Eny_y +  (weigt(3)^2)*Eny_z;
% DTD=ifftn(denom1);
 ht = htensor(sizeD);
 alpha = [2,8,4,8,8];
%  alpha = alpha / sum(alpha);
%% main loop

relChgZPath = zeros(maxIter, 1);
relChgXPath = zeros(maxIter, 1);


gamma    = 2;
display  = 1;
tic;
for iter = 1 : maxIter
    
    fprintf('\n*****************************iter: %d ******************************\n', iter');
    Z_pre=Z; X_pre=X; Y_pre=Y;
    
    %-Z subproblem
     alpha = alpha / sum(alpha);

    temp_beta=beta(1)+beta(2);
    tau=alpha/temp_beta;
    temp_Z=(beta(1)*X+v1+beta(2)*Y-v2)/temp_beta;
    [Z,alpha,k] = update_Z(temp_Z, tau, ht);
    
    Z=full(Z);
     %- Y subproblem
     temp_y=beta(2)*Z+v2;
     diffT_L = diffT3( beta(3)*L - v3, sizeD ,weigt);
     numer1   = reshape( diffT_L + temp_y(:), sizeD);
    Y = real( ifftn( fftn(numer1) ./ (beta(3)*denom1 + beta(2)) ) );
    
    % -L subproblem
    diff_y = diff3(Y, sizeD,weigt); 
    L      = softThres( diff_y + v3/beta(3), rho/beta(3));  

    %-X subproblem
       
    X= Z-v1/beta(1);
    X(known)=T(known);


    %- updating multipliers
	v1 = v1 - gamma*beta(1)*(Z - X);
    v2 = v2 - gamma*beta(2)*(Y-Z);
	v3 = v3 - gamma*beta(3)*(L-diff_y);
    
    %- terminating the algorithm and saving some ralted information  
    relChgZ = norm(Z(:)- Z_pre(:),'fro')/norm(Z_pre(:));
    relChgY = norm(Y(:)- Y_pre(:),'fro')/norm(Y_pre(:));
    relChgX = norm(X(:)- X_pre(:),'fro')/norm(X_pre(:));
   
    relChgZPath(iter) = relChgZ;
    relChgXPath(iter) = relChgX;
    relChgYPath(iter) = relChgY;
	if  display
        fprintf('relChgX0:%4.4e,     relChgX1: %4.4e \n       relChgX2: %4.4e\n', relChgZ, relChgX, relChgY);
    end
    
    if (iter> 50) &&  (relChgX < tol ) 
          disp(' !!!stopped by termination rule!!! ');  break;
    end
    
    %- continuation
      if  isCont
            nr1 = norm(Z(:)-X(:), 'fro');
            nr2 = norm(Y(:)-Z(:), 'fro');
            nr3 = norm(L-diff_y, 'fro');
           
            if iter >1 && nr1 > contpar * nr1_pre
                beta(1) = contfactor*beta(1);
            end
            if iter>1 && nr2 > contpar * nr2_pre
                beta(2) = contfactor*beta(2);
            end
            if iter>1 && nr3 > contpar * nr3_pre
                beta(3) = contfactor*beta(3);
            end
            nr1_pre =nr1;    nr2_pre = nr2;  nr3_pre = nr3;
      end    
    
end


%% ouput
out.time     = toc;
out.iter     = iter;
out.nois     = nois;
out.L       = L; 
out.relChgZPath = relChgZPath(1:iter);
out.relChgXPath = relChgXPath(1:iter);
return;

end


