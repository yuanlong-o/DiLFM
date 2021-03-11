
%% this image is a TV denoiser for a given 2D image. anisotropic TV 
%  functionalized
%  last update: 4/11/2020. YZ

% %% parameters
% observ = loadtiff('gaussian_noise_m_0_var_1.tiff');
% observ = double(observ);
% observ = observ / max(observ(:));
% 
% 
% maxIter = 10;
% 
% lambda =1e-2;
% mu = 10;
% tau = 2;
% rho_1 = 10;
% rho_2 = 10;

function x_new = TV_2d_denoiser(observ, lambda, maxIter)

%% parameter
[Nx, Ny] = size(observ);
mu = 10;
tau = 2;
rho_1 = 10;
rho_2 = 10;
%% operator
pad2d = @(x)padarray(x,[Nx,Ny],'both'); 
crop2d = @(x)x(Nx+1:end-Nx,Ny+1:end-Ny); 

vec = @(X)reshape(X,numel(X),1);

% TV operator Psi
Psi = @(x)deal(-diff(x,1,1),-diff(x,1,2)); 

PsiT = @(P1,P2)cat(1,P1(1,:),diff(P1,1,1),-P1(end,:)) + ...
    cat(2,P2(:,1),diff(P2,1,2),-P2(:,end));

PsiTPsi = generate_laplacian_2d(pad2d(zeros(Nx, Ny))); % interesting here
%% initialization
x_new = observ;
[y_1_x, y_1_y] = Psi(x_new);
[z1_x_new, z1_y_new] = Psi(x_new);
z2_new = x_new;
y_2 = x_new;


%% pre-calculation
x_denominator = 1./(2 + rho_1 * PsiTPsi + rho_2);  

%% ADMM main part
for i = 1 : maxIter
    
    % update x
    x_old = x_new;
    x_numerator = 2 * observ + rho_1 * PsiT(z1_x_new,z1_y_new) + rho_2 * z2_new - PsiT(y_1_x, y_1_y) - y_2;
    x_new = crop2d(real(ifft2(x_denominator .* fft2(pad2d(x_numerator)))));
    
    % updage z1
    z1_x_old = z1_x_new;
    z1_y_old = z1_y_new;
%     z1_new = shrinkage(Psi(x_new) + y_1 / rho_1, lambda / rho_1);
    [grad_x, grad_y] = Psi(x_new);
    [z1_x_new, z1_y_new] = shrinkage_2d(grad_x + y_1_x/rho_1, grad_y+y_1_y/rho_1, lambda/rho_1);
    
	s1 = sqrt(norm(vec(z1_x_new - z1_x_old)).^2 + norm(vec(z1_y_new - z1_y_old)).^2); % dual residual
    r1 = sqrt(norm(vec( grad_x- z1_x_new )).^2 + norm(vec( grad_y- z1_y_new )).^2);   % primary residual

    % updage z2
    z2_old = z2_new;
    z2_new = max(y_2/rho_2 + x_new,0);
	s2 = norm(vec(z2_new - z2_old).^2); % dual residual
    r2 = norm(vec(x_new - z2_new).^2);   % primary residual
    
    % update y
    y_1_x = y_1_x + rho_1 * (grad_x - z1_x_new);
    y_1_y = y_1_y + rho_1 * (grad_y - z1_y_new);
    y_2 = y_2 + rho_2 * (x_new - z2_new);
    
    
    % cost info
    data_fidelity = norm(vec(x_new - observ).^2);
    regularizer_penalty = lambda * (sum(abs(vec(grad_x))) + sum(abs(vec(grad_x))));
    
    % update
    if r1 >= mu * s1
        rho_1 = rho_1 * tau;
    elseif s1>= mu * r1
        rho_1 = rho_1 / tau;
    end
    if r2 >= mu * s2
        rho_2 = rho_2 * tau;
    elseif s2>= mu * r2
        rho_2 = rho_2 / tau;
    end    
    
	if mod(i, round(maxIter / 10)) == 0
        fprintf('iter: %i \t total: %i \t data_fidelity: %.2g \t regularizer: %.2g \t prinmary1: %.2g  \t dual1: %.2g \t prinmary2: %.2g  \t dual2: %.2g\n',...
                i, maxIter,data_fidelity,regularizer_penalty, r1, s1, r2, s2)
    end
 
end % end of the while
end
% figure, subplot(1, 2, 1), imshow(observ, []), title('original')
% subplot(1, 2,2), imshow(x_new, []), title('denoised')
% saveastiff(im2uint16(x_new / max(x_new(:))), sprintf('denoised_lambda_%g.tiff', lambda))
%% utility functions
function PsiTPsi = generate_laplacian_2d(lapl) 
    lapl(1) = 4;
    lapl(1,2) = -1;
    lapl(2,1) = -1;
    lapl(1,end) = -1;
    lapl(end,1) = -1;
    PsiTPsi = abs(fft2(lapl));   %Compute power spectrum of laplacian. Note this is 3D version, so might this is the only way?
end


function threshed = shrinkage(x,tau)
    threshed = max(abs(x)-tau,0);
    threshed = threshed.*sign(x);
end


function [varargout] = shrinkage_2d(v, h, tau, varargin)
    if size(v,1) ~= 0  %If no v, h, or d passed in, skip gradient thresholding
        mag = sqrt(cat(1,v,zeros(1,size(v,2),size(v,3))).^2 + ...
            cat(2,h,zeros(size(h,1),1,size(h,3))).^2);
        magt = shrinkage(mag,tau);
        mmult = magt./mag;
        mmult(mag==0) = 0;
        varargout{1} = v.*mmult(1:end-1,:,:);
        varargout{2} = h.*mmult(:,1:end-1,:);
        if ~isempty(varargin)  %Fourth argument is native sparsity
            varargout{4} = shrinkage(varargin{1},tau);
        end
    else
        varargout{1} = shrinkage(varargin{1},tau);
    end

end