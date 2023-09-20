function [W, U, V, w, Qout] =  MTGPIC(Data, opts)
%------------------------------------------
% Date created: 03-10-2023
% @Northwestern Polytechnical University 
% Please contact Lei Du and Jin Zhang(jinzhang@mail.nwpu.edu.cn) for any comments or questions.
% -----------------------------------------

% data
X = Data.X{1}; p = size(X, 2);
Y = Data.X{2}; d = size(Y, 2);
Z = Data.X{3}; q = size(Z, 2);

% pre-calculate
% XX = X' * X;
YY = Y' * Y; ZZ = Z' * Z;
XY = X' * Y; YZ = Y' * Z; XZ = X' * Z;

% set parameters
lambda_1 = opts.lambda_1;
lambda_21u = opts.lambda_21u;
lambda_FGL21 = opts.lambda_FGL21;

lambda_21v = opts.lambda_21v;
lambda_2 = opts.lambda_2;

lambda_3 = opts.lambda_3;
lambda_GGL = opts.lambda_GGL;

lambda = opts.lambda;
alpha = opts.alpha;
% initialize u, v and w
u1 = ones(p, 1); u1= u1 / norm(X * u1);
u2 = ones(p, 1); u2= u2 / norm(X * u2);
v1 = ones(d, 1); v1 = v1 / norm(Y * v1);
v2 = ones(d, 1); v2 = v2 / norm(Y * v2);
w = ones(q, 1); w = w / norm(Z * w);

% initialize loss weights
Xu1 = X * u1; Xu2 = X * u2;
Yv1 = Y * v1; Yv2 = Y * v2;
Zw = Z * w;

% iteration
iter = 0; maxIter = 50;
tall = inf; tol = 1e-5;

Iinter = [];
I = [];
% creat interactions
for i = 1:p
    I = X(:,i).*Y;
    Iinter = [Iinter,I];
end

Q_ini = rand(p,d);
Qvec = vectorization(Q_ini);
Qvec = Qvec./norm(Qvec);

% get the structure
block = 50;
Xall = X;
%% solve the MT-GPIC problems along the sequence of parameter values
while (iter < maxIter && tall > tol)
    iter = iter + 1;
    
    % update u
    nblock = ceil(p/block);
    u1t = [];
    u2t = [];
    QB = [];
    su1 = 0;
    su2 = 0;
    s2 = 0;
    u1_old = u1;
    u2_old = u2;
    for iu = 1:nblock
        if iu*block <= p
            X = Xall(:,1+(iu-1)*block:iu*block);
            sub_u1 = u1(1+(iu-1)*block:iu*block);
            sub_u2 = u2(1+(iu-1)*block:iu*block);
            Iinterblock = Iinter(:,1+(iu-1)*block*d:iu*block*d);
            Qvecblock = Qvec(1+(iu-1)*block*d:iu*block*d);
        else
            X = Xall(:,1+(iu-1)*block:end);
            sub_u1 = u1(1+(iu-1)*block:end);
            sub_u2 = u2(1+(iu-1)*block:end);
            Iinterblock = Iinter(:,1+(iu-1)*block*d:end);
            Qvecblock = Qvec(1+(iu-1)*block*d:end);
        end
        XX = X'*X;
        % solve u1 u2
        % updata Du
        Du21 = updateDs(sub_u1, sub_u2);%L21
        DuFGL21 = updateD_FGL21(sub_u1,sub_u2);
        D1 = updateD(sub_u2);
        F1 = alpha * XX  + lambda_1 * D1 + lambda_21u * Du21 + lambda_FGL21 * DuFGL21;
        b1 = X' * Yv2 ;
        sub_u2 = F1\b1;
        u2t = [u2t; sub_u2];
        su2 = su2 + sub_u2' * XX * sub_u2;
        % update Du
        D1 = updateD(sub_u1);
        F1 = alpha * XX  + lambda_1 * D1 + lambda_21u * Du21 + lambda_FGL21 * DuFGL21;
        b1 = X' * Zw - X' * Yv1 - X' * Iinter * Qvec;
        sub_u1 = F1\b1;
        u1t = [u1t; sub_u1];
        su1 = su1 + sub_u1' * XX * sub_u1;

    % update Q
    D3 = updateD(Qvecblock);
    F3 = Iinterblock' * Iinterblock + lambda * D3;
    b3 = Iinterblock' * (Zw - Xall * u1 -Yv1);
    Qvecblock = F3 \ b3;
    QB = [QB; Qvecblock];
    s2 = s2 + Qvecblock' * Qvecblock;
    end
    % scale u1 u2 Q
    u1 = u1t./su1;
    Xu1 = Xall * u1;
    u2 = u2t./su2;
    Xu2 = Xall * u2;
    Qvec = QB./norm(QB);
    
    % solve v
    Dv21 = updateDs(v1, v2);% L21
    % update v
    v1_old = v1;
    v2_old = v2;
    D2 = updateD(v2);
    F2 = alpha * YY + lambda_2 * D2 + lambda_21v * Dv21;
    b2 = alpha * Y'* Xu2;
    v2 = F2 \ b2;
    v2 = v2 / norm(Y * v2);
    Yv2 = Y * v2;
    D2 = updateD(v1);
    F2 =  YY + lambda_2 * D2+ lambda_21v*Dv21;
    b2 =  Y' * Zw - Y'* Xu1 - Y' * Iinter * Qvec;
    v1 = F2 \ b2;    v1 = v1 / norm(Y * v1);
    Yv1 = Y * v1;
   
    % update w
    w_old = w;
    % ----------------------------------------
    D3 = updateD(w);
    % update diagnal matrix Dw
    DGGL = updateD(w, 'GGL');    
    D2 = diag(D2);
    F3 =  ZZ + lambda_3 * D3 + lambda_GGL * DGGL;
    b3 =  Z' * Xu1 +  Z' * Yv1 + Z' * Iinter * Qvec;
    w = F3 \ b3;
    w = w / norm(Z * w);
    Zw = Z * w;
    % stopping condition
    % ----------------------------------------
    tu1 = max(abs(u1 - u1_old));
    tv1 = max(abs(v1 - v1_old));
    tu2 = max(abs(u2 - u2_old));
    tv2 = max(abs(v2 - v2_old));
    tw = max(abs(w - w_old));
    tall = max([tu1,tv1,tu2,tv2,tw]);
end
U = [u1 u2];
V = [v1 v2];
W{1} = U; W{2} = V; W{3} = w; 
Qout = de_vectorization(Qvec,p,d);

end

function D = updateD(w, type)

        if nargin == 1
            % for L1-norm
            d = 1 ./ sqrt(w .^ 2 + eps);
        elseif strcmpi(type, 'FGL')
            % for FGL-norm
            [n_features, ~] = size(w);
            structure = updateGraph(n_features, 'FGL');
            Gp = 1 ./ sqrt(structure * (w .^ 2) + eps);
            d = [Gp(1), sum(reshape(Gp(2 : end - 1), 2, [])), Gp(end)];
        elseif strcmpi(type, 'GGL')
            % for GGL-norm
            [n_features, ~] = size(w);
            structure = updateGraph(n_features, 'GGL');
            Gp = 1 ./ sqrt(structure * (w .^ 2) + eps);
            d = sum(reshape(Gp, n_features - 1, []));
        else
            error('Error type.');
        end

        D = diag(d);

    function E = updateGraph(n, type)

        if strcmpi(type, 'FGL')
            E = zeros(2 * (n - 1), n);
            for i = 1 : n - 1
                j = i + 1;
                E(2 * i - 1, i) = 1;
                E(2 * i - 1, j) = 1;
                E(2 * i, i) = 1;
                E(2 * i, j) = 1;
            end
        elseif strcmpi(type, 'GGL')
            num = 0;
            E = zeros(n * (n - 1), n);
            for i = 1 : n
                for j = 1 : n
                    if i ~= j
                        num = num + 1;
                        E(num, i) = 1;
                        E(num, j) = 1;
                    end
                end
            end
        else
            error('Error type.');
        end
        end
end

function Dout = updateD_FGL21(u1, u2)

    ulen = length(u1);
    for i = 1:ulen

    if i == 1
            d(i) = sqrt(u1(i).^2+u2(i).^2+u1(i+1).^2+u2(i+1).^2+eps);
    d(i) = 0.5 ./ d(i);

    elseif i==ulen
            d(i) = sqrt(u1(i-1).^2+u2(i-1).^2+u1(i).^2+u2(i).^2+eps);
    d(i) = 0.5 ./ d(i);

    else 
            d(i) = 0.5./(sqrt(u1(i-1).^2+u2(i-1).^2+u1(i).^2+u2(i).^2+eps))+0.5./(sqrt(u1(i).^2+u2(i).^2+u1(i+1).^2+u2(i+1).^2+eps));
    end

    D = d;
    Dout = diag(D);
    end

end

function Dout = updateDs(v1, v2, v3, v4)

    vlen = length(v1);

    if nargin < 4
        for i = 1:vlen
            d(i) = sqrt(v1(i).^2+v2(i).^2+eps);
        end
    else
        for i = 1:vlen
            d(i) = sqrt(v1(i).^2+v2(i).^2+v3(i).^2+v4(i).^2+eps);
        end
    end
    D = 0.5 ./ d;
    Dout = diag(D);
    % pen_Value = sum(d);
    
end
