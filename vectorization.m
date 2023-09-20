function z = vectorization(Q)
%------------------------------------------
% Date created: 03-10-2023
% @Northwestern Polytechnical University 
% Please contact Lei Du and Jin Zhang(jinzhang@mail.nwpu.edu.cn) for any comments or questions.
% -----------------------------------------
% reshape a d*1 vector w and d*d matrix to d*d+d vector z
    p = size(Q,1);
    d = size(Q,2);
    z = zeros(p*d,1);
    z(1:end) = reshape(Q,p*d,1);
end