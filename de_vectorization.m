function [Q] = de_vectorization(z,p,d)
%------------------------------------------
% Date created: 03-10-2023
% @Northwestern Polytechnical University 
% Please contact Lei Du and Jin Zhang(jinzhang@mail.nwpu.edu.cn) for any comments or questions.
% -----------------------------------------
% reshape a d*d+d vector z to a d*1 vector w and d*d matrix Q
Q = reshape(z(1:end,1),p,d);
end