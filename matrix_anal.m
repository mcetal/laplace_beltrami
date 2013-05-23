close all; clear;
k = 1;
load amat_stereo.dat;
load amat_sphere.dat;
n = size(amat_stereo);
NX = sqrt(n(1)); NY = sqrt(n(1));
np = (NX-k)/k;
mat_st = zeros(NX,NY); mat_sp = zeros(NX,NY);
mat_st(:) = amat_stereo(:); mat_st = mat_st';
mat_sp(:) = amat_sphere(:); mat_sp = mat_sp';
disp(['Condition number of stereo system matrix = ',num2str(cond(mat_st))])
prec = mat_st;
prec(1:k*np,1:k*np) = 0.5*eye(k*np);
prec_mat = inv(prec)*mat_st;
cond_prec = cond(prec_mat);
disp(['Condition number of preconditioned stereo system = ',num2str(cond_prec)])
disp(['Condition number of sphere system matrix = ',num2str(cond(mat_sp))])
prec = mat_sp;
prec(1:k*np,1:k*np) = 0.5*eye(k*np);
prec_mat = inv(prec)*mat_sp;
cond_prec = cond(prec_mat);
disp(['Condition number of preconditioned sphere system = ',num2str(cond_prec)])
