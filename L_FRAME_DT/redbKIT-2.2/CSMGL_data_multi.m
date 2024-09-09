data.datafile      = 'CSMGL_data_multi.m';

%data.model     = 'CSM';
data.Material_Model   = 'LinearSpatial';

%data.plane_analysis='planestrain';
%data.plane_analysis='planestress';

% Number of generated instances
data.how_many_snaps       = 400;
data.how_many_train       = 10000;
data.how_many_test        = 4000;
data.how_many_identify    = 40;

% Source term
data.force{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.force{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Inertia term
data.inertia{1} = @(t, param)(0);
data.inertia{2} = @(t, param)(0);
data.inertia{3} = @(t, param)(0);

% Dirichlet
data.bcDir{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcDir{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Neumann
data.bcNeu{1} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{2} = @(x, y, z, t, param)(0.*x.*y.*z);
data.bcNeu{3} = @(x, y, z, t, param)(0.*x.*y.*z);

% Normal Pressure
data.bcPrex   = @(x, y, z, t, param)((-1.*(z==0.4).*(x>=3.7).*(x<=4).*(y>=3.7).*(y<=4))+0.*x.*y.*z);

% BC flag
data.flag_dirichlet{1}    =  [1];
data.flag_neumann{1}      =  [];
data.flag_pressure{1}     =  [];
data.flag_robin{1}     = [];

data.flag_dirichlet{2}    =  [1];
data.flag_neumann{2}      =  [];
data.flag_pressure{2}     =  [];
data.flag_robin{2}     = [];

data.flag_dirichlet{3}    =  [1];
data.flag_neumann{3}      =  [];
data.flag_pressure{3}     =  [3];
data.flag_robin{3}     = [];

data.u0{1} = @(x, y, z, t, param)(0.*x.*y);
data.u0{2} = @(x, y, z, t, param)(0.*x.*y);
data.u0{3} = @(x, y, z, t, param)(0.*x.*y);

data.du0{1} = @(x, y, z, t, param)(0.*x.*y);
data.du0{2} = @(x, y, z, t, param)(0.*x.*y);
data.du0{3} = @(x, y, z, t, param)(0.*x.*y);

data.density_spatial = @(x,y,z) (2500+0.*x.*y.*z);
%data.Thickness = 0.1;

% Damping da Luca
% sono implementati due possibili sistemi di smorzamento:
% 'None'     - nessun termine di smorzamento è implementato
% 'Rayleigh' - defined by damping coefficient typical for for that
%              structural typology as in Bathe, Finite Element Procedure
%              https://dianafea.com/manuals/d94/Examples/node165.html
%              (vedi dispense Perotti). Nota- è possibile definire un
%              livello di smorzamento per ciascun modo: stiamo considerando
%              che sia lo stesso (pari a data.damp_coeff) per i primi due
%              modi-
data.damp_type  = 'None';
data.damp_coeff = 0.05;
data.n_modi = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.mu=1;
data.lambda=1;

data.muSpatial    = @(x, y, z, t, param)...
                   (((param(5)==0).*...
                    (1 + 0.*x.*y))+...
                    ((param(5)==1).*...
                    (1.*(1-(param(6).*...
                    (x>=0).*(x<=1).*(y>=0).*(y<=0.3).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==2).*...
                    (1.*(1-(param(6).*...
                    (x>=1).*(x<=2).*(y>=0).*(y<=0.3).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==3).*...
                    (1.*(1-(param(6).*...
                    (x>=2).*(x<=3).*(y>=0).*(y<=0.3).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==4).*...
                    (1.*(1-(param(6).*...
                    ((x>=3).*(x<=4).*(y>=0).*(y<=0.3).*(z>=0).*(z<=0.4)+(x>=3.7).*(x<=4).*(y>0.3).*(y<=1).*(z>=0).*(z<=0.4))))))+...
                    ((param(5)==5).*...
                    (1.*(1-(param(6).*...
                    (x>=3.7).*(x<=4).*(y>=1).*(y<=2).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==6).*...
                    (1.*(1-(param(6).*...
                    (x>=3.7).*(x<=4).*(y>=2).*(y<=3).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==7).*...
                    (1.*(1-(param(6).*...
                    (x>=3.7).*(x<=4).*(y>=3).*(y<=4).*(z>=0).*(z<=0.4)))))+...
                    (0.*x.*y)); % * param(1)
data.lambdaSpatial= @(x, y, z, t, param)...
                   (((param(5)==0).*...
                    (1 + 0.*x.*y))+...
                    ((param(5)==1).*...
                    (1.*(1-(param(6).*...
                    (x>=0).*(x<=1).*(y>=0).*(y<=0.3).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==2).*...
                    (1.*(1-(param(6).*...
                    (x>=1).*(x<=2).*(y>=0).*(y<=0.3).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==3).*...
                    (1.*(1-(param(6).*...
                    (x>=2).*(x<=3).*(y>=0).*(y<=0.3).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==4).*...
                    (1.*(1-(param(6).*...
                    ((x>=3).*(x<=4).*(y>=0).*(y<=0.3).*(z>=0).*(z<=0.4)+(x>=3.7).*(x<=4).*(y>0.3).*(y<=1).*(z>=0).*(z<=0.4))))))+...
                    ((param(5)==5).*...
                    (1.*(1-(param(6).*...
                    (x>=3.7).*(x<=4).*(y>=1).*(y<=2).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==6).*...
                    (1.*(1-(param(6).*...
                    (x>=3.7).*(x<=4).*(y>=2).*(y<=3).*(z>=0).*(z<=0.4)))))+...
                    ((param(5)==7).*...
                    (1.*(1-(param(6).*...
                    (x>=3.7).*(x<=4).*(y>=3).*(y<=4).*(z>=0).*(z<=0.4)))))+...
                    (0.*x.*y)); % * param(2)   
                
data.options.LinSolver.solver            = 'backslash';

% Time options
% t0 - initial time
data.time.t0          = 0;
% dt - temporal step
data.time.dt          = 0.005;
% tf_long  - tempo di analisi utilizzato per il record delle instnaces
data.time.tf_long     = 1;
% tf_short  - tempo di analisi utilizzato per il record delle instnaces
data.time.tf_short    = 1;
% t_start  - tempo di inizio record delle instnaces
data.time.t_start     = 0;

gen_alpha_method=0;
if gen_alpha_method
%rho_inf = 1 reduces gamma and theta to the undamped trapezoidal rule
    rho_inf=1;  
    data.time.gamma      = (3-rho_inf)/(2*(1+rho_inf));
    data.time.beta       = 1/((1+rho_inf)^2);
    data.time.alpha_m    = (2*rho_inf-1)/(rho_inf+1);
    data.time.alpha_f    = rho_inf/(rho_inf+1);    
else
    data.time.gamma      = 1/2;
    data.time.beta       = 1/4;
    data.time.alpha_m    = 0;
    data.time.alpha_f    = 0;
end

% data.tol_POD_U: tolleranza sul contenuto energetico degli
% snapshots degli spostamenti
data.tol_POD_U_local    = 1e-3;
data.tol_POD_U          = 1e-3;