function [A_in, F_in, u_D] =  CSM_ApplyBC(A, F, FE_SPACE, MESH, DATA, t, zero_Dirichlet)
%CSM_APPLYBC apply boundary conditions for CSM problem in 2D/3D
%
%   [A_IN, F_IN, U_DIRICHLET] = CSM_APPLYBC(A, F, FE_SPACE, MESH, DATA) given an
%   assembled matrix A, righ-hand side vector F, a FE_SPACE, a MESH data structure and
%   a DATA structure, applies Neumann, Normal Pressure and Dirichlet boundary
%   conditions. It returns the matrix A_IN (matrix A + BCs then restricted
%   to internal dofs), the vector F_IN (vector F + BCs then restricted
%   to internal dofs) and the vector U_DIRICHLET containing the
%   Dirichlet datum evaluated in the Dirichlet dofs.
%
%   [A_IN, F_IN, U_DIRICHLET] = CSM_APPLYBC(A, F, FE_SPACE, MESH, DATA, T) as
%   before, but with the additional input T (time) for time-dependent
%   problems.
%
%   [A_IN, F_IN, U_DIRICHLET] = CSM_APPLYBC(A, F, FE_SPACE, MESH, DATA, T, ZERO_DIRICHLET)
%   If ZERO_DIRICHLET = 1, applies homogeneous Dirichlet boundary
%   conditions (useful for Newton iterations). ZERO_DIRICHLET = 0 by
%   default.

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch> 

if nargin < 6
    t = [];
end

if isempty(A)
    A = sparse(MESH.numNodes*MESH.dim, MESH.numNodes*MESH.dim);
end

if isempty(F)
    F = sparse(MESH.numNodes*MESH.dim, 1);
end

if nargin < 7
    zero_Dirichlet = 0;
end

param = DATA.param;

u_D = [];

switch MESH.dim
    case 2
        %% Pressure condition
        for k = 1 : MESH.dim
            if ~isempty(MESH.Pressure_side{k})
                %resitituisce i punti e i pesi di integrazione di
                %Gauss-Legendre all'interno dell'intervallo 0 e 1
                [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
                %L: raccoglie le funzioni di forma e le loro derivate
                %valutate nei punti di quadratura
                [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; 0*csi], 1);
                %L: coordinate dei punti di integrazione rispetto all'altro
                %estremo dell'elemento (1D) rispetto a cui sono definiti
                eta            =  1 - csi;
                %L: number of quadrature points
                nqn            =  length(csi);
                
                %L: number of nodes belonging to the sides where pressure
                %conditions are applied
                nof         = length(MESH.Pressure_side{k});
                %L: numero di gradi di libertà dell'elemento boundary
                %(compatibili al tipo di elemento finito utilizzato nella
                %mesh)
                nbn         = MESH.numBoundaryDof;
                
                Rrows       = zeros(nbn*nof,1);
                Rcoef       = Rrows;
                
                %L: coordinate totali necessarie per l'applicazione delle
                %condizioni al contorno (numero di righe = gdl boundaries
                %                        numero di colonne = pti quadratura
                xlt = zeros(nof,nqn); ylt = xlt;
                coord_ref = [eta; csi];
                for j = 1 : 2
                    %L: gradi di libertà associati a ciascuno dei lati in
                    %cui sono applicate condizioni al contorno di pressione
                    dof = MESH.boundaries(j,MESH.Pressure_side{k});
                    %L: coordinate dei punti della mesh corrispondenti ai
                    %gradi di libertà dove sono applicate le condizioni al
                    %contorno
                    vtemp = MESH.vertices(1,dof);
                    %L: somma il contributo che gli viene da entrambi i
                    %sistemi di coordinate considerate (xsi e il suo
                    %"complementare" rispetto all'elemento base per cui sono
                    %definiti)
                    xlt = xlt + vtemp'*coord_ref(j,:);
                    vtemp = MESH.vertices(2,dof);
                    ylt = ylt + vtemp'*coord_ref(j,:);
                end
                
                %L: DATA.bcPrex è un function handle (definito all'interno
                %di datafile.m) con il quale si definiscono le condizioni
                %al contorno di pressione in dipendenza delle coordinate
                %dei punti appartenenti al contorno in cui sono state
                %imposte questo tipo di condizioni (oltre che dal tempo e
                %dal tipo di legame costitutivo è stato introdotto per
                %l'integrazione
                pressure = DATA.bcPrex(xlt,ylt,t,param);
                one       = ones(nof,nqn);
                pressure = pressure.*one;
                
                %L: restituisce le coodirnate dei punti relativi alla
                %condizione al contorno k
                x    =  MESH.vertices(1,MESH.boundaries(1:2, MESH.Pressure_side{k}));
                y    =  MESH.vertices(2,MESH.boundaries(1:2, MESH.Pressure_side{k}));
                
                %L: i nodi non hanno una numerazione che segue l'andamento
                %delle coordinate 'x' o 'y'. Questo non cambia niente ai
                %fini del calcolo, neanche nel caso in cui vi sia
                %'applicazione di carichi non lineari. Si osservi che
                %determinare le componenti nodali è equivalente a
                %determinare, per ogni elemento, le reazioni vincolari
                %ipotizzando che questo sia semplicemente appoggiato (i
                %nodi non trasmettono alcuna continuità nel campo della
                %rotazione: diverse considerazioni andranno fatte se la
                %tipologia di elemento finito utilizzato è differente). Per
                %determinare i carichi nodali è necessario effettuare
                %un'integrazione del carico applicato: questo viene fatto
                %con l'applicazione dello schema di Gauss-Legendre.
                side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
                
                %L: l indice utilizzato per la numerazione locale del rhs
                for l = 1 : nof
                    face = MESH.Pressure_side{k}(l);
                    
                    pressure_loc  = pressure(l,:).*wi;
                    pressure_loc  = pressure_loc(1,:)';
                    
                    Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                    Rcoef(1+(l-1)*nbn:l*nbn)    = MESH.Normal_Faces(k,face)*side_length(l)*phi*pressure_loc;
                end
                F = F + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,MESH.dim*MESH.numNodes,1).* DATA.Thickness;
                             
            end
        end
        
        %% Dirichlet condition
        for k = 1 : 2
            if ~isempty(MESH.Dirichlet_dof_c{k})
                x           = MESH.nodes(1,MESH.Dirichlet_dof_c{k});
                y           = MESH.nodes(2,MESH.Dirichlet_dof_c{k});
                u_Dirichlet{k} = DATA.bcDir{k}(x,y,t,param);
                
            else
                u_Dirichlet{k}        = [];
            end
            u_D = [u_D; u_Dirichlet{k}'];
        end
        
        %% Neumann condition
        for k = 1 : MESH.dim
            if ~isempty(MESH.Neumann_side{k})
                [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
                [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; 0*csi], 1);
                eta            =  1 - csi;
                nqn            =  length(csi);
                
                nof         = length(MESH.Neumann_side{k});
                nbn         = MESH.numBoundaryDof;
                
                Rrows       = zeros(nbn*nof,1);
                Rcoef       = Rrows;
                
                xlt = zeros(nof,nqn); ylt = xlt;
                coord_ref = [eta; csi];
                for j = 1 : 2
                    dof = MESH.boundaries(j,MESH.Neumann_side{k});
                    vtemp = MESH.vertices(1,dof);
                    xlt = xlt + vtemp'*coord_ref(j,:);
                    vtemp = MESH.vertices(2,dof);
                    ylt = ylt + vtemp'*coord_ref(j,:);
                end
                
                u_Neumann = DATA.bcNeu{k}(xlt,ylt,t,param);
                one       = ones(nof,nqn);
                u_Neumann = u_Neumann.*one;
                
                x    =  MESH.vertices(1,MESH.boundaries(1:2, MESH.Neumann_side{k}));
                y    =  MESH.vertices(2,MESH.boundaries(1:2, MESH.Neumann_side{k}));
                
                side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);
                
                for l = 1 : nof
                    
                    face = MESH.Neumann_side{k}(l);
                    
                    u_Neumann_loc  = u_Neumann(l,:).*wi;
                    u_Neumann_loc  = u_Neumann_loc(1,:)';
                    
                    Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                    Rcoef(1+(l-1)*nbn:l*nbn)    = side_length(l)*phi*u_Neumann_loc;
                end
                F = F + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,MESH.dim*MESH.numNodes,1).* DATA.Thickness;
            end
        end
        
    case 3
        
        %% Neumann condition
        for k = 1 : MESH.dim
            if ~isempty(MESH.Neumann_side{k})
                
                [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
                csi = quad_points(1,:);
                eta = quad_points(2,:);
                [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
                eta1           =  1-csi-eta;
                nqn            =  length(wi);
                
                nof         = length(MESH.Neumann_side{k});
                nbn         = MESH.numBoundaryDof;
                
                Rrows       = zeros(nbn*nof,1);
                Rcoef       = Rrows;
                
                xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
                coord_ref = [eta1; csi; eta];
                for j = 1 : 2
                    dof = MESH.boundaries(j,MESH.Neumann_side{k});
                    vtemp = MESH.vertices(1,dof);
                    xlt = xlt + vtemp'*coord_ref(j,:);
                    vtemp = MESH.vertices(2,dof);
                    ylt = ylt + vtemp'*coord_ref(j,:);
                    vtemp = MESH.vertices(3,dof);
                    zlt = zlt + vtemp'*coord_ref(j,:);
                end
                
                u_Neumann = DATA.bcNeu{k}(xlt,ylt,zlt,t,param);
                one       = ones(nof,nqn);
                u_Neumann = u_Neumann.*one;
                
                x    =  MESH.vertices(1,MESH.boundaries(1:3, MESH.Neumann_side{k}));
                y    =  MESH.vertices(2,MESH.boundaries(1:3, MESH.Neumann_side{k}));
                z    =  MESH.vertices(3,MESH.boundaries(1:3, MESH.Neumann_side{k}));
                
                areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                  [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                
                for l = 1 : nof
                    
                    area   = 0.5*norm(areav(:,l));
                    detjac = 2*area;
                    
                    face = MESH.Neumann_side{k}(l);
                    
                    u_Neumann_loc  = u_Neumann(l,:).*wi;
                    u_Neumann_loc  = u_Neumann_loc(1,:)';
                    
                    Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                    Rcoef(1+(l-1)*nbn:l*nbn)    = detjac*phi*u_Neumann_loc;
                end
                F = F + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,MESH.dim*MESH.numNodes,1);
                
            end
        end        
        
        %% Pressure condition
        for k = 1 : MESH.dim
            if ~isempty(MESH.Pressure_side{k})
                
                [quad_points, wi] = quadrature(MESH.dim-1, FE_SPACE.quad_order);
                csi = quad_points(1,:);
                eta = quad_points(2,:);
                [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; eta; 0*eta], 1);
                eta1           =  1-csi-eta;
                nqn            =  length(wi);
                
                nbn         = MESH.numBoundaryDof;
                
                for flag = 1 : length(DATA.flag_pressure{k})
                    
                    nof         = length(MESH.Pressure_side_CompFlag{k,flag});
                    
                    Rrows       = zeros(nbn*nof,1);
                    Rcoef       = Rrows;
                    
                    xlt = zeros(nof,nqn); ylt = xlt; zlt = xlt;
                    coord_ref = [eta1; csi; eta];
                    for j = 1 : 2
                        dof = MESH.boundaries(j,MESH.Pressure_side_CompFlag{k,flag});
                        vtemp = MESH.vertices(1,dof);
                        xlt = xlt + vtemp'*coord_ref(j,:);
                        vtemp = MESH.vertices(2,dof);
                        ylt = ylt + vtemp'*coord_ref(j,:);
                        vtemp = MESH.vertices(3,dof);
                        zlt = zlt + vtemp'*coord_ref(j,:);
                    end
                    
                    if length(DATA.bcPrex) == 1
                        pressure = DATA.bcPrex(xlt,ylt,zlt,t,param);
                    else
                        pressure = DATA.bcPrex{DATA.flag_pressure{k}(flag)}(xlt,ylt,zlt,t,param);
                    end
                    
                    one      = ones(nof,nqn);
                    pressure = pressure.*one;
                    
                    x    =  MESH.vertices(1,MESH.boundaries(1:3, MESH.Pressure_side_CompFlag{k,flag}));
                    y    =  MESH.vertices(2,MESH.boundaries(1:3, MESH.Pressure_side_CompFlag{k,flag}));
                    z    =  MESH.vertices(3,MESH.boundaries(1:3, MESH.Pressure_side_CompFlag{k,flag}));
                    
                    areav = cross(  [x(2:3:end)-x(1:3:end);  y(2:3:end)-y(1:3:end);  z(2:3:end)-z(1:3:end)], ...
                        [x(3:3:end)-x(1:3:end);  y(3:3:end)-y(1:3:end);  z(3:3:end)-z(1:3:end)]);
                    
                    for l = 1 : nof
                        
                        area   = 0.5*norm(areav(:,l));
                        detjac = 2*area;
                        
                        face = MESH.Pressure_side_CompFlag{k,flag}(l);
                        
                        pressure_loc  = pressure(l,:).*wi;
                        pressure_loc  = pressure_loc(1,:)';
                        
                        Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                        Rcoef(1+(l-1)*nbn:l*nbn)    = MESH.Normal_Faces(k,face)*detjac*phi*pressure_loc;
                    end
                    F = F + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,MESH.dim*MESH.numNodes,1);
                end
            end
        end
        
        %% Dirichlet condition
        for k = 1 : 3
            if ~isempty(MESH.Dirichlet_dof_c{k})
                
                x           = MESH.nodes(1,MESH.Dirichlet_dof_c{k});
                y           = MESH.nodes(2,MESH.Dirichlet_dof_c{k});
                z           = MESH.nodes(3,MESH.Dirichlet_dof_c{k});
                u_Dirichlet{k} = DATA.bcDir{k}(x,y,z,t,param);
                
            else
                u_Dirichlet{k}        = [];
            end
            u_D = [u_D; u_Dirichlet{k}'];
        end
        
end

u_D  = u_D * (1 - zero_Dirichlet);

if ~isempty( MESH.Dirichlet_dof )
    %L: contributo al rhs relativo ai nodi non vincolati dovuto alle
    %reazioni vincolari dei nodi su cui è applicata una condizione di
    %Dirichet
    F_in = F(MESH.internal_dof) - A(MESH.internal_dof,MESH.Dirichlet_dof)*u_D;
    %L: restrizione (dopo averne considerato il contributo relativo alle
    %reazioni vincolari) della matrice di rigidezza ai soli gdl non
    %vincolati
    A_in = A(MESH.internal_dof,MESH.internal_dof);
    
else
    
    F_in = F(MESH.internal_dof);
    
    A_in = A(MESH.internal_dof,MESH.internal_dof);
    
end

end
