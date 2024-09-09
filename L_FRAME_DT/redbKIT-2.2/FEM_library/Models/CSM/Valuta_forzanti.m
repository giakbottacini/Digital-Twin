function Forzanti =  Valuta_forzanti(FE_SPACE, MESH, DATA, param, t)

F_pressure = sparse(MESH.numNodes*MESH.dim, 1);
F_Neumann  = sparse(MESH.numNodes*MESH.dim, 1);

switch MESH.dim
    case 2
        %% Pressure condition
        for k = 1 : MESH.dim
            if ~isempty(MESH.Pressure_side{k})
                [csi,wi]       =  xwgl(FE_SPACE.quad_order, 0, 1);
                [phi]          =  fem_basis(MESH.dim, FE_SPACE.fem, [csi; 0*csi], 1);
                eta            =  1 - csi;
                nqn            =  length(csi);
                
                nof         = length(MESH.Pressure_side{k});
                nbn         = MESH.numBoundaryDof;
                
                Rrows       = zeros(nbn*nof,1);
                Rcoef       = Rrows;
                
                xlt = zeros(nof,nqn); ylt = xlt;
                coord_ref = [eta; csi];
                for j = 1 : 2
                    dof = MESH.boundaries(j,MESH.Pressure_side{k});
                    vtemp = MESH.vertices(1,dof);
                    xlt = xlt + vtemp'*coord_ref(j,:);
                    vtemp = MESH.vertices(2,dof);
                    ylt = ylt + vtemp'*coord_ref(j,:);
                end
                
                pressure = DATA.bcPrex(xlt,ylt,t,param);
                one       = ones(nof,nqn);
                pressure = pressure.*one;

                x    =  MESH.vertices(1,MESH.boundaries(1:2, MESH.Pressure_side{k}));
                y    =  MESH.vertices(2,MESH.boundaries(1:2, MESH.Pressure_side{k}));

                side_length = sqrt((x(2:2:end)-x(1:2:end-1)).^2+(y(2:2:end)-y(1:2:end-1)).^2);

                for l = 1 : nof
                    face = MESH.Pressure_side{k}(l);
                    
                    pressure_loc  = pressure(l,:).*wi;
                    pressure_loc  = pressure_loc(1,:)';
                    
                    Rrows(1+(l-1)*nbn:l*nbn)    = MESH.boundaries(1:nbn,face);
                    Rcoef(1+(l-1)*nbn:l*nbn)    = MESH.Normal_Faces(k,face)*side_length(l)*phi*pressure_loc;
                end
                F_pressure = F_pressure + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,MESH.dim*MESH.numNodes,1).* DATA.Thickness;
                             
            end
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
                F_Neumann = F_Neumann + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,MESH.dim*MESH.numNodes,1).* DATA.Thickness;
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
                for j = 1 : 3
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
                F_Neumann = F_Neumann + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,MESH.dim*MESH.numNodes,1);
                
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
                    for j = 1 : 3
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
                    F_pressure = F_pressure + sparse(Rrows+(k-1)*MESH.numNodes,1,Rcoef,MESH.dim*MESH.numNodes,1);
                end
            end
        end
        
end
    
Forzanti.F_pressure = F_pressure;
Forzanti.F_Neumann  = F_Neumann;
    
end
