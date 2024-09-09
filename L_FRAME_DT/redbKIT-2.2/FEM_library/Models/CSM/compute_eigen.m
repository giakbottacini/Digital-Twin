function [n, omega, f, T]=compute_eigen(DATA,vertices)

%caso - MENSOLA
    h=max(vertices(2,:))-min(vertices(2,:));
    l=max(vertices(1,:))-min(vertices(1,:));
    %sezione rettangolare
    area=DATA.Thickness*h;
    m=DATA.Density*area;
    J=(DATA.Thickness*h^3)/12;
    alpha=[1.875, 4.694, 7.855];
    %sono richiesti i primi tre modi
    n=3;
    omega=zeros(1,n);
    f=zeros(1,n);
    T=zeros(1,n);
    for i1=1:n
        omega(i1)=(alpha(i1)^2)*sqrt(DATA.Young*J/(m*l^4));
        f(i1)=omega(i1)/(2*pi);
        T(i1)=1/f(i1);
    end
end