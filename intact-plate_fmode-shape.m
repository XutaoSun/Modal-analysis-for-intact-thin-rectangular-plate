clear all
clc
format long
n=80;
a=0.8; 
b=2*a; %a, b are the plate dimensions in x and y directions
h=0.008; %thickness
al=a/b;
nu=0.33; %the Poisson's ratio
E=69e9; %Young's modulus
Rho=2700; %density
D=E*h^3/(12*(1-nu^2)); %bending rigidity of the plate
%penalty terms for edges of each plate segment kp1 kr1(left edge) kp2 kr2(right edge) kp3 kr3(lower edge) kp4 kr4(upper edge)
kp1=0; %1e9; %left edge penalty
kp2=1e10; %upper edge
kp3=0; %1e9; %right edge
kp4=1e10; %lower edge
kr1=0; %left edge
kr2=1e10; %upper edge
kr3=0; %right edge
kr4=1e10; %lower edge
%%%%%%%%%%%%%%%%%%%%%% define integrals %%%%%%%%%%%%%%%%%%%%
% Product of intergrals of admissible functions Emi00 stands for the product of two zeroth derivatives, Emi02 stands for the product of the zeroth derivative and second derivative (the order should be consistent with arguments in brackets)
Emi00(1,1)=1;
Emi00(1,2)=1/2;
Emi00(2,1)=1/2;
Emi00(3,1)=1/3;
Emi00(1,3)=1/3;
Emi00(2,2)=1/3;
Emi00(2,3)=1/4;
Emi00(3,2)=1/4;
Emi00(3,3)=1/5;
for i=4:n
    Emi00(2,i)=-(1-cos((i-3)*pi))/((i-3)^2*pi^2);
    Emi00(i,2)=-(1-cos((i-3)*pi))/((i-3)^2*pi^2);
    Emi00(3,i)=(2*cos((i-3)*pi))/((i-3)^2*pi^2);
    Emi00(i,3)=(2*cos((i-3)*pi))/((i-3)^2*pi^2);
end
for i=4:n
    Emi00(i,i)=1/2;
end
%%%%
Emi22(3,3)=4;
for i=4:n
    Emi22(i,i)=0.5*(i-3)^4*pi^4;
end
%%%%
Emi02(1,3)=2;
Emi02(2,3)=1;
Emi02(3,3)=2/3;
for i=4:n
    Emi02(2,i)=1-cos((i-3)*pi);
    Emi02(3,i)=-2*cos((i-3)*pi);
    Emi02(i,i)=-0.5*(i-3)^2*pi^2;
end
%%%%
Emi11(2,2)=1;
Emi11(2,3)=1;
Emi11(3,2)=1;
Emi11(3,3)=4/3;
for i=4:n
    Emi11(2,i)=cos((i-3)*pi)-1;
    Emi11(i,2)=cos((i-3)*pi)-1;
    Emi11(3,i)=2*cos((i-3)*pi);
    Emi11(i,3)=2*cos((i-3)*pi);
    Emi11(i,i)=0.5*(i-3)^2*pi^2;
end
%%%%
Emi20=Emi02';
Fni00=Emi00;
Fni20=Emi20;
Fni02=Emi02;
Fni11=Emi11;
Fni22=Emi22;
%%%%
%%%%%%%%%%%%%%%%%% mass matrix %%%%%%%%%%%%%%%%%%
for m=1:n 
	for i=1:n 
		u=i+(m-1)*(n); 
		for nn=1:n
			for j=1:n
				v=j+(nn-1)*(n);
				M(u,v)=Emi00(m,nn)*Fni00(i,j)*a*b*Rho*h;
			end
		end
	end
end
%KP1是板的左边界，KP2是板的右边界，KP3是下边界，KP4是上边界，KC1是左边界转角，KC2是右边界转角，KC3是下边界转角，KC4是上边界转角。
for i=1:n
	KP1(i)=1;
	KP3(i)=1;
	KC1(i)=0;
	KC2(i)=0;
end
KP1(2)=0;
KP1(3)=0;
KP3(2)=0;
KP3(3)=0;
KC1(2)=1;
KC2(1)=0;
KC2(2)=1;
KC2(3)=2;
KC3=KC1;
KC4=KC2;
KP2(1)=1;
KP4(1)=1;
KP2(2)=1;
KP4(2)=1;
KP2(3)=1;
KP4(3)=1;
for i=4:n
	KP2(i)=cos((i-3)*pi);
	KP4(i)=KP2(i);
end
%%%%%%%%%%%%%%%%%%% stiffness matrix %%%%%%%%%%%%%%%%%%%
K=zeros(n^2,n^2);
for m=1:n
	for i=1:n
		u=i+(m-1)*n;
		for nn=1:n
			for j=1:n
				v=j+(nn-1)*n;
				%for K11
				K(u,v)=K(u,v)+(Emi22(m,nn)*Fni00(i,j)+al^4*Emi00(m,nn)*Fni22(i,j)+nu*al^2*Emi02(m,nn)*Fni20(i,j)+nu*al^2*Emi20(m,nn)*Fni02(i,j)+2*(1-nu)*al^2*Emi11(m,nn)*Fni11(i,j))*D*b/a^3; %the segment 1 itself
				K(u,v)=K(u,v)+kp2*KP4(i)*KP4(j)*Emi00(m,nn)*a; %the translational boundary condition for the upper edge of segment 1
				K(u,v)=K(u,v)+kr2*KC4(i)/b*KC4(j)/b*Emi00(m,nn)*a; %the rotational boundary condition for the upper edge of segment 1
				K(u,v)=K(u,v)+kp4*KP3(i)*KP3(j)*Emi00(m,nn)*a; %the translational boundary condition for the lower edge of segment 1
				K(u,v)=K(u,v)+kr4*KC3(i)/b*KC3(j)/b*Emi00(m,nn)*a; %the rotational boundary condition for the lower edge of segment 1
				K(u,v)=K(u,v)+kp1*KP1(m)*KP1(nn)*Emi00(i,j)*b; %the translational boundary condition for the left edge of segment 1
				K(u,v)=K(u,v)+kr1*KC1(m)/a*KC1(nn)/a*Emi00(i,j)*b; %the rotational boundary condition for the left edge of segment 1
				K(u,v)=K(u,v)+kp3*KP2(m)*KP2(nn)*Emi00(i,j)*b; %the translational boundary condition for the right edge of segment 1
				K(u,v)=K(u,v)+kr3*KC2(m)/a*KC2(nn)/a*Emi00(i,j)*b; %the rotational boundary condition for the right edge of segment 1
			end
		end
	end
end
[X,eigenvalues] = eig(K,M);
eigenvalues = diag(eigenvalues);
Omega=eigenvalues.^0.5*a^2*(Rho*h/D)^0.5;
Omega=real(Omega);
fHz=eigenvalues.^0.5/(2*pi);
fHz(1:10)
%%%%%%%%%%%%%%%%%% building the eigenvectors %%%%%%%%%%%%%%
x_nodes=linspace(0,a,31);
y_nodes=linspace(0,b,61); %coordinates of nodes in the half segment
for MODE_NUMBER=1:6
	EIGENVECTOR=real(X(:,MODE_NUMBER))';
	for i=1:size(x_nodes,2)
		x=x_nodes(i);
		for j=1:size(y_nodes,2)
			y=y_nodes(j);
			SET_X(1)=1;
			SET_X(2)=x/a;
			SET_X(3)=(x/a)^2;
			SET_Y(1)=1;
			SET_Y(2)=y/b;
			SET_Y(3)=(y/b)^2;
			for k=1:n-3
				SET_X(k+3)=cos(k*pi*x/a);
				SET_Y(k+3)=cos(k*pi*y/b);
			end
			counter=0;
			for p=1:n
				for q=1:n
					counter=counter+1;
					EVALUATED_FUNCTIONS(counter)=SET_X(p)*SET_Y(q);
				end
			end
			EVALUATED_MODE(i,j)=-dot(EVALUATED_FUNCTIONS,EIGENVECTOR); %这一个点的模态
		end
	end
	figure
	Modeshape=surf(y_nodes,x_nodes,EVALUATED_MODE);
	axis off;
	zdir=[0,90];
	rotate(Modeshape,zdir,180)
end
