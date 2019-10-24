% #############################################################################
%							TP Calculo Numerico
%							2 cuatrimestre 2017
%								
%	Crucianelli, Carla
%	Giardino, Lucas
%	Tracchia, Marcos
%
% #############################################################################
% Definiciones a tener en cuenta
% T             >> Tiempo Total
% k             >> Paso Temporal
% 0 < n < T/k   >> Iteracion Temporal
% N             >> Nro Pasos Espaciales
% -N/2<= i,j <=N/2 
%    i          >> Iteracion de X
%    j          >> Iteracion de Y
% 2L            >> Largo membrana
% h             >> Paso Espacial
% 
% c             >> Constante 
%  
% r := c**2 * k**2 / h**2
clear all;

%DECLARACIONES DE PARAMETROS VARIOS
%-------------------------------------------------------------------------------
N=100;
L1=0.5;
h1=2*L1/N;
Nt=200;
T=2;
c1=1;
Y=X=linspace(-L1,L1,N+2);		  %Discretizo la malla (1)
L2=2.5;
Y2=X2=linspace(-L2,L2,N+2);		%Discretizo la malla 2
t=linspace(0,T,Nt);			      %Discretizo el tiempo
k=T/Nt;
r1=c1**2 * k**2 / h1**2;
a=2;                          %Si a>7 se rompe

aux1=pi/(L1*sqrt(2));	% Las escribo aca afuera
aux2=pi/(2*L1);        % Asi ahorro tiempo en la iteracion
sol_posta=@(x,y,t) cos(aux1*t).*sin(aux2*(x+L1)).*sin(aux2*(y+L1));
cond_init=@(x,y) sin(aux2*(x+L1)).*sin(aux2*(y+L1)); %Condicion inicial (1) (t=0)
cond_init2=@(x,y) -1.*exp(a*(x.**2+y.**2));         %Condicion inicial (2) (t=0)
%-------------------------------------------------------------------------------


% Creo Matriz "pentadiagonal"
%-------------------------------------------------------------------------------
v1(1:N**2,1)=-r1;          %Diags 
v2(1:N-1)=-r1;		          %Diags
v2(N)=0;   		            %Agrego el 0 por chorizear al vector incognitas
diag(1:N**2,1)=1+4*r1;     %Diagonal principal

%Matriz de diagonales:
mat_diags = [repmat(v2,1,N)',v1,diag,v1,repmat(v2,1,N)'] ;
clear v1;
clear v2;
clear diag;

%Matriz del sistema A*u(t+1,x,y)=u(t,x,y)
mat = spdiags(mat_diags,[-N,-1,0,1,N],N**2,N**2);
clear mat_diags;
%-------------------------------------------------------------------------------

%Armo malla a t=0 con la formula que tiran
malla_interior=Y(2:N+1);					%Perdon pero necesito el interior...
malla_interior2=Y2(2:N+1);	
u_0_1=cond_init (malla_interior,malla_interior(:));	%Lo devuelve en una matriz de N*N
u_0_2=cond_init2(malla_interior2,malla_interior2(:));	%Me abuso fuertemente de la simetria

u(:,2,1)=u(:,1,1)=reshape(u_0_1,N**2,1);		%Vuelvo los u0 a un vector y agarro los dos anteriores
u(:,2,2)=u(:,1,2)=reshape(u_0_2,N**2,1);
																							  
error_cosmico(1)=error_cosmico(2)=eps;	%Error inicial = Epsilon de maquina	

%Comienzo la iteracion polenta
%-------------------------------------------------------------------------------										
for n=3:Nt 	%arranco en 3 porque ya tengo el 1,2
	%Calculo el valor posta de la funcion
	la_vera_sol = reshape( sol_posta(malla_interior,malla_interior(:),t(n)) ,N**2,1);
	u(:,n,1)= mat \ (2.*u(:,n-1,1) - u(:,n-2,1));  
	error_cosmico(n) = max(abs(u(:,n,1)-la_vera_sol)); %Agarro el mayor error en este t
endfor
%----------------------------------------------------------------------------------
clear la_vera_sol;

%Segunda Parte tengo que rearmar la matriz
c2=10;
h2=2*L2/N;
r2=c2**2 * k**2 / h2**2;
v1(1:N**2,1)=-r2;         %Diags 
v2(1:N-1)=-r2;		        %Diags
v2(N)=0;   		            %Agrego el 0 por chorizear al vector incognitas
diag(1:N**2,1)=1+4*r2;    %Diagonal principal

%Matriz de diagonales:
mat_diags = [repmat(v2,1,N)',v1,diag,v1,repmat(v2,1,N)'] ;
clear v1;
clear v2;
clear diag;

%Matriz del sistema A*u(t+1,x,y)=u(t,x,y) (2)
mat = spdiags(mat_diags,[-N,-1,0,1,N],N**2,N**2);
clear mat_diags;

%Iteracion del segundo item
for n=3:Nt
  u(:,n,2)= mat \ (2.*u(:,n-1,2) - u(:,n-2,2));
endfor

%-------------------------------------------------------------------------------
clear mat;

%Ploteo los errores maximos para cada t
figure('name','Error Item a)');
xlabel("x");
ylabel("y");
xlabel("tiempo");
ylabel("error");
plot(t,error_cosmico);
clear error_cosmico;

%Ploteo de la peli 1
f1=figure('name','Item a)');
matriz_plot(1:N+2,1:N+2)=0;
matriz_plot(2:N+1,2:N+1)=u_0_1;
fig_1=mesh(X,Y,matriz_plot);
axis([-0.8 0.8 -0.8 0.8 -1 1]);	%Fijo los ejes si no se descontrola todo
xlabel("x");
ylabel("y");
zlabel("u(x,y)");

%Ploteo de la peli 2
f2=figure('name','Item b)');
matriz_plot(2:N+1,2:N+1)=u_0_2;
fig_2=mesh(X2,Y2,matriz_plot);
xlabel("x");
ylabel("y");
zlabel("u(x,y)");
clear u_0_1;
clear u_0_2;

%Arma dos videos continuos simultaneamente
while (ishandle(f1) || ishandle(f2))	%Mientras mantengas algun grafico abierto
	for n=1:Nt
		if ishandle(f1)				%Si algun grafico se cierra que el otro no pare
			matriz_plot(2:N+1,2:N+1)=reshape(u(:,n,1),N,N);
			set(fig_1,"ZData",matriz_plot);
			matriz_plot(:,:)=0;
		endif
		if ishandle(f2)				%Si algun grafico se cierra que el otro no pare
			matriz_plot(2:N+1,2:N+1)=reshape(u(:,n,2),N,N);
			set(fig_2,"ZData",matriz_plot);
			matriz_plot(:,:)=0;
		endif
		pause(0.02);
		if !(ishandle(f1) || ishandle(f2))
			break;
		endif	
	endfor
endwhile
clear all;