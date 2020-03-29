% not acheived the objective
weights=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];          %7 3 6 
cy=[0 0 1 0 -1 1 1 -1 -1];          % 4 1 2 
q=9;                                %8 5 9
nx=20;
ny=20;
u0=0.05;
omega=0.56; % taken from rattans code
ts =30000;

rho0 = 1; 
feq = zeros(q);
f= zeros(q,nx,ny);
fs = zeros(q,nx,ny);

ux = zeros(nx,ny);
uy = zeros(nx,ny);
rho = ones(nx,ny);

%initialise the distributive function
for x=1:nx
    for y = 1:ny
        for k = 1:q
        feq(k)=weights(k)*rho0;
        f(k,x,y)=feq(k);
        end
    end
end

%main for ts time step

for t=1:ts
    for x = 1:nx
        for y =1:ny
         %update rho ux uy
          dense=0;
          vx=0;
          vy=0;
          for k=1:q
              dense=dense+f(k,x,y);
              vx=vx+cx(k)*f(k,x,y);
              vy=vy+cy(k)*f(k,x,y);
          end
          
          if(dense~=0)
          rho(x,y)=dense;
          ux(x,y)=vx/dense;
          uy(x,y)=vy/dense;
          else
          rho(x,y)=0;
          ux(x,y)=0;
          uy(x,y)=0;    
          end
           for k=1:q
                %feq update
                feq(k) = weights(k)*rho(x,y)*(1+3*(cx(k)*ux(x,y)+cy(k)*uy(x,y))+4.5*(cx(k)*ux(x,y)+cy(k)*uy(x,y))^2-1.5*(ux(x,y)^2+uy(x,y)^2));% c^2 = 1/3 is not considered review
                %collision
                f(k,x,y)=f(k,x,y)*(1-omega)+feq(k)*omega;
                %streaming
                newx = 1+mod(x-1+cx(k)+nx,nx);
                newy = 1+mod(y-1+cy(k)+ny,ny);
                fs(k,newx,newy) = f(k,x,y);
            end
        end
    end
     f =fs;
    %boundary conditions
            %bounceback i.e no slip for walls
    for y=1:ny
         for x=1:nx
             if y==ny % top wall
              rhow = f(1,x,y)+f(2,x,y)+f(4,x,y)+2*(f(7,x,y)+f(3,x,y)+f(6,x,y));
               % uy(x,y)= ((f(1,x,y)+f(2,x,y)+f(4,x,y)+2*(f(3,x,y)+f(6,x,y)+f(7,x,y)))/rho0)-1;
               % f(5,x,y)= f(3,x,y)-(2/3)*rho0*uy(x,y);
                f(5,x,y)= f(3,x,y);
               % f(8,x,y)= f(6,x,y)+0.5*(f(2,x,y)-f(4,x,y))-(1/6)*rho0*uy(x,y)-0.5*rho0*ux(x,y);
                f(8,x,y)= f(6,x,y)-0.5*rhow*u0+0.5*(f(2,x,y)-f(4,x,y));
               % f(9,x,y) =f(7,x,y)-0.5*(f(2,x,y)-f(4,x,y))-(1/6)*rho0*uy(x,y)+0.5*rho0*ux(x,y);
                f(9,x,y) =f(7,x,y)+0.5*rhow*u0-0.5*(f(2,x,y)-f(4,x,y));
             end
             
             if y==1 % bottom wall
                 f(3,x,y)=f(5,x,y);
                 f(6,x,y)=f(8,x,y);
                 f(7,x,y)=f(9,x,y);
             end
             
           % if x==1 % symmetry
            %    f(7,x,y)=f(7,nx,y);
             %   f(4,x,y)=f(4,nx,y);
              %  f(8,x,y)=f(8,nx,y);
                
            %end
            %if x==nx % Outlet symmetry
             %   f(6,x,y)=f(6,1,y);
             %   f(2,x,y)=f(2,1,y);
              %  f(9,x,y)=f(9,1,y);
            %end
         end
     end
    i = 1:1:ny;
    plot(i,ux(15,i));
    hold on;
    
end


    
                    
			
        
                 
      
   
            
            
    
    