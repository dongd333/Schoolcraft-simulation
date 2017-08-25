ext=.2; npts=101; dx=2*ext/(npts-1);
x=linspace(-ext,ext,npts);
y=x; z=x;
[X,Y]=meshgrid(x,y);
[x3, y3, z3]=meshgrid(x,y,z);

dt=[.05 .1 .5 1 5]

D=.0001;

for i=1:length(dt)
    figure(1)
    dens=exp(-x.^2/(8*D*dt(i)))./sqrt(8*pi*D*dt(i));
    semilogy(x, dens ,'-o'); hold on;
    int=sum(dens)*dx

    axis([-ext ext 10^-15 10^15]);
    blah=exp(-(X.^2+Y.^2)/(8*D*dt(i)))./(8*pi*D*dt(i));
    
    figure(2)
    semilogy(x,blah(:,floor(npts/2+1)),'-xr'); hold on;
    int2=sum(sum(blah))*dx*dx
    axis([-ext ext 10^-15 10^15]);
    blah3=exp(-(x3.^2+y3.^2+z3.^2)/(8*D*dt(i)))./(8*pi*D*dt(i)).^3/2;
    
    figure(3)
    semilogy(x,blah3(:,floor(npts/2+1),floor(npts/2+1)),'-db'); hold on;
    int2=sum(sum(blah))*dx*dx
    axis([-ext ext 10^-15 10^15]);

end
   hold off
   
   
 