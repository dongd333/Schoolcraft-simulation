N=1000;
D=0.0001;
L=1;
C0=1;
k=1;
mfact=0.5
A=ones(N,2);B=ones(N/mfact,2); NA=N;
A(:,1)=L*rand(N,1); B(:,1)=L*rand(N/mfact,1);
mppA=C0*L/N; mppB=L*C0*mfact/N;        % particle masses

ttime=100; dt=min([0.1/k  sqrt(D/L/10)]); time=0;

timestore=linspace(.1/k, 100/k, 100); cstore=zeros(length(timestore),3);
tstep=1
while time<ttime/k/C0
var=4*D*dt; 

dist=4*sqrt(var);
 
    BS=max(50,N/100);
    [idx r]=rangesearch(B(:,1),A(:,1),dist,'BucketSize',BS);

    for jj=1:NA     %  A particle loop 

       ma=A(jj,2);
       
       blah=idx{jj};     % This is the indices list of nearby particles
       blahr=r{jj};      % Associated radii

       if(length(blah)<1) 
            continue 
       end 

       % Per Dani, separately calculate and execute the probability for A
       % and B particles
       
         for ii=1:length(blah)   %  B particle loop
            pidx=blah(ii);       % index of B particle
            existB=B(pidx,2);    % the binary exist flag for the particle - 
                                 % need this because B particle might have just been killed
            
            v_s=existB*(1/sqrt(2*pi*var))*exp(blahr(ii)^2/(-2*var));
            probB=k*dt*B(pidx,2)*v_s;
            probA=k*dt*A(jj,2)*v_s;
            if probB>rand
               B(pidx,2)=0;      % kill particle
            end
            if probA>rand
              A(jj,2)=0;
               continue          % kill A particle and stop looking
            end
         end   %  B particle loop end
         
    end        %  A particle loop end
  apple = find(A(:,2)==0);
  A(apple,:)=[];
  apple = find(B(:,2)==0);
  B(apple,:)=[];
  
  NA=length(A(:,1));
  NB=length(B(:,1));

  time=time+dt
  A(:,1)=A(:,1)+sqrt(2*D*dt)*randn(NA,1);
  B(:,1)=B(:,1)+sqrt(2*D*dt)*randn(NB,1);
  % periodic boundaries
  apple=find(A(:,1)<0);
  A(apple,1)=L+A(apple,1);
  apple=find(A(:,1)>L);
  A(apple,1)=-L+A(apple,1);
  apple=find(B(:,1)<0);
  B(apple,1)=L+B(apple,1);
  apple=find(B(:,1)>L);
  B(apple,1)=-L+B(apple,1);
  

  if time>timestore(tstep)
     cstore(tstep,1)=time; cstore(tstep,2)=mppA*sum(A(:,2)); cstore(tstep,3)=mppB*sum(B(:,2));
     tstep=tstep+1;
  end
end
loglog(cstore(:,1),cstore(:,2),'+'); hold on;
loglog(cstore(:,1),cstore(:,3),'o');


