function pad=Padding(Pi)
  L=size(Pi);
  L=L(1);
  pad=zeros(L+2,L+2);
  pad(2:L+1,2:L+1)=Pi;
  pad(1,2:L+1)=Pi(L,:);
  pad(L+2,2:L+1)=Pi(1,:);
  pad(2:L+1,1)=Pi(:,L);
  pad(2:L+1,L+2)=Pi(:,1);
  pad(1,L+2)=Pi(L,1);
  pad(1,1)=Pi(L,L);
  pad(L+2,1)=Pi(1,L);
  pad(L+2,L+2)=Pi(1,1);
end
