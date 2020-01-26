function model=simdata(x,loc,scale,sigma,shift,nsize)

T=size(x,1);

Y1=NaN(nsize,T);
Y2=NaN(nsize,T);

basis=NaN(size(loc,1),T);
for i=1:size(basis,1)
    basis(i,:)=normpdf(x',loc(i),scale(i));
end

for i=1:nsize
    wgt1=normrnd(30,1,4,1);
    wgt1(1)=1.75*trnd(2)+30;
    wgt1(3)=normrnd(30,0.4)+0.5;
    mean1=sum(bsxfun(@times,basis,wgt1));
    Y1(i,:)=normrnd(mean1,sigma)+shift;
end

for i=1:nsize
    wgt2=normrnd(30,1,4,1);
    wgt2(3)=30+1/gamrnd(1,1/0.35);
    mean2=sum(bsxfun(@times,basis,wgt2));
    Y2(i,:)=normrnd(mean2,sigma)+shift;
end

model.Y=vertcat(Y1,Y2);
model.n=size(model.Y,1);
model.T=size(model.Y,2);

model.x0=x;