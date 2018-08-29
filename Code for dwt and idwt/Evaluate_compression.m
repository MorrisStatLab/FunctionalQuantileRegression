  clear X x;
    alpha=[.95,.975,.99,.995,.999,.9999];
    xticks=[200,300,400,500,1000,2000,3000,4000,5000,7500,10000,20000];
    x_tick=[.95+(0:49)*.001,(.999+(1:9)*.0001)];
    x_length=length(x_tick)-1;
    En=single(D);
    total_Energy=zeros(1,size(D,1));
    for i=1:n
        [Csort,ix]=sort(-abs(D(i,:)));
        total_Energy(i)=sum(Csort.^2);
        energy=cumsum(Csort.^2)/total_Energy(i);
        En(i,ix)=energy;
    end;
    'Done Computing Relative Energy ',

    Ncoeffs=zeros(size(D,1),x_length+1);
    for j=1:(x_length+1);
        temp=sum(En<x_tick(j));
        for i=1:n;    
            Ncoeffs(i,j)=sum(temp>(i-1));
        end;
    end;

    'Done Computing # coefficients ', toc

    En_min=zeros(size(Ncoeffs));
    for j=1:(x_length+1)
        temp=sum(En<x_tick(j));
        for i=1:n
            keep=(temp>(i-1));
            if (sum(keep)>1)
                En_min(i,j)=min(sum(D(:,keep)'.^2)./total_Energy);
            end;
            if (sum(keep)==1)
                En_min(i,j)=min(D(:,keep)'.^2./total_Energy);
            end;
            if (sum(keep)==0)
                En_min(i,j)=0;
            end;
        end;
    end;

    'Done Computing Minimum Energies', toc

    settings=[x_tick',x_tick',x_tick',x_tick'];
    for k=1:length(x_tick)
        min_energy=x_tick(k);
        rows=1:n;
        i=rows(sum(En_min'>min_energy)>0);
        j=(length(x_tick)+1)-sum(En_min(i,:)'>min_energy);
        compression=(diag(Ncoeffs(i,j)));
        istar=min(rows(compression==min(compression)));
        settings(k,3)=istar;
        settings(k,4)=x_tick(j(istar));
        settings(k,2)=compression(istar);
    end;

    color_wheel=[1 0 0;0 1 0;1 0 1;0 1 1;1 1 0];
    figure(4)
    ncoeffs=settings(:,2);
    plot(ncoeffs,settings(:,1),'b-','LineWidth',4)
    ylabel('Minimum % Energy Preserved','Fontsize',16)
    xlabel('Number of Coefficients (K*)','Fontsize',16)
    title('Minimum % Energy Preserved vs. Number of Coefficients','Fontsize',16)
    set(gca,'XLim',[min(ncoeffs)*0.9,max(ncoeffs)*1.05],'Xscale','log','Xtick',xticks,'Fontsize',16)
    %set(gca,'YLim',[min(x_tick)-.001,max(x_tick)],'Yscale','log','Ytick',yticks,'Fontsize',16)
    hold on
    set(0,'DefaultAxesLineStyleOrder','--|--|--|--|--|--|--')
    alpha_row=zeros(length(alpha),1);
    for i=1:length(alpha)
        alpha_row(i)=sum(settings(:,1)<alpha(i))+1;
        a=[ncoeffs(1:alpha_row(i));wrev(repmat(ncoeffs(alpha_row(i)),alpha_row(i),1))];
        b=[(repmat(alpha(i),alpha_row(i),1));wrev(settings(1:alpha_row(i),1))];
        plot(a,b,'Linewidth',2,'Color',color_wheel(mod(i-1,size(color_wheel,1))+1,:))
    end;
    text(min(ncoeffs),0.98,['K*=',num2str(size(D,2))],'FontSize',16)

    hold off

    toc

    results=settings(alpha_row,:);
    keep=zeros(size(results,1),size(D,2));
    for i=1:size(results,1)
        keep(i,:)=sum(En<results(i,4))>=results(i,3)==1;
    end;

    wavespecs.keep=keep(sum(alpha<=wavespecs.compress),:);
    %wavespecs.compress=compress;
    D_compressed=D(:,wavespecs.keep==1);

    [wavespecs.Kj,wavespecs.Kj_all]=get_Kj_compress(wavespecs);
    wavespecs.J_all=length(wavespecs.Kj_all);
    wavespecs.J=length(wavespecs.Kj);
    wavespecs.K=sum(wavespecs.Kj);
end