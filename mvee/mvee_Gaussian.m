Num_rep=2;
Dim_array=[10,50,100];
Num_vtx=500;
Red_array=[5000,15000,30000];
Eps_array=[0.01,0.001,0.0001];
time_ta=zeros(length(Eps_array),length(Dim_array));
time_scale_ta=zeros(length(Eps_array),length(Dim_array));
time_qh=zeros(length(Eps_array),length(Dim_array));
iter_ta=zeros(length(Eps_array),length(Dim_array));
iter_scale_ta=zeros(length(Eps_array),length(Dim_array));
for ii=1:length(Dim_array)
    Num_dim=Dim_array(ii);
    Num_pts=Red_array(ii);
    for jj=1:length(Eps_array)
        epsilon=Eps_array(jj);
        
        for kk=1:Num_rep
            ii
            jj
            kk
            vtx_A=Random_pts(Num_dim,Num_vtx,'normal');
            inhulldata=Random_cvx(vtx_A,Num_pts,'dir');
            matA=[inhulldata,vtx_A];
            [m,n]=size(matA);
            normal_val_index=(Num_pts+1):n;

            rndindx=randperm(n,n);
            [or_val or_ind]=sort(rndindx);
            vertices_index=or_ind(normal_val_index);
            matA=matA(:,rndindx);
            
            disp('ta')
            tic;
            [index,iter_num_ta]=AVTA_anti(matA',0.001);
            iter_ta(jj,ii)=iter_ta(jj,ii)+iter_num_ta;
            [avta_A , avta_c] = MinVolEllipse(matA(:,index)',epsilon);
            disp('the running time of MVE with  AVTA is')
            ta_end=toc
            time_ta(jj,ii)=time_ta(jj,ii)+ta_end;
            disp('end ta')
            disp('vtx_num')
            length(index)
            disp('spta')
            tic;
            [index iter_sp_ta]=AVTA_SP(matA',0.001);
            iter_scale_ta(jj,ii)=iter_scale_ta(jj,ii)+iter_sp_ta;
            disp('vtx_num')
            length(index)
            [spavta_A , spavta_c] = MinVolEllipse(matA(:,index)',epsilon);
            disp('the running time of MVE with spherical AVTA is')
            scale_ta_end=toc
            time_scale_ta(jj,ii)=time_scale_ta(jj,ii)+scale_ta_end;
            disp('end spta')
            length(index)
            disp('mvee')
            tic;
            [A , c] = MinVolEllipse(matA',epsilon);
            disp('the running time of MVE without AVTA is')
            qh_t=toc
            time_qh(jj,ii)=time_qh(jj,ii)+qh_t;
            disp('end mvee')
            
            
            
        end
        
    end
    
    
end





names=strings(6,1);
for i=1:length((time_ta(1,:)))
    this_size=['Dim:',num2str(Dim_array(i)), ', Redundant pts:', num2str(Red_array(i))];
    names(i)=this_size;
end
        
for i=1:length((time_ta(:,1)))

    this_title=['vertice generated using Gaussian ', 'epsilon ',num2str(Eps_array(i))];
    figure(i);
    hold on
    title(this_title);

    plot(log(time_ta(i,:)),'DisplayName','AVTA and MVEE','LineWidth',1.5)
    plot(log(time_scale_ta(i,:)),'DisplayName','AVTA+ and MVEE','LineWidth',1.5)  
    plot(log(time_qh(i,:)),'DisplayName','MVEE','LineWidth',1.5)  ;
    legend('show','Location','northwest')%,'Orientation','horizontal')
    set(gca,'xtick',[1:length((time_ta(i,:)))],'xticklabel',names)
    xlabel ("Problem Size");
    ylabel ("Running time (secs in log scale)");
    save_title= strrep(this_title,' ','_');
    saveas(gcf,['mvee_',save_title,'.png'])
    hold off;
end

        


