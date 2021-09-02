Num_rep=5;
Vtx_array=[50,100,500];
Dim_array=[500,1000,2500];
Eps_array=[0.5];
time_ta=zeros(length(Eps_array),length(Dim_array));
time_scale_ta=zeros(length(Eps_array),length(Dim_array));
time_lp=zeros(length(Eps_array),length(Dim_array));
iterrr=0;
for ii=1:length(Dim_array)
    Num_dim=Dim_array(ii);
    Num_vtx=Vtx_array(ii);
    
    for jj=1:length(Eps_array)
        epsilon=0.001;
        gap_val=Eps_array(jj);
        for kk=1:Num_rep
            ii
            jj
            kk
            
            iterrr=iterrr+1;
            A=Random_pts(Num_dim,Num_vtx,'normal');
%             A=Random_pts(Num_dim,Num_vtx,'unif');
            X=rand(Num_vtx,1);
            b=A*X+gap_val*ones(Num_dim,1);
            [M_con,N_var]=size(A);
            


            data_mat=sparse([A',zeros(N_var,1);b',1]);
            p=zeros(Num_dim+1,1);
            p=zeros(N_var+1,1);
            disp('ta')
            tic;
            [inorout,p_prime_1,alpha_coe,dist,ta_iter]=ta_anti(data_mat,p,epsilon);
            ta_end=toc;
            time_ta(jj,ii)=time_ta(jj,ii)+ta_end;
            disp('end ta')
            disp('spta')
            tic;
            [inorout,p_prime_2,alpha_coe,dist,iter_num]=Spherical_TA11(data_mat,epsilon,[0,0],p);
            scale_ta_end=toc;
            time_scale_ta(jj,ii)=time_scale_ta(jj,ii)+scale_ta_end;
            
            disp('end spta')
            
            
            
            
            
            disp('lp')
            tic;
            [n_1,n_2]=size(A);
            x_y_size=n_1+n_2;
            c_coe=[zeros(n_2,1);-1*ones(n_1,1)];

            Aeq=[A,eye(n_1)];

            beq=b;
            lb=zeros(x_y_size,1);
            
            [bb,fval] = linprog(c_coe,[],[],Aeq,beq,lb,[]);
            lp_t=toc;
            time_lp(jj,ii)=time_lp(jj,ii)+lp_t;
            disp('end lp')
            
            
            
        end
        
    end
    
    
end
        

names=strings(3,1);
for i=1:length((time_ta(1,:)))
    this_size=[num2str(Dim_array(i)), 'x', num2str(Vtx_array(i))];
    names(i)=this_size;
end
        
for i=1:length(length(Eps_array))
    i
    this_title=['Gaussian feasible case '];
    figure(i);
    hold on
    title(this_title);
    plot(log(time_ta(i,:)),'DisplayName','Triangle Algorithm','LineWidth',1.5)
    plot(log(time_scale_ta(i,:)),'DisplayName','Spherical TA','LineWidth',1.5)  
    plot(log(time_lp(i,:)),'DisplayName','Lp Solver','LineWidth',1.5)  ;
    legend('show','Location','northwest')%,'Orientation','horizontal')
    set(gca,'xtick',[1:length((time_ta(1,:)))],'xticklabel',names)
    xlabel ("Problem Size");
    ylabel ("Running time (secs in log scale)");
    save_title= strrep(this_title,' ','_');
    saveas(gcf,['strict_lp_feasibility_',save_title,'.png'])
    hold off;
end

