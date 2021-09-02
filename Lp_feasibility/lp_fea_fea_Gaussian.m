Num_rep=10;
Dim_array=[50,100,200];
Vtx_array=[500,1000,2000];
Eps_array=[0.000001,0.0000005,0.0000001];
time_ta=zeros(length(Eps_array),length(Dim_array));
time_scale_ta=zeros(length(Eps_array),length(Dim_array));
time_lp=zeros(length(Eps_array),length(Dim_array));
time_mebalgo=zeros(length(Eps_array),length(Dim_array));
time_fw=zeros(length(Eps_array),length(Dim_array));
time_eg=zeros(length(Eps_array),length(Dim_array));
iterrr=0;
for ii=1:length(Eps_array)
    Num_dim=Dim_array(ii);
    Num_vtx=Vtx_array(ii);
    
    for jj=1:length(Eps_array)
        epsilon=Eps_array(jj);
        for kk=1:Num_rep
            ii
            jj
            kk
            iterrr=iterrr+1;
            A=Random_pts(Num_dim,Num_vtx,'normal');
            X=rand(Num_vtx,1);
            b=A*X;
            [M_con,N_var]=size(A);
            
            M=1200;

            tmp_mat=[A,zeros(M_con,1);ones(1,N_var),1];
            tmp_b=[-b;-M];

            data_mat=[tmp_mat,tmp_b;zeros(1,N_var+1),1];

            p=[zeros(Num_dim+1,1);1/(1+M)];
            disp('ta')
            tic;
            [inorout,p_prime_1,alpha_coe,dist,ta_iter]=ta_anti(data_mat,p,epsilon);
            ta_end=toc;
            time_ta(ii,jj)=time_ta(ii,jj)+ta_end;
            ta_iter
            disp('end ta')
            disp('spta')
            tic;
            [inorout,p_prime_2,alpha_coe,dist,iter_num]=Spherical_TA11(data_mat,epsilon,[0,0],p);
            scale_ta_end=toc;
            time_scale_ta(ii,jj)=time_scale_ta(ii,jj)+scale_ta_end;
            iter_num
            disp('end spta')
            
            
            
            
            disp('lp')
            tic;
            [n_1,n_2]=size(A);
            x_y_size=n_2;
            c_coe=[zeros(n_2,1)];

            Aeq=[A];

            beq=b;
            lb=zeros(x_y_size,1);
            
            [bb,fval] = linprog([],[A;-A],[b;-b],[],[],lb,[]);
            lp_t=toc;
            time_lp(ii,jj)=time_lp(ii,jj)+lp_t;
            disp('end lp')
            
            
            
        end
        
    end
    
    
end
        

names=strings(3,1);
for i=1:length((time_ta(:,1)))
    this_size=[num2str(Dim_array(i)), 'x', num2str(Vtx_array(i))];
    names(i)=this_size;
end
        
for i=1:length((time_ta(:,1)))
    this_title=['Gaussian feasible case ','epsilon ',num2str(Eps_array(i))];
    figure(i);
    hold on
    title(this_title);
    plot(log(time_ta(:,i)),'DisplayName','Triangle Algorithm','LineWidth',1.5)
    plot(log(time_scale_ta(:,i)),'DisplayName','Spherical TA','LineWidth',1.5)  
    plot(log(time_lp(:,i)),'DisplayName','Lp Solver','LineWidth',1.5)  ;
    legend('show','Location','northwest')%,'Orientation','horizontal')
    set(gca,'xtick',[1:length((time_ta(:,1)))],'xticklabel',names)
    xlabel ("Problem Size");
    ylabel ("Running time (secs in log scale)");
    save_title= strrep(this_title,' ','_');
    saveas(gcf,['lp_feasibility_',save_title,'.png'])
    hold off;
end
