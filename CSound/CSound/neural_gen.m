%-------------------------Neural excitation model------------------------
%------------------------------------------------------------------------
%----------------------ubm= time domain u_bm data------------------
%----------------------thr =  thr u_bm motion for excitation-------
%---------------------nrt = neural relaxation time-> graded according to CF
%---------------------dt= t(2)-t(1)--------------------------------



function neuron_res=neural_gen(ubm,thr,nrt,dt)

x_len=size(ubm,1);
t_len=size(ubm,2);

neuron_ac=false(x_len,2);
neuron_res=false(x_len,t_len);
neuron_relax=zeros(x_len,1);

for t_indx=1:t_len
    
    neuron_relax(neuron_relax<=0)=0;
    
    neuron_ac(:,1)=(neuron_relax==0);                                   % Available neural excitation spots
    
    neuron_ac(:,2)=(ubm(:,t_indx)>thr).*(neuron_ac(:,1)==1);            % Excite from bm
    act_bm= (neuron_ac(:,2)==1);
    
    neuron_res(:,t_indx)=neuron_ac(:,2);                                % Record activity
    
    neuron_relax(act_bm)=neuron_relax(act_bm)+nrt(act_bm);                % Exciting neurons
    neuron_relax=neuron_relax-dt;                                         % Relaxing in time
    
   
end