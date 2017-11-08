%************************************************************************
%%-----------------This file implements the Post-Proc--- ----------------
%------------------------------------------------------------------------
%------------------Started by Aritra Sasmal (7/28/16)--------------------
%-----------------------asasmal@umich.edu--------------------------------
%************************************************************************

close all
clear all


%---------------------------flags----------------------------------


tecplot_flag=0;                 % Write a tecplot file
ascii_flag=1;                   % Read from ascii file

freq_sweep_flag=1;
x_sweep_flag=1;
%-----------------------Read data file for post-proc-------------------
% Construct a questdlg with three options
choice=[];

while(isempty(choice)==1)
choice = questdlg('Would you like to open data from folder or load default data?', ...
	'Menu', ...
	'Select folder','Default data','Select folder');
% Handle response
switch choice
    case 'Select folder'
        choice = 1;
    case 'Default data'
        choice = 2;
end

end

if choice==1

folder_sel=[];

while(isempty(folder_sel)==1)
foldername=uigetdir;
folder_sel=1;
end
clear folder_sel

trigger=1;
run('splash_init.m')
    
   

if(exist(foldername)~=7)
    if(exist(sprintf('%s.zip',foldername))==2)
        fprintf('\nUnzipping folder.\n');
        unzip(sprintf('%s.zip',foldername));
        fprintf('\nUnzipping complete.\n');
    else
        fprintf('\nNo such folder.\n');
        return
    end
end


fpdata=fopen(sprintf('%s/data_param.dat',foldername),'r');
if(fpdata==-1)
    fprintf('\n Cannot read dataparam.dat\n')
    return
end

tmp=str2num(fgetl(fpdata));
L=tmp(1);
H_T=tmp(2);
H_B=tmp(3);
L_h=tmp(4);

tmp=str2num(fgetl(fpdata));
ndymode=tmp(1);
n_y_mode=tmp(2);
n_z_node=tmp(3);
x_node_length=tmp(4);

tmp=str2num(fgetl(fpdata));
Freq_0=tmp(1);
Del_Freq=tmp(2);
num_freq=tmp(3);
fclose(fpdata);
%%
s.ProgressRatio=0.1;
%------------------Read all data from files------------------------------
n_x_node=floor(L/x_node_length)+1;
P_fluid=zeros(n_x_node,(n_z_node+1)*n_y_mode,num_freq);
num_y_nodes=50;



if(ascii_flag==1)                           % Read data from ascii file
    
    
    Freq=zeros(num_freq,1);
    
    for freq_indx=1:num_freq
    
        Freq(freq_indx)=Freq_0+(freq_indx-1)*Del_Freq;
      
        
        fpdata=fopen(sprintf('%s/Freq_indx_%d.dat',foldername,freq_indx),'r');
        
        if(fpdata==-1)
        
            fprintf(strcat('\n Cannot read Freq_indx_',num2str(freq_indx),'.dat'))
            return
        end

        X = textscan(fpdata, '%f', 'HeaderLines',0);
        X = X{1};
         
        %-----------------x axis map-------------
        
        X_dist=linspace(x_node_length,(n_x_node-1)*x_node_length,n_x_node-1);
       
        %-----------------Stapes------------------
        
        u_stapes(freq_indx)=X(2);
        
        %-----------------RW----------------------
        
        u_rw(freq_indx)=X(1);
        
        %-----------------Organ of Corti----------
        
        for cross_sec_indx=1:n_x_node-1
            
            
            x_gauss=cross_sec_indx*1.85/n_x_node;
            Mesh_x(cross_sec_indx,1:num_y_nodes)=x_gauss;
           
            b(cross_sec_indx) = 0.008+(0.018-0.008)*x_gauss/1.85;     %local width of the BM
            Mesh_y(cross_sec_indx,:)=linspace(-b(cross_sec_indx)/2,b(cross_sec_indx)/2,num_y_nodes);
            
            
            % first_ooc_node refers to the first OoC node in each cross-section. Right
            % now it starts with u_bm in each cross-section. 2 is added to
            % compensate for the stapes and RW nodes.
           
            first_ooc_node=2+(cross_sec_indx-1)*(9+(n_z_node+1)*n_y_mode)+0.5*(n_z_node+1)*n_y_mode+1;
           
            
            P(cross_sec_indx,n_y_mode*0.5*(n_z_node+1):1,freq_indx)=X(first_ooc_node-1:first_ooc_node-(0.5*(n_z_node+1)*n_y_mode));
            P(cross_sec_indx,n_y_mode*0.5*(n_z_node+1)+1:n_y_mode*(n_z_node+1),freq_indx)=X(first_ooc_node+9 ...
                :first_ooc_node+8+(0.5*(n_z_node+1)*n_y_mode));
 
            u_bm(cross_sec_indx,freq_indx)=X(first_ooc_node);
            
        
        end
        
    s.ProgressRatio=0.1+0.9*freq_indx/num_freq;
    fclose(fpdata);    
    end
end
        

Freq=Freq./1000;

%%------------------------Close splash------------------------------

trigger=2;
run('splash_init.m')


else
 
trigger=1;
run('splash_init.m')    
s.ProgressRatio=0.1;   
load defaultdata.mat     
    
    
end

%% ------------------------Input waveform-----------------------------
trigger=1;
while(1)

run('maingui.m')
handles=guihandles(maingui);

uiwait(handles.figure1);

if(trigger==0)
    close(gcf)
    return
    
end

[s,f_sig,t_sig,p]=spectrogram(sig,128,120,128,(max(sigtime)-min(sigtime))\(length(sigtime)));



[f,ufft]=FFT(sigtime,sig,'nodisplay');
f=f(3:end);
f=f-5;
ufft=ufft(3:end);

%%--------------------------Inverse FFT-------------------------------



Freq_stimulation=[0];            % in Hz, 0->Impulse response

u_bm_td=zeros(n_x_node-1,num_y_nodes,length(Freq));
P_bm_td=zeros(n_x_node-1,num_y_nodes,length(Freq));
u_bm_tmp=zeros(length(Freq),1);
P_td1=zeros(length(Freq),1);

for cross_sec_indx=1:n_x_node-1
        
            
            first_ooc_node=2+(cross_sec_indx-1)*(9+(n_z_node+1)*n_y_mode)+0.5*(n_z_node+1)*n_y_mode+1;
            
            
            [time,u_bm_tmp(1:end,1)]=FR2TR(Freq*1000,(u_bm(cross_sec_indx,:).*ufft).',Freq_stimulation,'nodisplay');
            
   
               z_indx=0.5*(n_z_node+1)*n_y_mode-n_y_mode+1;               
               P_tmp=squeeze(P(cross_sec_indx,z_indx,:))-squeeze(P(cross_sec_indx,z_indx+n_y_mode,:));                
               [~,P_td1]= FR2TR(Freq*1000,(P_tmp.*ufft.'),Freq_stimulation,'nodisplay');
               
               if(n_y_mode>1)
               z_indx=0.5*(n_z_node+1)*n_y_mode-n_y_mode+2;               
               P_tmp=squeeze(P(cross_sec_indx,z_indx,:))-squeeze(P(cross_sec_indx,z_indx+n_y_mode,:));     
               [~,P_td2]= FR2TR(Freq*1000,(P_tmp.*ufft.'),Freq_stimulation,'nodisplay');
               end
               
               if(n_y_mode>2)
               z_indx=0.5*(n_z_node+1)*n_y_mode-n_y_mode+3;               
               P_tmp=squeeze(P(cross_sec_indx,z_indx,:))-squeeze(P(cross_sec_indx,z_indx+n_y_mode,:));     
               [~,P_td3]= FR2TR(Freq*1000,(P_tmp.*ufft.'),Freq_stimulation,'nodisplay');
               end

               
               %-----------------Recreating whole cross-section-----------
               
               for y_indx=1:num_y_nodes
                   
                   y_gauss=-b(cross_sec_indx)/2+(y_indx-1)*b(cross_sec_indx)/num_y_nodes;
                   
                   u_bm_td(cross_sec_indx,y_indx,:)=u_bm_tmp(1:end,1).*sin(pi*(b(cross_sec_indx)/2+y_gauss)/b(cross_sec_indx));
                   
                   P_bm_td(cross_sec_indx,y_indx,:)=P_td1;
                   
                   if(n_y_mode>1)
                       P_bm_td(cross_sec_indx,y_indx,:)=P_bm_td(cross_sec_indx,y_indx,:)+P_td2*cos(2*pi*(0.05+y_gauss)/0.1);
                   end
                   if(n_y_mode>2)
                       P_bm_td(cross_sec_indx,y_indx,:)=P_bm_td(cross_sec_indx,y_indx,:)+P_td3*cos(4*pi*(0.05+y_gauss)/0.1);
                   end
                   
               end
            
end


%%


t_indx=1;
%--------------------Scaling and mesh------------------------
motion_amp=0.02/max(max(max(u_bm_td)));
umin=(min(min(min(u_bm_td))));
umax=(max(max(max(u_bm_td))));

Pmin=min(min(min(min(P_bm_td))));
Pmax=min(max(max(max(P_bm_td))));


%-------------------Neural-------------------------

[thr,I]=max(abs(u_bm),[],2);

neural_relax_time=1./2./Freq(I)./1000;

R=squeeze(u_bm_td(:,floor(num_y_nodes/2),:));
YU=neural_gen(R,thr./5,neural_relax_time,time(2)-time(1));
    
%---------------------------------------------------
h=figure('visible','on','DoubleBuffer','On','Color',[1 1 1]);%,'units','normalized','outerposition',[0 0 1 1]);

caxis manual;          % allow subsequent plots to use the same color limits
colormap hsv
an = linspace(0,2*pi);
y_ow = 1/5*cos(an);
z_ow = 1/5*sin(an);
x_ow = 0*an;



F= struct('cdata',[],'colormap',[]);
set(h, 'nextplot', 'replacechildren');  


for t_indx=2:(length(Freq)*time_frac)
    
    P_t=squeeze(P_bm_td(:,:,t_indx))./(max(abs(Pmax),abs(Pmin)));
    u_t=squeeze(u_bm_td(:,:,t_indx))./thr;
    
    subplot(2,2,[1,3])
    plot3(x_ow,y_ow,z_ow+1/3,'r')
    hold on
    plot3(x_ow,y_ow,z_ow-1/3,'r')
    line([Mesh_x(end),Mesh_x(end)],[-b(end)/b(1)/2, b(end)/b(1)/2],[0,0],'Color','blue')
    
    plot3(x_ow,y_ow,z_ow-4+1/3,'r')
    plot3(x_ow,y_ow,z_ow-4-1/3,'r')
    line([Mesh_x(end),Mesh_x(end)],[-b(end)/b(1)/2, b(end)/b(1)/2],[-4,-4],'Color','blue')
    hold on
    surf(Mesh_x,Mesh_y./b(1),u_t,'edgecolor','none')
    shading interp
    caxis([-1 1]);         % set the color axis scaling to your min and max color limits
    zlim([-6 1.5])
    set(gca,'ztick',[])
    set(gca,'ytick',[])

 
    [~,t]=contourf(Mesh_x,Mesh_y./b(1),P_t,50,'EdgeColor','none');
    t.ContourZLevel = -4;
    view([20 10])
    hold off
    zlabel('Pressure on BM                 BM disp.','FontSize',15)
    title('Cochlear Profile','FontSize',10)
    
    t_lapse=10*(t_indx>10)+(t_indx-1)*(t_indx<=10);
    
    t=subplot(2,2,2)    ;
    map=[1 1 1; 1 0 0];
    contourf(time(t_indx-t_lapse:t_indx),Mesh_x(:,1),YU(:,t_indx-t_lapse:t_indx),2,'edgecolor','none')
    colormap(t,map);
    caxis([0 1])
    xlabel('Time (s)','FontSize',10)
    ylabel('Position','FontSize',10)
    title('Neural activity','FontSize',10)
    axis on
    set(gca,'fontsize',10) 
    
    
    t=subplot(2,2,4);
    contourf(t_sig(1:find(t_sig>time(floor(length(Freq)/20)),1)),f_sig,p(:,1:find(t_sig>time(floor(length(Freq)/20)),1)))
    line_pos=find(t_sig>time(t_indx),1);
    hold on
    line([t_sig(line_pos),t_sig(line_pos)],[0, f_sig(end)],'Color','green')
    hold off
    ylabel('Freq (Hz)','FontSize',10)
    xlabel('Time (s)','FontSize',10)
    title('Input Spectrogram','FontSize',10)
    axis on
    set(gca,'fontsize',10) 
   
    F(t_indx-1)=getframe(h);
    
end

close(h);

h=figure;
flag=0;

while(flag==0)
choice = questdlg('Would you like to view or save movie?', ...
	'Menu', ...
	'View','Save','Continue','View');
% Handle response
switch choice
    case 'View'
        movie(gcf,F);
        
    case 'Save'
        prompt = {'Enter filename'};
        dlg_title = 'Input';
        num_lines = 1;
        defaultans = {'Sim.avi'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        vidObj=VideoWriter(answer{1});
        open(vidObj);
        writeVideo(vidObj,F);
        close(vidObj);
        
    case 'Continue'
            flag=1;
end
end
close(h);
end
