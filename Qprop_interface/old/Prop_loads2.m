%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Codigo - Calculo de distribuicao de Forcas na pa da helice     %
%               Codigo Desenvolvido Por Livia Souza e Pedro Luta         %  
%                         Integracao de Projeto 2021                     %
%                                                                        %      
%                                                                        %
%                                              Ref:QPROP theory document %         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs Operacionais  - Coloque os seus dados de:
clear all
motorfile = 'motorfile';     %Txt do motor

rho = 1.225;                 %Densidade do Ar
mu = 0.17800*10^-4;          %Viscosidade cinética

%Inputs 
propfile = '25_4_Subida';        %Txt da helice
D_inches = 25;
B = 4;                         %Numero de Pas
V = 12;                        %Velocidade 
propfile2 = [propfile 'emSubida'];
rpm = 2972;                  %RPM


%% Operacional
R = D_inches*0.0254/2;%0.3;                     %Raio [m]

run = Qprop(1, propfile, motorfile, V, rpm);                  %Função rodar o Qprop
type ( [propfile '_out.txt'])                                 %Output - Verificar se o arquivo esta correto 

                                                              %Processar os dados para o Workspace
Propfile = [propfile '_out.txt'];                             %Nome do arquivo da Helice
delimiterIn = ' ';
headerlinesIn = 20;                                           %Pular o header
qprop = importdata(Propfile,delimiterIn,headerlinesIn);       %Importar a data 


%% Calculo dos carregamentos

n = length(squeeze(qprop.data(:,1)));   %numero de seções ao longo do raio



for i = 2:(n+1)
    
    
  switch (i)       
      case n + 1
          
          r = R ;
          r_before = qprop.data(i-1,1);
          r_mean = (r + r_before)/2;
          dr = r - r_before;
          
          w_before = (qprop.data(i-1,6)*mu)/( rho*( qprop.data(i-1,2)));
          w = 0;
          W = (w_before + w)/2;
          
          phi_before= asind(w_a_before/w_before);
          phi = 0;
          phi_mean = (phi_before + phi)/2;
          
          cl_before=qprop.data(i-1,4);                              %Coeficiente anterior a ponta de pá
          cd_before=qprop.data(i-1,5);
          c_before=qprop.data(i-1,2);       
          cl= 0;                                                    %Coeficientes atuais                        
          cd= 0;
          c= 0;         
          cl = (cl + cl_before)/2;                                  % media do cl local para a seção de raio r e r+1
          cd =(cd + cd_before)/2;                                   % media do cd local para as seção de raio r e r+1  
          c = (c + c_before)/2;                                     % media da corda 
               
      otherwise
          
          r = qprop.data(i,1) ;
          r_before = qprop.data(i-1,1);
          r_mean = (r + r_before)/2;
          dr = r - r_before;
          
          w_before= (qprop.data(i-1,6)*mu)/( rho*( qprop.data(i-1,2)));
          w = (qprop.data(i,6)*mu)/( rho*( qprop.data(i,2)));
          W = (w_before + w)/2;
  
          w_a_before = qprop.data(i-1,10);
          w_a = qprop.data(i,10);
          phi_before =asind(w_a_before/w_before);
          phi =(asind(w_a/w)); 
          phi_mean = (phi_before + phi)/2;
          
          cl_before=qprop.data(i-1,4);                                  %Coeficiente anterior a ponta de pá
          cd_before=qprop.data(i-1,5);
          c_before=qprop.data(i-1,2);       
          cl=qprop.data(i,4);                                           %Coeficientes atuais                        
          cd=qprop.data(i,5);
          c=qprop.data(i,2);  
          cl = (cl + cl_before)/2;                                      %media do cl local para a seção de raio r e r+1
          cd =(cd + cd_before)/2;                                       %media do cd local para as seção de raio r e r+1  
          c = (c + c_before)/2;                                         %media da corda 

  end

r_vector(i-1, 1) = r_mean; 
dL(i-1,1) = 1/2*B*rho*(W^2)*cl*c*dr;                                      %calculo sustentacao infinitesimal para secao r_mean
dD(i-1,1) = 1/2*B*rho*(W^2)*cd*c*dr;                                      %calculo arrasto infinitesimal para secao de r_mean
 
dT(i-1,1) = dL(i-1,1)*cosd(phi_mean) - dD(i-1,1)*sind(phi_mean);              %calculo da derivada da tracao para a secao r_mean
dQ(i-1,1) = (dL(i-1,1)*sind(phi_mean) + dD(i-1,1)*cosd(phi_mean))*r_mean;     %calculo da derivada do torque  para a secao r_me
dQr(i-1,1) = dL(i-1,1)*sind(phi_mean) + dD(i-1,1)*cosd(phi_mean);
  
end

dT_blade= dT./B;                                                         %calculo da derivada da tracao para uma pa
dQ_blade= dQr./B;                                                         %calculo da derivada do torque para uma pa
rR_vec = r_vector./R;

T = sum (dT);                                                            %Integrar todos os dT
Q = sum (dQ);                                                            %Integrar todos os dQ

Cl_vec = squeeze(qprop.data(:,4));
Cd_vec = squeeze(qprop.data(:,5));
Re_vec = squeeze(qprop.data(:,6));
Wa_vec = squeeze(qprop.data(:,10));

T = table(rR_vec, r_vector, dT_blade, dQ_blade, dT, dQ, Cl_vec, Cd_vec, Re_vec, Wa_vec);              %escrever no excel
writetable(T,'BladeLoads.xlsx','Sheet', propfile2,'WriteRowNames',true);
