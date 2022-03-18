%% Inputs Operacionais  - Coloque os seus dados de:
clear all
propfile = 'propfile';        %Txt da helice
motorfile = 'motorfile';     %Txt do motor
runfile = 'runfile';


%% Operacional
system_command_string2 = ['qprop ' propfile '.txt ' , motorfile '.txt ', runfile '.txt ', '> ' 'out.txt'];
status2 = system(system_command_string2);
type ( 'out.txt')                                 %Output - Verificar se o arquivo esta correto 

                                                              %Processar os dados para o Workspace
delimiterIn = ' ';
headerlinesIn = 17;                                           %Pular o header
qprop = importdata('out.txt',delimiterIn,headerlinesIn);       %Importar a data 

plot(qprop.data(:,1), qprop.data(:,4))




%T = table(rR_vec, r_vector, dT_blade, dQ_blade, dT, dQ, Cl_vec, Cd_vec, Re_vec, Wa_vec);              %escrever no excel
%writetable(T,'BladeLoads.xlsx','Sheet', propfile2,'WriteRowNames',true);
