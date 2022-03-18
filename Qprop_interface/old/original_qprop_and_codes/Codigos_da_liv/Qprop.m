%Função Para Automatizar o qprop
function run = Qprop(x,propfile, motorfile, V, rpm)



system_command_string2 = ['qprop ' propfile '.txt ' , motorfile '.txt ',num2str(V) ' ' num2str(rpm) '> ' propfile '_out.txt'];%'qprop Prop24x15.txt TmotorKV190.txt 15 2450 > Helice24_out.txt';
status2 = system(system_command_string2);

run = x;

end
