% close all
clear
clc

tic
options = odeset('MaxStep',1, 'outputfcn', @odewbar);

%% Initial conditions for the basal cell model
Y=[];

Y(1) = -60.1; % potential
Y(2) = 0.135; % m activation gate INa
Y(3) = 0.03123; %h1 inactivation gate INa
Y(4) = 1;%4.894; %h2 inactivation gate INa
Y(5) = 0.135; % ms activation gate INas
Y(6) = 0.03123; %h1s inactivation gate INas
Y(7) = 1;%4.894; %h2s inactivation gate INas

Y(9)  = 0.02248; % d_L activation gate ICaL Cav1_2
Y(10) = 0.43; %f_L inactivation gate ICaL Cav1_2

Y(11) = 0.02248; %d_D activation gate ICaL Cav1_3
Y(12) = 0.531; %f_D inactivation gate ICaL Cav1_3

Y(13) = 0.02217; %d_T activation gate ICaT
Y(14) = 0.06274; % f_T inactivation gate ICaT

Y(15) = 0.002359413; %P_af activation fast gate Ikr
Y(16) = 0.09102082; %P_as activation slow gate Ikr
Y(17) = 0.9977152; % P_i inactivation gate Ikr

Y(18) = 0.06609; % q inactivation gate Ito
Y(19) = 0.05733; % r activation gate Ito and I_sus

Y(20) = 0.007645; %P hyperpolarisation_activated_current_y_gate

Y(21) = 0.0015225; %d_s sustained_inward_current_d_gate
Y(22) = 0.283; %f_s sustained_inward_current_f_gate

Y(23) = 0; %xs gate slow_delayed_rectifying_potassium_current

Y(24) = 0.99979196;%0; %m_cak

Y(25) = 0.64484881; %R1
Y(26) = 0.000012240116; %O1
Y(27) = 0.0000067409369; %I1

Y(28) = 7; % Update to the SS value
Y(29) = 159.43881; %Ki
Y(30) = 0.0000326892; %Cai
Y(31) = 0.0000534741; %Ca_sub
Y(32) = 4.5887879; %CaNSR
Y(33) = 3.431262; %CaJSR

Y(34) = 0.45367981;  %fTC
Y(35) = 0.93724623;  %fMC
Y(36) = 0.05542706;  %fTMM
Y(37) = 0.6364781;   %fCMi
Y(38) = 0.29785784;  %fCMs
Y(39) = 0.84628533;  %fCQ
Y(40) = 1;            %jh
Y(41) = 1;            %jhs

%% Integration method

model = @AVNmodel_new;
nsec = 10; %s
CL = 1000*nsec;

bpm_stim = 450; %CL (ms) 

flagsim2 = ['lastMDP_stim', num2str(bpm_stim),'bpm_'];
flagsim  = ['data_stim', num2str(bpm_stim),'bpm_'];
% flagsim2 = ['lastEnd_ICaDBlock_'];
% flagsim  = ['dataEndICaDBlock_'];
flag_save = 0;
% load CI_old.mat
% Y=Yfinal;
% % loading of steady state IC
load CIss_lastMDP_300s
Y=IC';
% Y(1)=Y(1);
[t,Yc] = ode15s(@(x,y) model(x, y, bpm_stim),[0 CL], Y , options);
Vm   = Yc(:,1);
dVm  = [0; diff(Vm)./diff(t)];
Yfinal = Yc(end,:);

%%
for i= 1:size(Yc,1)
    [~, dati]    = model(t(i), Yc(i,:), bpm_stim);
    INa(i)   = dati(1);
    INas(i)  = dati(2);
    ICaL(i)  = dati(3);
    ICaD(i)  = dati(4);
    ICaT(i)  = dati(5);
    Ikr(i)   = dati(6);
    Ik1(i)   = dati(7);
    Ito(i)   = dati(8);
    Isus(i)  = dati(9);
    If(i)    = dati(10);
    IbNa(i)  = dati(11);
    IbK(i)   = dati(12);
    IbCa(i)  = dati(13);
    Ip(i)    = dati(14);
    INaCa(i) = dati(15);
    Ist(i)   = dati(16);
    IKs(i)   = dati(17);
    h1(i)    = dati (18);
    h2(i)    = dati (19);
    m(i)     = dati (20);
    d_D(i)   = dati (21);
    f_D(i)   = dati (22);
    d_T(i)   = dati (23);
    f_T(i)   = dati (24);
    Istim(i) = dati (25);
    d_s(i)   = dati (26);
    f_s(i)   = dati (27);

    %gates INas
    m_inf(i)  = dati(28);
    tau_m(i)  = dati(29);
    h1_inf(i) = dati(30);
    tau_h1(i) = dati(31);
    tau_h2(i) = dati(32);

    %Gates of ICaD
    d_D_inf (i) = dati (33);
    f_D_inf (i) = dati (34);
    tau_d_D (i) = dati (35);
    tau_f_D (i) = dati (36);

    %Gates of ICaT
    d_T_inf (i) = dati (37);
    f_T_inf (i) = dati (38);
    tau_d_T (i) = dati (39);
    tau_f_T (i) = dati (40);

    %INas
    ms (i)  = dati (41);
    hs1 (i) = dati (42);
    hs2 (i) = dati (43);

    %ISK
    ISK (i)   = dati(44);
    m_cak (i) = dati(45);

    %Ca flux
    j_Ca_dif (i) = dati(46);
    j_up (i)     = dati(47);
    j_tr (i)     = dati(48);
    j_SRCarel (i)= dati(49);

    %IKACh
    IKACh (i) = dati(50);

    %intracellular currents
    Cai(i) = dati(51);
    Nai(i) = dati(52);
    Ki(i)  = dati(53);

    %calcium currents
    Ca_sub(i) = dati(54);
    CaNSR(i)  = dati(55);
    CaJSR(i)  = dati(56);
    h(i)  = dati(62);
    hs(i)  = dati(63);

     % Ito/IKur currents
     r(i)=dati(64);
     q(i)=dati(65);

     %If current
     P(i)=dati(66);

end
result       = [INa;INas; ICaL; ICaD; ICaT; Ikr; Ik1; Ito ; Isus; If; IbNa; IbK; IbCa; Ip; INaCa; Ist; IKs; h1; h2;m; d_D; f_D; d_T;f_T; Istim; d_s; f_s; m_inf; tau_m ; h1_inf; tau_h1; tau_h2; d_D_inf; f_D_inf; tau_d_D; tau_f_D; d_T_inf; f_T_inf; tau_d_T; tau_f_T; ms; hs1; hs2; ISK; m_cak; j_Ca_dif; j_up; j_tr; j_SRCarel; IKACh; r;q; P];
mat_correnti = [INa;INas; ICaL; ICaD; ICaT; Ikr; Ik1; Ito; Isus; If; IbNa; IbK; IbCa; Ip; INaCa; Ist; IKs; Istim; ISK; IKACh];
I_tot = sum(mat_correnti);

data= datestr(datetime('today'),'yyyymmdd_');
file = [data,flagsim,num2str(nsec),'s.mat'];
file_IC = [data,flagsim2,num2str(nsec),'s.mat'];
[mdp index_mdp] = findpeaks(-Vm);
toc
% % % 
if flag_save == 1
%     Yfinal = Yc(index_mdp(end),:);
    IC = Yfinal;
    save(fullfile(file_IC), 'IC')
    save(fullfile(file),'t','Vm', 'Cai', 'Ca_sub', 'CaNSR', 'CaJSR', 'Nai', 'I_tot', 'result' ,'Yfinal')
end

%% Plot
t = t/1e3; % ms to s
figure;
nrighe = 4;
ncolonne= 1;
pos = 1:1:nrighe;
subplot(nrighe,ncolonne, pos(1))
hold on
plot(t,Vm,'b','LineWidth', 1.2)
% xlim([nsec-2 nsec])
hold on
% plot(t(index_mdp(end)),Vm(index_mdp(end)),'r*','LineWidth', 1.2)
ylabel('Vm (mV)')
subplot(nrighe,ncolonne, pos(2))
hold on
plot(t,Cai,'b','LineWidth', 1.2)
ylabel('Ca_i(mM)')
subplot(nrighe,ncolonne, pos(3))
hold on
plot(t,Ca_sub,'b','LineWidth', 1.2)
ylabel('Ca_{sub} (mM)')
subplot(nrighe,ncolonne, pos(4))
hold on
plot(t,Nai,'b','LineWidth', 1.2)
ylabel('Na_{i} (mM)')
xlabel('t (s)')
% figure;
% plot(t,Ito ,'b',t, Isus, 'r')
figure;
plot(t,I_tot )

% figure;
% plot(t,INa)
% hold on 
% plot(t,INas,'r')
% figure;
% plot(t,ICaL)
% hold on 
% plot(t,ICaD,'r')
% figure
% subplot(411)
% plot(t,If)
% subplot(412)
% plot(t,Ip)
% subplot(413)
% plot(t,INaCa)
% subplot(414)
% plot(t,Ist)
% 
% figure
% plot(t,IKs)

%% 

% tstart = find(abs((t-299.27)) <0.001);
% tstart(2) = [];
% 
% tend = find(abs((t-299.65)) <0.001);
% tend(2) = [];
% 
% t_last = t(tstart:tend);
% Vm_last = Vm(tstart:tend);

% figure
hold on
plot(t, Vm)

writematrix([t, Vm], 'INas_100block_AP.xlsx')

