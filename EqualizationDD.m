% Equ
% 选择Eq 类型
Eqmode = eqmodei;
xn=real(xn);
switch Eqmode
    case 0
        sigRx_E = xn;
    case 1
        EQ=struct();
        EQ.u=0.01; %  0.0001
        EQ.k1=31;
        EQ.k2=15;
        EQ.ref=8;
        EQ.sps=2;
        EQ.lamda=0.9999;
        EQ.delta=0.01;
        % FFE_LMS
        [sigRx_E,en,w] = FFE_LMS(EQ, xn, label);
    case 11
        EQ=struct();
        EQ.u=0.01;
        EQ.k1=31;
        EQ.k2=15;
        EQ.ref=8;
        EQ.sps=2;
        EQ.lamda=0.9999;
        EQ.delta=0.01;
        % FFE_RLS
        [sigRx_E,en,w] = FFE_RLS(EQ, xn, label);

    case 2
        EQ=struct();
        EQ.u=0.01;
        EQ.k1=31;
        EQ.k2=15;
        EQ.ref=8;
        EQ.sps=2;
        EQ.lamda=0.9999;
        EQ.delta=0.01;
        % DFE_LMS
        [sigRx_E,en,w] = DFE_LMS(EQ, xn, label);

    case 21

        EQ=struct();
        EQ.u=0.01;
        EQ.k1=31;
        EQ.k2=15;
        EQ.ref=8;
        EQ.sps=2;
        EQ.lamda=0.9999;
        EQ.delta=0.01;
        % DFE_RLS
        [sigRx_E,en,w] = DFE_RLS(EQ, xn, label);
    case 22
        EQ=struct();
        EQ.u=0.01;
        EQ.k1=31;
        EQ.k2=15;
        EQ.ref=8;
        EQ.sps=2;
        EQ.lamda=0.9999;
        EQ.delta=0.01;
        % APA
        [sigRx_E,en,w]=APA(EQ, xn, label);
    case 3
        % % Volterra FFE
        sps=2;
        ref=8;
        taps_list = [31 11 0 0 0 0];
        step_len = 0.001;
        lamda=0.9999;
        % volterra_lms
        [sigRx_E,en,w]=volterra_dfe_lms(xn,label,sps,ref,taps_list,step_len);
    case 31
        % % Volterra  FFE
        sps=2;
        ref=8;
        taps_list = [31 11 0 0 0 0];
        step_len = 0.001;
        lamda=0.9999;
        % volterra_rls
        [sigRx_E,en,w]=volterra_dfe_rls(xn,label,sps,ref,taps_list,lamda);
    case 4
        % % Volterra DFE
        sps=2;
        ref=8;
        taps_list = [31 11 0 15 0 0];
        step_len = 0.001;
        lamda=0.9999;
        % volterra_lms
        [sigRx_E,en,w]=volterra_dfe_lms(xn,label,sps,ref,taps_list,step_len);
    case 41
        % % Volterra DFE
        sps=2;
        ref=8;
        taps_list = [31 11 0 15 0 0];
        step_len = 0.001;
        lamda=0.9999;
        % volterra_rls
        [sigRx_E,en,w]=volterra_dfe_rls(xn,label,sps,ref,taps_list,lamda);
    case 5
        % % ABF
        sps=2;
        ref=8;
        taps_list_ABS = [31 0 0 15 5 3];
        step_len = 0.001;
        lamda=0.9999;
        % abf_rls
        [sigRx_E,en,w]=abf_lms(xn,label,sps,ref,taps_list_ABS,step_len);
    case 6
        % % PWL
        sps=2;
        ref=8;
        % pwl
        [sigRx_E,en,w] = pwl(xn,label,sps,ref,15,0.005,[-0.7661,0,0.7661]);

end