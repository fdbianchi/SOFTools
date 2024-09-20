function sys = unstandardizeSys(sys,sysInfo)

% UNSTANDARDIZESYS serves to unstandardize the systems in design functions

% fbianchi - 2021-07-31
    
switch sysInfo.typ
    
    case 'lti'
        % nothing to do
        
    case 'mat'
        % to mat format
        sys = lti2mat(sys);
        
    case 'psys'
        % to psys
        auxsys = [];
        nv = size(sys,3);
        for ii = 1:nv
            auxsys = [auxsys, lti2mat(sys(:,:,ii))];
        end
        sys = psys(auxsys);
        
    case 'lpv'
        % to lpvss
        Pset = sysInfo.sys.parset;
        sys = ppss(sys,Pset);
        
    otherwise
        error('UNSTANDARDIZESYS:InputError',...
            'SYS is not a valid system description');
end
