function nM = update_bloch(M, Beff, T1, T2, Mz0, gamma, dt)
% UPDATE_BLOCH - 
%   

if size(Beff,1) == 3
    Beffmag = sqrt(sum(Beff.^2));
    dalpha = (gamma*dt)*Beffmag;
    Beffmag(Beffmag==0) = 1; %fix to avoid divide by zero below
    nM = rotvec(M, Beff./([1;1;1]*Beffmag), -dalpha);
    R2 = exp(-dt./T2);
    nM(1,:) = nM(1,:).*R2;
    nM(2,:) = nM(2,:).*R2;
    nM(3,:) = Mz0 - (Mz0 - nM(3,:)).*exp(-dt./T1);
else
    dalpha = (-gamma*dt)*Beff;
    cosmdalpha = cos(dalpha);
    sinmdalpha = sin(dalpha);
    
    R2 = exp(-dt./T2);
    nM(1,:) = (cosmdalpha.*M(1,:) - sinmdalpha.*M(2,:)).*R2;
    nM(2,:) = (sinmdalpha.*M(1,:) + cosmdalpha.*M(2,:)).*R2;
    nM(3,:) = Mz0 - (Mz0 - M(3,:)).*exp(-dt./T1);
end




