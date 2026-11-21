function tf = isLinearFaultLaw(rcv)
% ISLINEARFAULTLAW  True if the numerical integration uses a linear viscous fault.
%
%   tf = isLinearFaultLaw(rcv)
%
%   Reads rcv.faultLaw: 'linear' (linear viscous fault, state v) gives
%   true; 'ratestrengthening' (the default, state log(v)) gives false. Any
%   other value is an error. Objects without the property (older saved
%   receivers) count as rate-strengthening.
%
% INPUTS
%   rcv     geometry.receiver object
%
% OUTPUTS
%   tf      logical scalar
%
% SEE ALSO
%   odeViscoelastic, extractRatesFromOdeState,
%   runImposedEarthquakeCycles

if ~isprop(rcv, 'faultLaw')
    tf = false;
    return
end
switch lower(char(rcv.faultLaw))
    case 'linear'
        tf = true;
    case 'ratestrengthening'
        tf = false;
    otherwise
        error('isLinearFaultLaw:unknown', ...
            'rcv.faultLaw must be ''ratestrengthening'' or ''linear'', not ''%s''.', char(rcv.faultLaw));
end
end
