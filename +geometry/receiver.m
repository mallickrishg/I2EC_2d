classdef receiver < geometry.fpatch
% RECEIVER  Fault or boundary mesh of straight patches with rheology.
%
%   rcv = geometry.receiver(filename, earthModel)
%
%   Reads one or more .seg files (format in geometry.fpatch), splits the
%   segments into patches and stores their geometry, unit vectors and
%   default rheological properties. I2EC uses this class for the
%   megathrust, for the model boundary of the boundary-element correction
%   and for the kinematic hinge fault.
%
%   For the model boundary, Vpl is read as the boundary-condition type of
%   each patch (0: displacement, Dirichlet; 1: traction-free, Neumann) and
%   Vx, Vz as the prescribed displacements. The BEM kernels solve the
%   homogeneous problem and do not use Vx, Vz yet; see
%   computeFaultTractionKernelsBem.
%
% PROPERTIES
%   Vo, a, b, l     [N x 1] rate-and-state friction parameters: reference
%                   velocity (m/s), a, b, and characteristic slip (m);
%                   defaults only, not used by the I2EC solvers
%   sigma           [N x 1] fault-normal stress (MPa); default only
%   Asigma          [N x 1] rate-dependence (a-b)*sigma (MPa) for the
%                   rate-strengthening fault law; for faultLaw = 'linear'
%                   and in the semi-analytical solvers, a viscosity per
%                   unit length eta' (MPa s/m), with tau = eta' * v
%   faultLaw        fault law of the numerical integration (char):
%                   'ratestrengthening' (default), d log(v)/dt =
%                   (dtau/dt)/Asigma, state log(v); or 'linear' (linear
%                   viscous, as in the semi-analytical solvers),
%                   dv/dt = (dtau/dt)/Asigma, state v
%   mu0             [N x 1] friction coefficient at Vo; default only
%   Vpl             [N x 1] long-term slip rate (m/s), or the
%                   boundary-condition type for the model boundary
%   Vx, Vz          [N x 1] boundary-condition values in x and z
%   pinnedPosition  [N x 1] 1 for locked patches (zero slip rate between
%                   earthquakes), 0 otherwise
%   dgf             [1 x 1] degrees of freedom per patch in the state
%                   vector of the numerical integration (1: log(v), or v
%                   for faultLaw = 'linear')
%   plus the geometry inherited from geometry.fpatch (N, x, xc, W, dip,
%   slip, dv, nv, earthModel)
%
% SEE ALSO
%   geometry.fpatch, geometry.shearZoneReceiver, geometry.LDhs
%
% AUTHORS
%   Rishav Mallick, 2023
    properties
        % rate-and-state friction parameters [N x 1]: reference velocity Vo
        % (m/s), a, b, characteristic slip l (m); defaults, unused by I2EC
        Vo;a;b;l;
        % fault-normal stress [N x 1] (MPa); default, unused by I2EC
        sigma;
        % rate dependence (a-b)*sigma [N x 1] (MPa); for faultLaw 'linear'
        % and in the semi-analytical solvers, a viscosity per unit length
        % (MPa s/m)
        Asigma;
        % fault law of the numerical integration: 'ratestrengthening'
        % (default; state log(v)) or 'linear' (linear viscous; state v)
        faultLaw = 'ratestrengthening';
        % friction coefficient at the reference velocity [N x 1]; default
        mu0;
        % long-term slip rate [N x 1] (m/s); boundary-condition type for the
        % model boundary (0: displacement, 1: traction-free)
        Vpl;
        % boundary-condition value in x [N x 1]
        Vx;
        % boundary-condition value in z [N x 1]
        Vz;
        % locked patches [N x 1]: 1 locked (zero slip rate), 0 otherwise
        pinnedPosition;
        % degrees of freedom per patch in the numerical integration [1 x 1]
        dgf;
    end

    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = receiver(filename,earthModel)
            % RECEIVER  Create a patch mesh from one or more .seg files.
            %
            %   rcv = geometry.receiver(filename, earthModel)
            %
            %   Reads the segments with geometry.fpatch.seg2flt, converts
            %   depths to z negative below the surface, computes the unit
            %   vectors and patch centres, and sets default properties:
            %   velocity-strengthening friction (a = 0.01, b = a - 0.004,
            %   Vo = 0.1 m/s, l = 1e-3 m, sigma = 100 MPa, mu0 = 0.6),
            %   Asigma = 1 MPa, no locked patches and dgf = 1.
            %
            % INPUTS
            %   filename    path to a .seg file, or a cell array of paths
            %               whose patches are concatenated
            %   earthModel  geometry.LDhs object with the elastic properties
            %
            % OUTPUTS
            %   obj         geometry.receiver object with N patches
            %
            % EXAMPLE
            %   earthModel = geometry.LDhs(30e3, 0.25);
            %   rcv = geometry.receiver('inputs/megathrust2d.seg', earthModel);
            %
            % SEE ALSO
            %   geometry.fpatch, geometry.flt2flt
            %
            % AUTHORS
            %   Rishav Mallick, 2023

            if ~iscell(filename)
                filename={filename};
            end

            obj.earthModel=earthModel;

            fm=[];
            for k=1:length(filename)
                assert(2==exist(filename{k},'file'),['patch:file not found ' filename{k}]);

                [lfm]=obj.seg2flt(filename{k});
                fm=[fm;lfm];

            end

            % set degrees of freedom by default
            obj.dgf = 1;

            % patch properties
            obj.N=size(fm,1);
            obj.slip=zeros(obj.N,1);
            obj.x=[fm(:,1),-fm(:,2)];
            obj.W=fm(:,3);
            obj.dip=fm(:,4);
            obj.Vpl=fm(:,5);
            obj.Vx=fm(:,6);
            obj.Vz=fm(:,7);

            % default friction properties (velocity strengthening)
            obj.a=obj.W*0+1e-2;
            obj.b=obj.a-4e-3;
            obj.Vo=obj.W*0+1e-1;
            obj.l=obj.W*0+1e-3;
            obj.sigma=obj.W*0+1e2;
            obj.Asigma=obj.W*0+1;
            obj.mu0=obj.W*0+0.6;
            obj.pinnedPosition=obj.W*0;

            % unit vectors in the dip direction
            obj.dv=[...
                -cosd(obj.dip), ...
                +sind(obj.dip)];

            % unit vectors in the normal direction
            obj.nv=[...
                +sind(obj.dip), ...
                +cosd(obj.dip)];

            % center of fault patch
            obj.xc=[...
                obj.x(:,1)-obj.W/2.*obj.dv(:,1),...
                obj.x(:,2)-obj.W/2.*obj.dv(:,2)];
        end % constructor




    end % methods

end
