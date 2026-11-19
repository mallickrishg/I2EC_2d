classdef LDhs < geometry.earthModel
% LDHS  Homogeneous elastic half-space for 2-D line dislocations.
%
%   earthModel = geometry.LDhs(G, nu)
%
%   Holds the elastic properties of a homogeneous, isotropic elastic
%   half-space in plane strain. Every fault (geometry.receiver) and shear
%   zone (geometry.shearZoneReceiver) object stores one, and the kernel
%   functions read G and nu from it.
%
% PROPERTIES
%   G       [1 x 1] rigidity (shear modulus) (MPa)
%   nu      [1 x 1] Poisson's ratio
%
% SEE ALSO
%   geometry.earthModel, geometry.receiver, geometry.shearZoneReceiver
%
% AUTHORS
%   Rishav Mallick, 2023
    properties
        % rigidity (shear modulus) [1 x 1] (MPa)
        G;
        % Poisson's ratio [1 x 1]
        nu;
    end
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function o=LDhs(G,nu)
            % LDHS  Create a homogeneous elastic half-space model.
            %
            %   earthModel = geometry.LDhs(G, nu)
            %
            %   Checks that G >= 0 and -1 <= nu <= 0.5 and stores them.
            %   With no arguments, returns an empty object.
            %
            % INPUTS
            %   G       [1 x 1] rigidity (shear modulus) (MPa)
            %   nu      [1 x 1] Poisson's ratio
            %
            % OUTPUTS
            %   o       geometry.LDhs object
            %
            % EXAMPLE
            %   earthModel = geometry.LDhs(30e3, 0.25);
            %
            % SEE ALSO
            %   geometry.receiver, geometry.shearZoneReceiver
            %
            % AUTHORS
            %   Rishav Mallick, 2023

            if (0==nargin)
                return
            end

            assert(0<=G,'LDhs::rigidity must be positive.')
            assert(nu<=0.5,'LDhs::Poisson''s ratio should be lower than 0.5.')
            assert(-1<=nu,'LDhs::Poisson''s ratio should be greater than -1.')

            o.G=G;
            o.nu=nu;
        end

    end % methods
end % class definition
