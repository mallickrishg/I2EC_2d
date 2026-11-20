classdef shearZoneReceiver < geometry.shzpatch
% SHEARZONERECEIVER  Triangular mantle mesh with viscous rheology.
%
%   shz = geometry.shearZoneReceiver(filename, earthModel)
%
%   Reads a triangular mesh (file format in geometry.shzpatch) and stores
%   the triangle vertices, centroids and areas, with zeroed rheological
%   properties and long-term strain rates for the caller to set. Each
%   triangle is a shear zone whose deviatoric viscous strain rate follows
%   the power law
%       de/dt = alpha * sigma * |sigma|^(n-1),
%       |sigma| = sqrt(sigma22'^2 + sigma23^2),
%   where sigma22' is the deviatoric normal stress and sigma23 the shear
%   stress.
%
% PROPERTIES
%   e22pl, e23pl  [M x 1] long-term deviatoric strain rates (1/s)
%   n             [M x 1] power-law stress exponent (1: Newtonian)
%   alpha         [M x 1] fluidity, 1/viscosity for n = 1 (1/(MPa s))
%   dgf           [1 x 1] degrees of freedom per triangle in the state
%                 vector of the numerical integration (2: sigma22', sigma23)
%   plus the geometry inherited from geometry.shzpatch (N, xc, tri, vert,
%   A, B, C, area, earthModel)
%
% SEE ALSO
%   geometry.shzpatch, geometry.receiver, getStrainratesLongterm
%
% AUTHORS
%   Sharadha Sathiakumar, 2023
    properties
        % long-term deviatoric strain rates e22 and e23 [M x 1] each (1/s)
        e22pl;e23pl

        % power-law rheology de/dt = alpha*sigma*|sigma|^(n-1), with
        % |sigma| = sqrt(sigma22'^2 + sigma23^2): stress exponent n and
        % fluidity alpha (1/(MPa s)), [M x 1] each
        n;alpha;
        % degrees of freedom per triangle in the numerical integration [1 x 1]
        dgf;
    end

    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = shearZoneReceiver(filename,earthModel)
            % SHEARZONERECEIVER  Create a shear-zone mesh from mesh files.
            %
            %   shz = geometry.shearZoneReceiver(filename, earthModel)
            %
            %   Reads <filename>_vertices.dat and <filename>_triangulation.dat
            %   with geometry.shzpatch.loadshztri, then computes the vertices
            %   A, B, C, the areas and the centroids of the triangles. Sets
            %   alpha, n, e22pl and e23pl to zero and dgf = 2.
            %
            % INPUTS
            %   filename    mesh prefix, e.g. 'inputs/shearzone' for
            %               inputs/shearzone_vertices.dat and
            %               inputs/shearzone_triangulation.dat; a cell array
            %               is accepted, but only its last mesh is kept
            %   earthModel  geometry.LDhs object with the elastic properties
            %
            % OUTPUTS
            %   obj         geometry.shearZoneReceiver object with M triangles
            %
            % EXAMPLE
            %   earthModel = geometry.LDhs(30e3, 0.25);
            %   shz = geometry.shearZoneReceiver('inputs/shearzone', earthModel);
            %
            % SEE ALSO
            %   geometry.shzpatch, createShearZoneMesh
            %
            % AUTHORS
            %   Sharadha Sathiakumar, 2023


            if ~iscell(filename)
                filename={filename};
            end

            obj.earthModel=earthModel;

            vert=[];
            tri=[];
            for k=1:length(filename)

                  p_filename = [filename{k} '_vertices.dat'];
                  t_filename = [filename{k} '_triangulation.dat'];
                  assert(2==exist(p_filename,'file'),['patch:file not found ' filename{k}]);
                  [vert, tri]=obj.loadshztri(p_filename,t_filename);
            end

            % set degrees of freedom by default
            obj.dgf = 2;

            % triangular mesh properties
            obj.tri = tri;
            obj.vert = vert;

            % define triangles as A,B,C - 3 vertices (x2,x3)
            obj.A = [vert(tri(:,1),1),vert(tri(:,1),2)];
            obj.B = [vert(tri(:,2),1),vert(tri(:,2),2)];
            obj.C = [vert(tri(:,3),1),vert(tri(:,3),2)];

            % calculate area and triangle centers
            obj.area = 0.5*abs(obj.A(:,1).*(obj.B(:,2) - obj.C(:,2)) + ...
                       obj.B(:,1).*(obj.C(:,2)-obj.A(:,2)) + obj.C(:,1).*(obj.A(:,2)-obj.B(:,2)));

            x2c = mean([obj.A(:,1), obj.B(:,1), obj.C(:,1)],2);
            x3c = mean([obj.A(:,2), obj.B(:,2), obj.C(:,2)],2);
            obj.xc = [x2c,x3c];

            obj.N = length(tri(:,1));

            % rheological properties
            obj.alpha = zeros(obj.N,1);
            obj.n     = zeros(obj.N,1);

            % long-term deviatoric strain rates
            obj.e22pl = zeros(obj.N,1);
            obj.e23pl = zeros(obj.N,1);

            end % constructor




    end % methods

end
