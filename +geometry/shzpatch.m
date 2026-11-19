classdef shzpatch < handle
% SHZPATCH  Base class for triangular shear-zone meshes.
%
%   obj = geometry.shzpatch()
%
%   Holds the geometry of a triangulated volume (the mantle cells of I2EC)
%   and the static reader for mesh files. geometry.shearZoneReceiver
%   subclasses it and fills the properties.
%
%   A mesh is stored in two comma-separated files that share a prefix:
%   <name>_vertices.dat has one "x,z" pair per vertex (m, z negative below
%   the surface) and <name>_triangulation.dat has three 1-based vertex
%   indices per triangle.
%
% PROPERTIES
%   N           [1 x 1] number of triangles, M
%   xc          [M x 2] triangle centroids, (x, z) (m)
%   tri         [M x 3] vertex indices of each triangle
%   vert        [Nv x 2] vertex coordinates, (x, z) (m)
%   A, B, C     [M x 2] coordinates of the three vertices of each
%               triangle (m)
%   area        [M x 1] triangle areas (m^2)
%   earthModel  geometry.LDhs object with the elastic properties
%
% SEE ALSO
%   geometry.shearZoneReceiver
%
% AUTHORS
%   Sharadha Sathiakumar, 2023

    properties
        % number of triangles [1 x 1]
        N;
        % triangle centroids [M x 2] (m), (x, z) with z negative below the
        % surface
        xc;
        % vertex indices of each triangle [M x 3] and vertex coordinates
        % [Nv x 2] (m)
        tri; vert;
        % vertices A, B, C of each triangle [M x 2] each (m) and triangle
        % area [M x 1] (m^2)
        A;B;C;area;
        % earth model: geometry.LDhs object with the elastic properties
        earthModel;
    end

    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = shzpatch()
            % SHZPATCH  Create an empty triangular-mesh object.
            %
            %   obj = geometry.shzpatch()
            %
            %   Returns an empty object; geometry.shearZoneReceiver fills it.
            %
            % INPUTS
            %   none
            %
            % OUTPUTS
            %   obj     geometry.shzpatch object
            %
            % SEE ALSO
            %   geometry.shearZoneReceiver
            %
            % AUTHORS
            %   Sharadha Sathiakumar, 2023

            if (0==nargin)
                return
            end

        end % constructor
    end % methods
    methods(Static)
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %        Load traingle indices and vertices        %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [vert,tri] = loadshztri(p_filename,t_filename)
            % LOADSHZTRI  Read the vertex and triangle files of a mesh.
            %
            %   [vert, tri] = geometry.shzpatch.loadshztri(p_filename, t_filename)
            %
            %   Reads both comma-separated files with readtable and returns
            %   their contents as numeric arrays.
            %
            % INPUTS
            %   p_filename  path to <name>_vertices.dat: one "x,z" pair per
            %               vertex (m, z negative below the surface)
            %   t_filename  path to <name>_triangulation.dat: three 1-based
            %               vertex indices per triangle
            %
            % OUTPUTS
            %   vert        [Nv x 2] vertex coordinates, (x, z) (m)
            %   tri         [M x 3] vertex indices of each triangle
            %
            % SEE ALSO
            %   geometry.shearZoneReceiver
            %
            % AUTHORS
            %   Sharadha Sathiakumar, 2023

            p = readtable(p_filename);
            t = readtable(t_filename);

            tri = t{:,:};
            vert = p{:,:};

        end
    end % methods (Static)

end
