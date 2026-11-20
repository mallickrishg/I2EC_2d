classdef fpatch < handle
% FPATCH  Base class for straight line patches (faults and boundaries).
%
%   obj = geometry.fpatch()
%
%   Holds the geometry shared by all line-patch meshes and the static
%   readers for segment files. geometry.receiver subclasses it and fills
%   the properties from a .seg file.
%
%   A .seg file lists straight segments, one per line; lines starting with
%   '#' are comments. Each data line has ten columns:
%       n  Vpl  x  z  Width  Dip  W0  qW  bcx  bcz
%   n is the segment index, Vpl the long-term slip rate, (x, z) the upper
%   end of the segment with z the depth, positive downward (m), Width the
%   down-dip length (m), Dip the dip (degrees), W0 the approximate width of
%   the first patch (m), qW the geometric growth factor of patch widths
%   (1 = uniform), and bcx, bcz the boundary-condition values in x and z.
%   seg2flt splits each segment into patches with geometry.flt2flt.
%
% PROPERTIES
%   N           [1 x 1] number of patches
%   id          identification number (not set by geometry.receiver)
%   x           [N x 2] upper end of each patch, (x, z) with z negative
%               below the surface (m)
%   xc          [N x 2] centre of each patch, same convention (m)
%   W           [N x 1] patch width, i.e. down-dip length (m)
%   dip         [N x 1] patch dip (degrees)
%   slip        [N x 1] patch slip (m)
%   dv          [N x 2] unit vector along dip, pointing up-dip
%   nv          [N x 2] unit normal vector
%   earthModel  geometry.LDhs object with the elastic properties
%
% SEE ALSO
%   geometry.receiver, geometry.flt2flt
%
% AUTHORS
%   Rishav Mallick, 2023

    properties
        % number of patches [1 x 1]
        N;
        % identification number
        id;
        % upper end of each patch [N x 2] (m), (x, z) with z negative below
        % the surface
        x;
        % centre of each patch [N x 2] (m), same convention as x
        xc;
        % patch width, i.e. down-dip length [N x 1] (m)
        W;
        % patch dip [N x 1] (degrees)
        dip;
        % patch slip [N x 1] (m)
        slip;
        % unit vectors along dip (pointing up-dip) and normal to each patch
        % [N x 2] each
        dv;nv;
        % earth model: geometry.LDhs object with the elastic properties
        earthModel;
    end

    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = fpatch()
            % FPATCH  Create an empty line-patch object.
            %
            %   obj = geometry.fpatch()
            %
            %   Returns an empty object; geometry.receiver fills it.
            %
            % INPUTS
            %   none
            %
            % OUTPUTS
            %   obj     geometry.fpatch object
            %
            % SEE ALSO
            %   geometry.receiver
            %
            % AUTHORS
            %   Rishav Mallick, 2023

            if (0==nargin)
                return
            end

        end % constructor

    end % methods

    methods(Static)
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %   convert segment definition to fault patches    %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [fm] = seg2flt(filename)
            % SEG2FLT  Read a .seg file and split its segments into patches.
            %
            %   fm = geometry.fpatch.seg2flt(filename)
            %
            %   Reads the ten-column segment file (see geometry.fpatch)
            %   and splits each segment into patches with geometry.flt2flt.
            %   Every patch inherits Vpl, bcx and bcz from its segment.
            %
            % INPUTS
            %   filename  path to the .seg file
            %
            % OUTPUTS
            %   fm        [N x 7] one row per patch:
            %             x, z, width, dip, Vpl, bcx, bcz, with (x, z) the
            %             upper end of the patch and z the depth, positive
            %             downward (m)
            %
            % SEE ALSO
            %   geometry.flt2flt, geometry.receiver
            %
            % AUTHORS
            %   Rishav Mallick, 2023


            fid=fopen(filename);  % open and count the number of columns
            line=strtrim(fgetl(fid));
            while (strcmp('#',line(1)))
                line=strtrim(fgetl(fid));  % line ignores all shell comments and check for the first line that contains data
            end
            fclose(fid);

            [~,Vpl,x1,x3,width,d,wo,alphaw,Vx,Vz]=...
                textread(filename,'%u %f %f %f %f %f %f %f %f %f',...
                'commentstyle','shell');

            fm=[];
            for k=1:length(x1)
                % list of patches for current segment
                flt=geometry.flt2flt([x1(k);x3(k)],width(k),d(k),wo(k),alphaw(k));
                % list of patches for all segments
                fm=[fm;[flt,Vpl(k).*ones(size(flt,1),1),Vx(k).*ones(size(flt,1),1),Vz(k).*ones(size(flt,1),1)]];
            end

        end
    end % methods (Static)
    methods(Static)
        %% % % % % % % % % % % % % % % % % % % % % % % % % %
        %                                                  %
        %   convert segment definition to fault patches    %
        %                                                  %
        % % % % % % % % % % % % % % % % % % % % % % % % % %%
        function [fm] = loadFlt(filename)
            % LOADFLT  Read a five-column patch file without splitting it.
            %
            %   fm = geometry.fpatch.loadFlt(filename)
            %
            %   Reads a file with the columns n, x, z, width, dip (lines
            %   starting with '#' are comments) and returns one row per
            %   line. Unlike seg2flt, the lines are not split into patches.
            %
            % INPUTS
            %   filename  path to the patch file
            %
            % OUTPUTS
            %   fm        [N x 4] one row per patch: x, z, width, dip, with
            %             (x, z) as written in the file (m) and dip in
            %             degrees
            %
            % SEE ALSO
            %   geometry.fpatch.seg2flt
            %
            % AUTHORS
            %   Rishav Mallick, 2023

            fid=fopen(filename);  % open and count the number of columns
            line=strtrim(fgetl(fid));
            while (strcmp('#',line(1)))
                line=strtrim(fgetl(fid));  % Open the file, get the first line
            end
            fclose(fid);

            [~,x1,x3,width,d]=...
                textread(filename,'%u %f %f %f %f',...
                'commentstyle','shell');
            fm = [x1,x3,width,d];

        end
    end % methods (Static)
end
