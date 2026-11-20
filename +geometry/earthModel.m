classdef earthModel < handle
% EARTHMODEL  Empty base class for elastic earth models.
%
%   o = geometry.earthModel()
%
%   Base class of the earth models used by the patch classes. It holds no
%   properties or methods; geometry.LDhs is the subclass used in I2EC.
%
% PROPERTIES
%   none
%
% SEE ALSO
%   geometry.LDhs
%
% AUTHORS
%   Rishav Mallick, 2023

    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function o=earthModel()
            % EARTHMODEL  Create an empty earth model.
            %
            %   o = geometry.earthModel()
            %
            %   Returns an empty object. Subclasses such as geometry.LDhs
            %   add the elastic properties.
            %
            % INPUTS
            %   none
            %
            % OUTPUTS
            %   o       geometry.earthModel object
            %
            % SEE ALSO
            %   geometry.LDhs
            %
            % AUTHORS
            %   Rishav Mallick, 2023

            if (0==nargin)
                return
            end

        end
    end % methods
end % class definition
