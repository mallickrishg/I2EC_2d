classdef shzpatch < handle
    
    properties
        N;
        % center point 
        xc;
        % triangle indices and vertices 
        tri; vert; 
        % triangles and area
        A;B;C;area;
        % earth model
        earthModel;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = shzpatch()            
            
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
            % loadshztri takes in the vertices (p_filename) and triangles (t_filename)  
            % to load the vert and tri x2,x3 coordinates
            % shz=loadshztri(p_filename, t_filename)
            
            p = readtable(p_filename);
            t = readtable(t_filename);
        
            tri = t{:,:};
            vert = p{:,:};
                       
        end
    end % methods (Static)
    
end 