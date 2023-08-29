classdef LDhs < geometry.earthModel
    properties
        % rigidity
        G;
        % Poisson's ratio
        nu;
    end
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function o=LDhs(G,nu)
            % LDhs is a class providing the necessary functions to
            % compute the stress interactions between sources and receivers
            % and among the receivers.
            %
            %   earthModel = LDhs(G,nu);
            %
            % where G is rigidity and nu is the Poisson's ratio of a
            % homogeneous elastic half space for 2-d line displacement.
            %
            
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