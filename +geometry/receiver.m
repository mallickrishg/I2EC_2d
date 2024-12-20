classdef receiver < geometry.fpatch
    properties
        % friction properties
        Vo;a;b;l;
        % confining pressure or fault-normal traction (MPa)
        sigma;
        % combined rate-dependent friction parameter (a-b)sigma (MPa)
        Asigma;
        % friction coefficient at reference velocity
        mu0;
        % plate velocity
        Vpl;  
        Vx;
        Vz;
        % pinned patch positions (index)
        pinnedPosition;        
        % degrees of freedom (number of parameters solved in numerical integration)
        dgf;
    end
    
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = receiver(filename,earthModel)
            % RECEIVER is a class representing the geometry and physical
            % properties of receiver fault patches, including position, 
            % orientation, dimension and friction properties.
            %
            %   src = receiver('filename')
            % default friction properties are velocity strengthening:           
            
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
