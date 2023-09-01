classdef shearZoneReceiver < geometry.shzpatch
    properties
        % long-term deviatoric shear zone strain
        e22pl;e23pl
        
        % rheological properties of the form dε/dt = ασ|σ|^(n-1) where 
        % |σ| = sqrt(σxx'^2 + σxz^2), and σxx' is the deviatoric component
        n;alpha;
        % degrees of freedom (number of parameters solved in numerical integration)
        dgf;
    end
  
    methods
        %% % % % % % % % % % % % % % % % % %
        %                                  %
        %      c o n s t r u c t o r       %
        %                                  %
        % % % % % % % % % % % % % % % % % %%
        function obj = shearZoneReceiver(filename,earthModel)
            % shearZoneReceiver is a class representing the geometry and physical
            % properties of receiver shear zone, including position
            % (vertices and centers of the triangular mesh patches), area, and rheological properties.
            %
            % shz = geometry.shearZoneReceiver('filename', earthModel)
            

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
