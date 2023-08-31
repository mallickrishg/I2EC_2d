function shz = load_viscous_wedges(filename,earthModel)

p_filename = [filename '_vertices.dat'];
t_filename = [filename '_triangulation.dat'];

p = readtable(p_filename);
t = readtable(t_filename);
vert = p{:,:};
tria = t{:,:};

% create a data structure for shear zones
shz = [];
shz.tri = tria;
shz.vert = vert;

% define triangles as A,B,C - 3 vertices (x2,x3)
shz.A = [vert(tria(:,1),1),vert(tria(:,1),2)];
shz.B = [vert(tria(:,2),1),vert(tria(:,2),2)];
shz.C = [vert(tria(:,3),1),vert(tria(:,3),2)];

% calculate area and triangle centers 
shz.area = 0.5*abs(shz.A(:,1).*(shz.B(:,2) - shz.C(:,2)) + ...
    shz.B(:,1).*(shz.C(:,2)-shz.A(:,2)) + shz.C(:,1).*(shz.A(:,2)-shz.B(:,2)));

x2c = mean([shz.A(:,1), shz.B(:,1), shz.C(:,1)],2);
x3c = mean([shz.A(:,2), shz.B(:,2), shz.C(:,2)],2);
shz.xc = [x2c,x3c];

shz.N = length(tria(:,1));

% define rheological quantities
shz.alpha = zeros(shz.N,1);
shz.n = zeros(shz.N,1);

% long-term deviatoric strain rates
shz.e22pl = zeros(shz.N,1);
shz.e23pl = zeros(shz.N,1);

% add elastic parameters
shz.earthModel = earthModel;


end
