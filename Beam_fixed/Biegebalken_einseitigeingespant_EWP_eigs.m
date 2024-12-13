% -----------------------------------------
% Lineare Schwingungen
% 
% Bsp.: einseitig eingespannter Balken
% ====================================
% 
% Berechnung Eigenfrequenzen + Moden
% Visualisierung
% 
% -----------------------------------------

% Quellen
% verschiedene... u.a.
% "Structural Dynamics af Tuning Fork" --> https://de.mathworks.com/help/pde/ug/structural-dynamics-of-tuning-fork.html


clear; close all;


%% parameters of the geometry
L = 0.1;
B = 0.005;
H = 0.005;
E = 210e9; nu = 0.3; rho = 8000;


%% define FE-problem
% geometry
gm = multicuboid(L,B,H, Zoffset=-H/2);

% setup fe-model
sModel = femodel(AnalysisType='structuralModal', Geometry=gm);
sModel.MaterialProperties = materialProperties(YoungsModulus=E, PoissonsRatio=nu, MassDensity=rho);

% Boundary Conditions
    % BC_vertices_indx = [];% rechtes Ende: BC_vertices_indx =[2 3 5 8];   % rechtes Ende: BC_vertices_indx = [1 4 6 7]
    
    % sModel.VertexBC(BC_vertices_indx) = vertexBC(Constraint="fixed");
    BC_faces_indx =[5];
    sModel.FaceBC(BC_faces_indx) = faceBC(Constraint="fixed");

% create mesh
sModel = generateMesh(sModel, Hmax = 0.005);



%% show geometry without mesh
ModelFig = figure(units="normalized",outerposition=[0 0.3 0.5 0.4])
    ModelPlot = pdegplot(sModel,FaceLabels="on", VertexLabels="on", FaceAlpha=0.3);
    title("Beam model")
    b=findobj(gca,'Type','Quiver');set(b,'Visible','off');    % remove coordinate-axes... they are often not nicely placed
    set(gca, 'XLim', [-0.5 +0.5]*L, 'YLim', [-1.1 1.1]*B, 'ZLim', [-1.1 1.1]*H)

%% use solve of PDE-Toolbox 
% RF = solve(sModel,FrequencyRange=[-Inf,12e4]);


%%  solve Eigenvalue-Problem using EIGS
if ~isempty(BC_faces_indx)   % --- with Dirichlet-BC: prescribed DoFs are removed --> Kc < K
    mat = assembleFEMatrices(sModel, 'nullspace');
    [U D]=eigs(mat.Kc,mat.M, 20, 'smallestabs');
else
    mat = assembleFEMatrices(sModel, 'MK');     % --- no Dirichlet-BC (free-body): full problem 
    [U D]=eigs(mat.K,mat.M, 20, 'smallestabs');
end;

EVP=struct;
    % .EVP      --> eigenvectors
    % .omega    --> eigenfrequencies
    % .NodalDisplacements --> displacements of nodes (=vertices)
    %   .ux,.uy, .uz      --> components of displacement
    %   .mag    --> absoulte value of displacement
    EVP.EV = U; EVP.omega = sqrt(diag(D));
    N = length(U(:,1))/3;
    EVP.NodalDisp.ux = U(    1:  N,:); 
    EVP.NodalDisp.uy = U(  N+1:2*N,:); 
    EVP.NodalDisp.uz = U(2*N+1:3*N,:); 
        % Dirichlet-BC:  for graphics add 0 to solution at nodes which have been dropped from problem due to condensation of matrices.
        BCindx = sort(findNodes(sModel.Geometry.Mesh, 'region', 'face', BC_faces_indx ));   % find node numbers where BC is imposed. here 0 must be added to solution vector for plotting!
        for i=1:length(indx)
            nn = length(U(1,:));
            EVP.NodalDisp.ux = [EVP.NodalDisp.ux(1:BCindx(i)-1,:); zeros(1,nn); EVP.NodalDisp.ux(BCindx(i):end,:)]; 
            EVP.NodalDisp.uy = [EVP.NodalDisp.uy(1:BCindx(i)-1,:); zeros(1,nn); EVP.NodalDisp.uy(BCindx(i):end,:)];
            EVP.NodalDisp.uz = [EVP.NodalDisp.uz(1:BCindx(i)-1,:); zeros(1,nn); EVP.NodalDisp.uz(BCindx(i):end,:)];
        end;   
    EVP.NodalDisp.mag = sqrt(EVP.NodalDisp.ux.^2 + EVP.NodalDisp.uy.^2 + EVP.NodalDisp.uz.^2);

   


%% Visualize
for ModeNr = 1:4

    meshData = sModel.Geometry.Mesh;
    nodalData = EVP.NodalDisp.mag(:,ModeNr);
    displacementData = [ EVP.NodalDisp.ux(:,ModeNr) ...
                            EVP.NodalDisp.uy(:,ModeNr) ...
                            EVP.NodalDisp.uz(:,ModeNr)];
 
    ColorLimits = [0 max(abs(nodalData))];
    dimensions = [L,B,H];    
    DefScaleFactor = 0.0007;

    %ResultFig = 
    figure(units="normalized",outerposition=[0.5 0.0 0.45 1]);
    VisualizeModes_3D_side_top(meshData, nodalData, displacementData, 1,  [-1/2 1/2]*L, [-1 1]*B/2, [-1 1]*H/2, DefScaleFactor, ColorLimits);




end;





