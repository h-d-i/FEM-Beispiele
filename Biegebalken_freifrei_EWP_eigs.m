% -----------------------------------------
% Lineare Schwingungen
% 
% Bsp.: Balken mit "frei-frei"-Randbedingungen
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
sModel.MaterialProperties = materialProperties(YoungsModulus=E, ...
                                              PoissonsRatio=nu, ...
                                              MassDensity=rho);
% Boundary Conditions
    BC_vertices_indx = [];% rechtes Ende: BC_vertices_indx =[2 3 5 8];   % rechtes Ende: BC_vertices_indx = [1 4 6 7]
    %BC_vertices_indx =[2 3 5 8];
    sModel.VertexBC(BC_vertices_indx) = vertexBC(Constraint="fixed");

% create mesh
sModel = generateMesh(sModel, Hmax = 0.005);



%% show geometry without mesh
ModelFig = figure(units="normalized",outerposition=[0 0.3 0.5 0.4])
    ModelPlot = pdegplot(sModel,FaceLabels="on", VertexLabels="on", FaceAlpha=0.3);
    title("Beam model")
    b=findobj(gca,'Type','Quiver');set(b,'Visible','off');    % remove coordinate-axes... they are often not nicely placed
    set(gca, 'XLim', [-0.5 +0.5]*L, 'YLim', [-1.1 1.1]*B, 'ZLim', [-1.1 1.1]*H)

%% use solve of PDE-Toolbox 
%RF = solve(sModel,FrequencyRange=[-Inf,12e4]);


%%  solve Eigenvalue-Problem using EIGS
if ~isempty(BC_vertices_indx)   % --- with Dirichlet-BC: prescribed DoFs are removed --> Kc < K
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
        % Dirichlet-BC:  add 0 to solution at nodes which have been dropped from problem due to condensation of matrices.
        indx = findNodes(sModel.Geometry.Mesh, 'region', 'vertex', BC_vertices_indx );   % find node numbers where BC is imposed. here 0 must be added to solution vector for plotting!
        for i=1:length(indx)
            nn = length(U(1,:));
            EVP.NodalDisp.ux = [EVP.NodalDisp.ux(1:indx(i)-1,:); zeros(1,nn); EVP.NodalDisp.ux(indx(i):end,:)]; 
            EVP.NodalDisp.uy = [EVP.NodalDisp.uy(1:indx(i)-1,:); zeros(1,nn); EVP.NodalDisp.uy(indx(i):end,:)];
            EVP.NodalDisp.uz = [EVP.NodalDisp.uz(1:indx(i)-1,:); zeros(1,nn); EVP.NodalDisp.uz(indx(i):end,:)];
        end;   
    EVP.NodalDisp.mag = sqrt(EVP.NodalDisp.ux.^2 + EVP.NodalDisp.uy.^2 + EVP.NodalDisp.uz.^2);

   


%% Visualize
for ModeNr = 1:6

        
        meshData = sModel.Geometry.Mesh;
        nodalData = EVP.NodalDisp.mag(:,ModeNr);
        displacementData = [ EVP.NodalDisp.ux(:,ModeNr) ...
                            EVP.NodalDisp.uy(:,ModeNr) ...
                            EVP.NodalDisp.uz(:,ModeNr)];
 
    ColorLimits = [0 max(abs(nodalData))];
    %dimensions = [L,B,H];    
    DefScaleFactor = 0.0007;

    %ResultFig = 
    figure(units="normalized",outerposition=[0.5 0.0 0.45 1]);
    VisualizeModes_3D_side_top(meshData, nodalData, displacementData, 1,  [-1/2 1/2]*L, [-1 1]*B/2, [-1 1]*H/2, DefScaleFactor, ColorLimits)




end;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % Animation
    % gifFile = ['Balken-eingespannt--ModeNr=',num2str(ModeNr),'.gif'];
    % exportgraphics(ResultFig, gifFile);
    % 
    % for ii = 0:0.01:2
    %     phaseData = cospi(ii) + 1i*sinpi(ii);   % exp(j alpha) 
    %     resultViz.NodalData = abs(real(nodalData*phaseData));   %--> colors
    %     resultViz.DeformationData = deformationData*phaseData;
    %     pause(0.01)
    %     exportgraphics(ResultFig, gifFile, Append=true);
    % end




    %     ResultFig = figure(units="normalized",outerposition=[0.5 0.0 0.5 1]);
%     meshData = RF.Mesh;
%     nodalData = RF.ModeShapes.Magnitude(:,ModeNr);
%     deformationData = [ RF.ModeShapes.ux(:,ModeNr) ...
%                         RF.ModeShapes.uy(:,ModeNr) ...
%                         RF.ModeShapes.uz(:,ModeNr)];
%     phaseData = cospi(0) + 1i*sinpi(0);   
%     DefScaleFactor = 0.0007;
%     XLimits = [-L/2 L/2]; YLimits = [-1 1]*3*B; ZLimits = [-1 1]*3*H;
% 
% 
% 
% % Create PDE result visualization
% 
% t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact'); 
% nexttile
%     resultViz = pdeviz(meshData, abs(real(nodalData*phaseData)), ...
%             "MeshVisible", 'on', ...
%             "DeformationData",deformationData*phaseData, ...    % --> displacement
%             "DeformationScaleFactor",DefScaleFactor, ...                
%             "AxesVisible",true, ...
%             "XLabel", 'x', "YLabel", 'y', "ZLabel", 'z', ...
%             "ColorbarVisible",false, ...
%             "ColorLimits",[0 19.48] ...
%         );
%         resultViz.XLimits = XLimits;
%         resultViz.YLimits = YLimits;
%         resultViz.ZLimits = ZLimits;
% 
% nexttile
%      resultViz2 = pdeviz(meshData, abs(real(nodalData*phaseData)), ...
%             "MeshVisible", 'on', ...
%             "DeformationData",deformationData*phaseData, ...    % --> displacement
%             "DeformationScaleFactor",DefScaleFactor, ...                
%             "AxesVisible",true, ...
%             "XLabel", 'x',   "ZLabel", 'z', ...
%             "ColorbarVisible",false, ...
%             "ColorLimits",[0 19.48] ...
%         );
%         resultViz2.XLimits = XLimits;
%         resultViz2.YLimits = YLimits;
%         resultViz2.ZLimits = ZLimits;
%         resultViz2.View = [0 0];
% 
% nexttile
%      resultViz3 = pdeviz(meshData, abs(real(nodalData*phaseData)), ...
%             "MeshVisible", 'on', ...
%             "DeformationData",deformationData*phaseData, ...    % --> displacement
%             "DeformationScaleFactor",DefScaleFactor, ...                
%             "AxesVisible",true, ...
%             "XLabel", 'x', "YLabel", 'y',   ...
%             "ColorbarVisible",false, ...
%             "ColorLimits",[0 19.48] ...
%         );
%         resultViz3.XLimits = XLimits;
%         resultViz3.YLimits = YLimits;
%         resultViz3.ZLimits = ZLimits;
%         resultViz3.View = [0 90];

