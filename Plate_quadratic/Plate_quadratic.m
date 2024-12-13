% ===============================================================================
% Lineare Schwingungen
% 
% Bsp.: quadratische Platte (frei-frei)
% ====================================
% 
% Berechnung Eigenfrequenzen + Moden
% Visualisierung
% 
% -----------------------------------------

% Quellen
% verschiedene... u.a.
% "Structural Dynamics af Tuning Fork" --> https://de.mathworks.com/help/pde/ug/structural-dynamics-of-tuning-fork.html
% ===============================================================================

clear; close all;

%% parameters of the geometry
L = 0.2;
B = 0.2;
H = 0.005;
E = 210e9; nu = 0.3; rho = 8000;


%% define FE-problem
% geometry
gm = multicuboid(L,B,H, Zoffset=-H/2);

% setup fe-model
sModel = femodel(AnalysisType='structuralModal', Geometry=gm);
sModel.MaterialProperties = materialProperties(YoungsModulus=E, PoissonsRatio=nu, MassDensity=rho);

% create mesh
sModel = generateMesh(sModel, Hmax = 0.01);



% %% show geometry without mesh
% ModelFig = figure(units="normalized",outerposition=[0 0.3 0.5 0.4])
%     ModelPlot = pdegplot(sModel,FaceLabels="on", VertexLabels="on", FaceAlpha=0.3);
%     title("Beam model")
%     b=findobj(gca,'Type','Quiver');set(b,'Visible','off');    % remove coordinate-axes... they are often not nicely placed
%     set(gca, 'XLim', [-0.5 +0.5]*L, 'YLim', [-1.1 1.1]*B, 'ZLim', [-1.1 1.1]*H)

%% use solve of PDE-Toolbox 
%RF = solve(sModel,FrequencyRange=[-Inf,12e4]);


%%  solve Eigenvalue-Problem using EIGS
    mat = assembleFEMatrices(sModel, 'MK');     % --- no Dirichlet-BC (free-body): full problem 
    [U D]=eigs(mat.K,mat.M, 30, 'smallestabs');


    EVP=struct;
    % .EVP      --> eigenvectors
    % .omega    --> eigenfrequencies
    % .NodalDisplacements --> displacements of nodes (=vertices)
    %   .ux,.uy, .uz      --> components of displacement
    %   .mag    --> absoulte value of displacement
    EVP.EV = U; 
    EVP.omega = sqrt(diag(D));
    N = length(U(:,1))/3;
    EVP.NodalDisp.ux = U(    1:  N,:); 
    EVP.NodalDisp.uy = U(  N+1:2*N,:); 
    EVP.NodalDisp.uz = U(2*N+1:3*N,:); 
    EVP.NodalDisp.mag = sqrt(EVP.NodalDisp.ux.^2 + EVP.NodalDisp.uy.^2 + EVP.NodalDisp.uz.^2);

   


%% Visualize 
%plot geometry without mesh
ModelFig = figure(units="normalized",outerposition=[0 0.5 0.5 0.4], color = 'white');
    ModelPlot = pdegplot(sModel,FaceLabels="off", VertexLabels="off", FaceAlpha=0.3);
    title("Quadratic Plate");
    b=findobj(gca,'Type','Quiver');set(b,'Visible','off');    % remove coordinate-axes... they are often not nicely placed
    %set(gca, 'XLim', [-0.5 +0.5]*L, 'YLim', [-1.1 1.1]*B, 'ZLim', [-1.1 1.1]*H);
    set(gca, 'XLim', [-0.5 +0.5]*L, 'YLim', [-0.5 0.5]*B, 'ZLim', [-0.5 0.5]*H)
    exportgraphics(ModelFig, "FreeQuadraticPlate.gif");

% plot mesh    
MeshFig = figure(units="normalized",outerposition=[0 0.1 0.5 0.6], color = 'white');
    MeshPlot = pdeplot3D(sModel.Geometry.Mesh);
    b=findobj(gca,'Type','Quiver');set(b,'Visible','off');    % remove coordinate-axes... they are often not nicely placed
    set(MeshPlot, 'FaceColor', [1 1 1]*0.9, 'FaceAlpha', 0.8, 'LineWidth', 1, 'EdgeColor', 'b', 'MarkerSize', 3, 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r') ;       
    set(gca, 'XLim', [-0.5 +0.5]*L, 'YLim', [-0.5 0.5]*B, 'ZLim', [-0.5 0.5]*H)
    exportgraphics(MeshFig, "FreeQuadraticPlate_mesh.gif");


% plot eigenfrequencies
FrequFig = figure(units="normalized",outerposition=[0.1 0.05 0.5 0.5], color = 'white');
    frequplot = plot(1:length(EVP.omega), real(EVP.omega/2/pi), 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
    grid on;  set(gca, 'XLim', [0 21]);
    xlabel('Mode Nr.'); ylabel('Eigenfrequenz / Hz');
    exportgraphics(FrequFig, "FreeQuadraticPlate_Frequencies.gif");

% plot modes
for ModeNr = 1:20
   
   % .... data for figures ........................
   meshData = sModel.Geometry.Mesh;
   nodalData = EVP.NodalDisp.mag(:,ModeNr);
   displacementData = [ EVP.NodalDisp.ux(:,ModeNr) ...
                            EVP.NodalDisp.uy(:,ModeNr) ...
                            EVP.NodalDisp.uz(:,ModeNr)];
   % .... graphics parameters ....................
   DefScaleFactor = 0.02;
   ColorLimits = [0 max(abs(nodalData))];
   phaseData = 1;
   scaledmaxDisplacements = DefScaleFactor*max(abs(displacementData), [], 1);    % max. Displacements in x,y,z - direction
   Nanimsteps = 20;    % number of pictures for animation.... NO ANIMATION for =0
   AnimCycleTime = 3;       % duration of one Animation cycle --> to calculate delay between .gif - figures

   % .... PLOT ............................
   ResultFig = figure(units="normalized",outerposition=[0.5 0.0 0.35 1]);   
   t=VisualizeModes_3D_side_top(meshData, nodalData*phaseData, displacementData*phaseData, 1, [-1 1]*(L/2+scaledmaxDisplacements(1)), [-1 1]*(B/2+scaledmaxDisplacements(2)), [-1 1]*(H/2 + scaledmaxDisplacements(3)), DefScaleFactor, ColorLimits);
   t.Title.String = [num2str(round(EVP.omega(ModeNr)/2/pi), '%.0f'), ' Hz'];

   % ... Animation --> save to animated .gif ................................
   gifFile = ['Plate--ModeNr=',num2str(ModeNr),'.gif']; 
   if Nanimsteps >0
       for ii = 1:Nanimsteps     
            phaseData = cos(2*pi/Nanimsteps *(ii-1)) + 1i*sin(2*pi/Nanimsteps *(ii-1));   % exp(j alpha) 
            t=VisualizeModes_3D_side_top(meshData, nodalData*phaseData, displacementData*phaseData, 1, [-1 1]*(L/2+scaledmaxDisplacements(1)), [-1 1]*(B/2+scaledmaxDisplacements(2)), [-1 1]*(H/2 + scaledmaxDisplacements(3)), DefScaleFactor, ColorLimits);
            t.Title.String = [num2str(round(EVP.omega(ModeNr)/2/pi), '%.0f'), ' Hz'];
            set(gcf, units="normalized",outerposition=[0.5 0.0 0.35 1]);
            drawnow        
            frame = getframe(gcf);
            im{ii} = frame2im(frame);       
       end
       for idx = 1:numel(im)
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,gifFile,"gif",LoopCount=Inf, DelayTime=AnimCycleTime/Nanimsteps)
            else
                imwrite(A,map,gifFile,"gif",WriteMode="append", DelayTime=AnimCycleTime/Nanimsteps)
            end
       end;   
   end;

end;





