% ==========================================================================================================================================================
% Lineare Schwingungen
% 
% Bsp.: Turbinenschaufel 
% ====================================
% 
% Berechnung Eigenfrequenzen + Moden
% Visualisierung
% 
% -----------------------------------------
%
% Geometrie: Beispiel aus Matlab PDEToolbox
%
% Quellen: 
% verschiedene... u.a.
% "Structural Dynamics af Tuning Fork" --> https://de.mathworks.com/help/pde/ug/structural-dynamics-of-tuning-fork.html
% ==========================================================================================================================================================

clear; close all;

%% material parameters 
E = 210e9; nu = 0.3; rho = 8000;


%% define FE-problem
% ... geometry ..............................
gm = importGeometry('Blade.stl');
    rotate(gm, -90);               % turn -90° about z-axis 
    x0 = min(gm.Vertices(:,1));    % min x-values = left end of blade
    y0 = mean(gm.Vertices(:,2));   % mean y-values = center line of blade
    translate(gm, -[x0,y0, 0.12]); % translate
L =  max(gm.Vertices(:,1)) - min(gm.Vertices(:,1)); 
B =  max(gm.Vertices(:,2)) - min(gm.Vertices(:,2)); 
H =  max(gm.Vertices(:,3)) - min(gm.Vertices(:,3)); 


% ... setup fe-model
sModel = femodel(AnalysisType='structuralModal', Geometry=gm);
sModel.MaterialProperties = materialProperties(YoungsModulus=E, ...
                                              PoissonsRatio=nu, ...
                                              MassDensity=rho);
% ...Boundary Conditions .............................. 
    % ... on faces
    BC_faces_indx = [3];    % face no. 3 is fixed
    sModel.FaceBC(BC_faces_indx) = faceBC(Constraint="fixed");

% ... create mesh ..............................
sModel = generateMesh(sModel);


%% use solver of PDE-Toolbox --- not used here to have more control over format of solution!
% RF = solve(sModel,FrequencyRange=[-Inf,12e4]);


%%  solve Eigenvalue-Problem using EIGS
if ~isempty(BC_faces_indx)    % --- with Dirichlet-BC: restricted DoFs are removed --> Kc < K
    mat = assembleFEMatrices(sModel, 'nullspace');
    [U D]=eigs(mat.Kc,mat.M, 10, 'smallestabs');
else
    mat = assembleFEMatrices(sModel, 'MK');                % --- no Dirichlet-BC (free-body): full problem 
    [U D]=eigs(mat.K,mat.M, 10, 'smallestabs');
end;

% -- save EVP-data ........... and determine nodal displacement vector from it 
EVP=struct;
        % .EVP      --> eigenvectors
        % .omega    --> eigenfrequencies
        % .NodalDisplacements --> displacements of nodes (=vertices)
        %   .ux,.uy, .uz      --> components of displacement
        %   .mag    --> absoulte value of displacement
    EVP.EV = U; 
    EVP.omega = sqrt(diag(D));
    % .... determine nodal displacement vectors from eigenvectors ............................
    N = length(U(:,1))/3;
    EVP.NodalDisp.ux = U(    1:  N,:); 
    EVP.NodalDisp.uy = U(  N+1:2*N,:); 
    EVP.NodalDisp.uz = U(2*N+1:3*N,:); 
        % Dirichlet-BC:  add 0 to solution at nodes which have been dropped from problem due to condensation of matrices.
        %indxV = findNodes(sModel.Geometry.Mesh, 'region', 'vertex', BC_vertices_indx );   % find node numbers where BC is imposed. here 0 must be added to solution vector for plotting!
        indx = sort(findNodes(sModel.Geometry.Mesh, 'region', 'face', BC_faces_indx ));   % find node numbers where BC is imposed. here 0 must be added to solution vector for plotting!
        % indx = [indxV, indxF];
        for i=1:length(indx)
            nn = length(U(1,:));
            EVP.NodalDisp.ux = [EVP.NodalDisp.ux(1:indx(i)-1,:); zeros(1,nn); EVP.NodalDisp.ux(indx(i):end,:)]; 
            EVP.NodalDisp.uy = [EVP.NodalDisp.uy(1:indx(i)-1,:); zeros(1,nn); EVP.NodalDisp.uy(indx(i):end,:)];
            EVP.NodalDisp.uz = [EVP.NodalDisp.uz(1:indx(i)-1,:); zeros(1,nn); EVP.NodalDisp.uz(indx(i):end,:)];
        end;   
    EVP.NodalDisp.mag = sqrt(EVP.NodalDisp.ux.^2 + EVP.NodalDisp.uy.^2 + EVP.NodalDisp.uz.^2);

   
%% PLOT RESULTS


%% Visualize 

%plot geometry without mesh
ModelFig = figure(units="normalized",outerposition=[0 0.5 0.3 0.4], color = 'white');
    ModelPlot = pdegplot(sModel,FaceLabels="off", VertexLabels="off", FaceAlpha=0.3);
    title("Turbine Blade");
    b=findobj(gca,'Type','Quiver');set(b,'Visible','off');    % remove coordinate-axes... they are often not nicely placed
    set(gca, 'XLim', [0 1]*L, 'YLim', [-1.1 1.1]*B/2, 'ZLim', [-1.1 1.1]*H/2);
    exportgraphics(ModelFig, "TurbineBlade.gif");

% plot mesh    
MeshFig = figure(units="normalized",outerposition=[0 0.1 0.3 0.4], color = 'white');
    MeshPlot = pdeplot3D(sModel.Geometry.Mesh);
    b=findobj(gca,'Type','Quiver');set(b,'Visible','off');    % remove coordinate-axes... they are often not nicely placed
    set(MeshPlot, 'FaceColor', [1 1 1]*0.9, 'FaceAlpha', 0.8, 'LineWidth', 1, 'EdgeColor', 'b', 'MarkerSize', 3, 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r') ;       
    set(gca, 'XLim', [0 1]*L, 'YLim', [-1.1 1.1]*B/2, 'ZLim', [-1.1 1.1]*H/2);
    exportgraphics(MeshFig, "TurbineBlade_mesh.gif");


% plot eigenfrequencies
FrequFig = figure(units="normalized",outerposition=[0.1 0.05 0.4 0.5], color = 'white');
    frequplot = plot(1:length(EVP.omega), real(EVP.omega/2/pi), 'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
    grid on;  set(gca, 'XLim', [0 10], 'ylim', [0 max(real(EVP.omega/2/pi))]);
    xlabel('Mode Nr.'); ylabel('Eigenfrequenz / Hz');
    exportgraphics(FrequFig, "TurbineBlade_Frequencies.gif");

% plot modes
B = max(B,H);
H = B;
for ModeNr=1:10
   
   % .... data for figures ........................
   meshData = sModel.Geometry.Mesh;
   nodalData = EVP.NodalDisp.mag(:,ModeNr);
   displacementData = [ EVP.NodalDisp.ux(:,ModeNr) ...
                            EVP.NodalDisp.uy(:,ModeNr) ...
                            EVP.NodalDisp.uz(:,ModeNr)];
   % .... graphics parameters ....................
   DefScaleFactor = 0.035;
   ColorLimits = [0 max(abs(nodalData))];
   phaseData = 1;
   scaledmaxDisplacements = DefScaleFactor*max(abs(displacementData), [], 1);    % max. Displacements in x,y,z - direction
   Nanimsteps =  20;    % number of pictures for animation.... NO ANIMATION for =0
   AnimCycleTime = 3;       % duration of one Animation cycle --> to calculate delay between .gif - figures

   % .... PLOT ............................
   ResultFig = figure(units="normalized",outerposition=[0.6 0.1 0.35 0.85], color='white');
   t=VisualizeModes_3D_side_top(meshData, nodalData*phaseData, displacementData*phaseData, 1, [0 1]*(L+scaledmaxDisplacements(1)), [-1 1]*(B/2+scaledmaxDisplacements(2)), [-1 1]*(H/2 + scaledmaxDisplacements(2)), DefScaleFactor, ColorLimits);
   t.Title.String = [num2str(round(EVP.omega(ModeNr)/2/pi), '%.0f'), ' Hz'];
   

   % ... Animation --> save to animated .gif ................................
   gifFile = ['TurbineBlade--ModeNr=',num2str(ModeNr),'.gif']; 
   if Nanimsteps >0
       for ii = 1:Nanimsteps     
            phaseData = cos(2*pi/Nanimsteps *(ii-1)) + 1i*sin(2*pi/Nanimsteps *(ii-1));   % exp(j alpha) 
            t=VisualizeModes_3D_side_top(meshData, nodalData*phaseData, displacementData*phaseData, 1, [0 1]*(L+scaledmaxDisplacements(1)), [-1 1]*(B/2+scaledmaxDisplacements(2)), [-1 1]*(H/2 + scaledmaxDisplacements(2)), DefScaleFactor, ColorLimits);
            t.Title.String = [num2str(round(EVP.omega(ModeNr)/2/pi), '%.0f'), ' Hz'];
            set(gcf, units="normalized",outerposition=[0.6 0.1 0.35 0.85]);
            drawnow        
            frame = getframe(gcf);
            im{ii} = frame2im(frame);       
       end
       for idx = 1:numel(im)
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                if exist(gifFile, 'file')==2
                    delete(gifFile);
                end
                imwrite(A,map,gifFile,"gif",LoopCount=Inf, DelayTime=AnimCycleTime/Nanimsteps)
            else
                imwrite(A,map,gifFile,"gif",WriteMode="append", DelayTime=AnimCycleTime/Nanimsteps)
            end
       end;   
   end;
end;







