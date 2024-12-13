clear;

load('BeamExample.mat');
close all;


xvals = sModel.Geometry.Mesh.Nodes(1,:);
% DD = (xvals-min(xvals)).^2;     % quadr. 'stoerung'




f1=figure(units="normalized",outerposition=[0.1 0.0 0.45 1]);
f2=figure(units="normalized",outerposition=[0.6 0.3 0.35 0.3]);
for ModeNr = 1

    for k=50%:50
        scal = 0.02;
        offset = 25;
        
        meshData = sModel.Geometry.Mesh;
        nodalData = EVP.NodalDisp.mag(:,ModeNr);
        deformationData = [ EVP.NodalDisp.ux(:,ModeNr) ...
                            EVP.NodalDisp.uy(:,ModeNr) ...                            
                            EVP.NodalDisp.uz(:,ModeNr)+ ((offset-k)*scal*EVP.NodalDisp.uz(:,ModeNr+2))];   %EVP.NodalDisp.uz(:,ModeNr)+(25-k)*scal*DD.'];
 
    
    
        ColorLimits = [0 max(abs(nodalData))];
        dimensions = [L,B,H];    
        DefScaleFactor = 0.0007;
    
        
        figure(f1);
        VisualizeModes_3D_side_top(meshData, nodalData, deformationData, 1,  [-1/2 1/2]*L, [-1 1]*B/2, [-1 1]*H/2, DefScaleFactor, ColorLimits);
        VisualizeModes_3D_side_top(meshData, nodalData, deformationData, 1,  [-1/2 1/2]*L, [-1 1]*B/2, [-1 1]*H/2, DefScaleFactor, ColorLimits);
        drawnow;
    
        figure(2);
            indx=1:1:length(deformationData(:,1)); 
            %indx(BC_vertices_indx)=[];  % remove node indices with prescribed 0-displacement 
            indx(findNodes(sModel.Geometry.Mesh, 'region', 'vertex', BC_vertices_indx )) = [];
            indx = indx';
            u = [deformationData(indx,1); deformationData(indx,2); deformationData(indx,3) ];
            R(k) = (u.'*(mat.Kc*u))/(u.'*(mat.M*u));
        plot(offset-k, R(k), 'ko'); 
        hold on;
        pause(0.1);
    


    end;


end;