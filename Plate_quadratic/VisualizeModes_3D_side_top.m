function   t=VisualizeModes_3D_side_top(meshData, nodalData, deformationData, phasor, XLimits, YLimits, ZLimits,   DefScaleFactor, ColorLimits)



% ResultFig = figure(units="normalized",outerposition=[0.5 0.0 0.5 1]);
    %meshData = RF.Mesh;
    %nodalData = RF.ModeShapes.Magnitude(:,ModeNr);
    %deformationData = [ RF.ModeShapes.ux(:,ModeNr) ...
    %                    RF.ModeShapes.uy(:,ModeNr) ...
    %                    RF.ModeShapes.uz(:,ModeNr)];
    % phaseData = cospi(0) + 1i*sinpi(0);   
    phaseData = phasor;
    %DefScaleFactor = 0.0007;



% Create PDE result visualization

t = tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact'); 

%    title('tst')

    
nexttile

    
    resultViz = pdeviz(meshData, abs(real(nodalData*phaseData)), ...
            "MeshVisible", 'on', ...
            "DeformationData",deformationData*phaseData, ...    % --> displacement
            "DeformationScaleFactor",DefScaleFactor, ...                
            "AxesVisible",true, ...
            "XLabel", 'x', "YLabel", 'y', "ZLabel", 'z', ...
            "ColorbarVisible",false, ...
            "ColorLimits",ColorLimits ...
        );
        resultViz.XLimits = XLimits;
        resultViz.YLimits = YLimits;
        resultViz.ZLimits = ZLimits;
        % resultViz.View = angles;
    % hold on;
    % plot3(1,2,3, 'ro')
 
    
   
nexttile
     resultViz2 = pdeviz(meshData, abs(real(nodalData*phaseData)), ...
            "MeshVisible", 'on', ...
            "DeformationData",deformationData*phaseData, ...    % --> displacement
            "DeformationScaleFactor",DefScaleFactor, ...                
            "AxesVisible",true, ...
            "XLabel", 'x',   "ZLabel", 'z', ...
            "ColorbarVisible",false, ...
            "ColorLimits",ColorLimits ...
        );
        resultViz2.XLimits = XLimits;
        resultViz2.YLimits = YLimits;
        resultViz2.ZLimits = ZLimits;
        resultViz2.View = [0 0];

nexttile
     resultViz3 = pdeviz(meshData, abs(real(nodalData*phaseData)), ...
            "MeshVisible", 'on', ...
            "DeformationData",deformationData*phaseData, ...    % --> displacement
            "DeformationScaleFactor",DefScaleFactor, ...                
            "AxesVisible",true, ...
            "XLabel", 'x', "YLabel", 'y',   ...
            "ColorbarVisible",false, ...
            "ColorLimits",ColorLimits ...
        );
        resultViz3.XLimits = XLimits;
        resultViz3.YLimits = YLimits;
        resultViz3.ZLimits = ZLimits;
        resultViz3.View = [0 90];
