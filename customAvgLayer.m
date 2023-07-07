classdef customAvgLayer < nnet.layer.Layer
    
    methods
        function layer=customAvgLayer(layerName)
            layer.Name=layerName;
        end
        
        function Y=predict(~,X)
            Y=avgpool(X,'global','DataFormat','SCBT','PoolFormat','T');
            if any(isnan(extractdata(Y)),'all')
                keyboard;
            end
        end
    end
end