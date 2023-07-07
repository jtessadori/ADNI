classdef customConvLayer < nnet.layer.Layer
    properties
        N; % Number of filters in layer
        filterFormat;
        stride;
    end
    
    properties (Learnable)
        weights;
        bias;
    end
    
    methods
        function layer=customConvLayer(filterSize,nFilters,nFiltersPrev,FF,layerName,s)
            layer.Name=layerName;
            layer.N=nFilters;
            layer.filterFormat=FF;
            if ~exist('s','var')
                s=1;
            end
            layer.stride=s;
            
            % Initialize learnable properties
            layer.weights=randn([filterSize,nFiltersPrev,nFilters])/100;
            layer.bias=randn(1,nFilters)/100;
        end
        
        function Y=predict(layer,X)
            Y=dlconv(X,dlarray(layer.weights,layer.filterFormat),layer.bias,'Stride',layer.stride,'DataFormat','SCBT');
            if any(isnan(extractdata(Y)),'all')
                keyboard;
            end
        end
    end
end