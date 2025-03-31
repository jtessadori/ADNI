classdef GRBFlayer < nnet.layer.Layer
    properties
        N; % Number of neurons in layer
        epsilon;
    end
    
    properties (Learnable)
        c;
        sigma;
    end
    
    methods
        function layer=GRBFlayer(siz,inputDims,layerName)
            % Construct a generalized radial basis function neuron layer
            %(http://www.cs.cmu.edu/~fmetze/interACT/Publications_files/publications/grbf.pdf)
            layer.Name=layerName;
            layer.N=siz;
            layer.NumOutputs=2;
            
            % Initialize learnable properties
            layer.c=dlarray(randn(inputDims,layer.N,'single'));
%             layer.c=dlarray(randn(inputDims,layer.N,'single')*sqrt(6)/(sqrt(inputDims+layer.N)));
            layer.sigma=dlarray(ones(layer.N,1,'single')*sqrt(inputDims));
%             layer.sigma=dlarray(ones(layer.N,1,'single'));
%             layer.sigma=dlarray(ones(layer.N,1,'single'));
            layer.epsilon=1e-10;
        end
        
        function [clustPDF,clusterP]=predict(layer,X1)
%             D=sqrt(sum((X1-layer.c').^2,2)+eps);
%             clusterP=dlarray(zeros(length(layer.sigma),1,size(X1,3),size(X1,4)));
%             for currN=1:length(layer.sigma)
%                 clusterP(currN,:,:,:)=normpdf(D(currN,:,:,:),0,max(eps,layer.sigma(currN)));
%             end
            D=sqrt(sum((X1-permute(layer.c,[1,3,4,2])).^2)+eps);
            clusterP=dlarray(zeros(length(layer.sigma),size(X1,2),size(X1,3)));
            for currN=1:length(layer.sigma)
                clusterP(currN,:,:)=normpdf(D(:,:,:,currN),0,max(eps,layer.sigma(currN)));
            end
            clusterP(clusterP<layer.epsilon)=0;
            clusterP(:,sum(abs(clusterP))==0)=eps;
            clustPDF=clusterP;
            clusterP=single(clusterP./(sum(clusterP)+eps));
            if any(isnan(extractdata(clusterP)),'all')
                keyboard;
            end
        end
    end
end