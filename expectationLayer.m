classdef expectationLayer < nnet.layer.Layer
    methods
        function layer=expectationLayer(layerName)
            layer.Name=layerName;
            layer.NumOutputs=2;
        end
        
        function [T1,T2]=predict(~,clustPDF)
            % Estimate probability of being in each cluster from clustPDFs
            clustP=single(clustPDF./(sum(clustPDF)+eps));

            % Compute TR explicitly from clustP
            softMaxFun=@(x)exp(x)./sum(exp(x));
            clustP=softMaxFun(clustP);
            currTR=squeeze(mean(permute(clustP(:,:,1:end-1),[4,1,2,3]).*permute(clustP(:,:,2:end),[1,4,2,3]),[3,4]));
            currTR=currTR./sum(currTR);

            % Apply Viterbi
            T1=dlarray(zeros(size(clustPDF)));
            T2=dlarray(zeros(size(clustPDF)));
            T1(:,:,1)=log(clustPDF(:,:,1)/size(clustPDF,1)+eps);
            for currT=2:size(T1,3)
                currOut=T1(:,:,currT-1)+permute(log(currTR+eps),[1,3,2])+log(clustPDF(:,:,currT)+eps);
                [T1(:,:,currT),T2(:,:,currT)]=max(currOut,[],3);
            end
        end
    end
end