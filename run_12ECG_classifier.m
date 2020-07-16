function [score, label,classes] = run_12ECG_classifier(data,header_data, loaded_model)


	model=loaded_model.model;
	classes=loaded_model.classes;

    num_classes = length(classes);

    label = zeros([1,num_classes]);
    score = ones([1,num_classes]);
  [out_scores,out_labels]=get_12ECG_cls_ivo(data,header_data,model,classes);
  
    % Use your classifier here to obtain a label and score for each class.
    %features = get_12ECG_features(data,header_data);

    
   % score = mnrval(model,features);		
%     [~,idx] = max (score);

%     label(idx)=1;
    score=out_scores;
    label=out_labels;
    
    
end



