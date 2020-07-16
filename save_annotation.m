function save_annotation(testo,file_ann)
if(nargin<2),file_ann='error_log.txt';end
fprintf(' %s -> %s\n',datetime,testo);
fid=fopen(file_ann,'a');
fprintf(fid,'%s-> %s\n',datetime,testo);
fclose(fid);
end

