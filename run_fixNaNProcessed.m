

paths = {'/Volumes/CZ_FRNU/micro_lfp/NIH029' ...
         '/Volumes/CZ_FRNU/micro_lfp/NIH030' ...
         '/Volumes/CZ_FRNU/micro_lfp/NIH039' ...
         '/Volumes/CZ_FRNU/micro_lfp/NIH042' ...
         '/Volumes/CZ_FRNU/micro_lfp/NIH059' ...
         '/Volumes/CZ_FRNU/micro_lfp/NIH062' ...
         '/Volumes/CZ_FRNU/micro_lfp/NIH064' ...
         '/Volumes/CZ_FRNU/micro_lfp/NIH066' ...
         '/Volumes/CZ_FRNU/micro_lfp/NIH069'};
     
for iPath = 1:length(paths)
    
   current_path = paths{iPath};
   path_ls = dir(current_path);
   
   for iDir = 1:length(path_ls)
       
      dir_name = path_ls(iDir).name;
      
      dir_path = [current_path '/' dir_name];
      processed_path = [ dir_path '/raw/*processed.mat' ];
      
      processed_ls = dir(processed_path);
      
      if ~isempty(processed_ls) 
          found_processed_fpath = [processed_ls.folder '/' processed_ls.name];
          fixNaNProcessed(found_processed_fpath);
      end
   end
    
end