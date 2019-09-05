function convert_double2int16()

test = 1;
rootPath = '/Volumes/56PROC/private/micro_processing';

subjs = {'NIH029' ...
         'NIH030' ...
         'NIH034' ...
         'NIH036' ...
         'NIH037' ...
         'NIH039' ...
         'NIH042' ...
         'NIH046' ...
         'NIH047' ...
         'NIH048' ...
         'NIH049' ...
         'NIH050' ...
         'NIH051' ...
         'NIH052' ...
         'NIH053' ...
         'NIH054' ...
         'NIH057' ...
         'NIH059' ...
         'NIH060' ...
         'NIH062' ...
         'NIH063' ...
         'NIH064' ...
         'NIH066' ...
         'NIH067' ...
         'NIH069' };
     
for iSubj = 1:length(subjs)
    
    subj_ms_path = [rootPath '/' subjs{iSubj} '/mountain_sorts'];
        
    fprintf('looking in %s\n', subj_ms_path);

    % check for mountainsorts folder
    if exist(subj_ms_path, 'dir')
    
        session_ls = dir(subj_ms_path);
                
        for iSess = 1:length(session_ls)
            
            current_sess = session_ls(iSess).name;
            
            fprintf('\t%s\n', current_sess);
            
            sess_noreref_path = [subj_ms_path '/' current_sess '/raw/*noreref.mat'];
            sess_processed_path = [subj_ms_path '/' current_sess '/raw/*processed.mat'];
            
            
            
            
            sess_noreref_ls = dir(sess_noreref_path);
            
            if ~isempty(sess_noreref_ls)   
                convert_double2int16_([sess_noreref_ls.folder '/' sess_noreref_ls.name], test);
            end
            
            
                        
            sess_processed_ls = dir(sess_processed_path);
            
            if ~isempty(sess_processed_ls)   
                convert_double2int16_([sess_processed_ls.folder '/' sess_processed_ls.name], test);
            end
            
            
            
        end
        
    end
    
end

end

function convert_double2int16_(mat_path, test)

    fprintf('converting %s\n', mat_path);
    
    if ~exist(mat_path, 'file')
       keyboard;
       error('looking for file that does not exist'); 
    end
    
    if ~test

        lfpStruct = load(mat_path);
        lfpStruct = lfpStruct.lfpStruct;

        for iCell = 1:length(lfpStruct.lfp)

            if isa(lfpStruct.lfp{iCell}, 'double')
                lfpStruct.lfp{iCell} = int16(lfpStruct.lfp{iCell});
            end
        end

        save(mat_path,  '-v7.3', 'lfpStruct');
        
    end

end