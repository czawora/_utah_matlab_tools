function shrunk_str = shrink_num_array(arr)
    
    sorted_arr = sort(arr);
    
    last_num = sorted_arr(1);
    last_written_num = sorted_arr(1);
    shrunk_str = sprintf('%d', sorted_arr(1));
    
    for i = 2:length(sorted_arr)
       
        if abs(last_num - sorted_arr(i)) > 1
            
            if last_num ~= last_written_num
                
                shrunk_str = sprintf('%s:%d %d', shrunk_str, last_num, sorted_arr(i));
            
            else
                
                shrunk_str = sprintf('%s %d', shrunk_str, sorted_arr(i));
            end
            
            last_written_num = sorted_arr(i);
        
        end
        
        last_num = sorted_arr(i);  
    end
   
    if last_num ~= last_written_num
        shrunk_str = sprintf('%s:%d', shrunk_str, last_num);
    end
    
end