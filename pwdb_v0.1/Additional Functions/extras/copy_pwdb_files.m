function copy_pwdb_files(pwdb_no)

PATHS = setup_paths_for_post_processing(pwdb_no);

source_dir = '/Volumes/pete_data/201901_PWDB_complete/';
source_dir = '/Volumes/pete_data/201902_PWDB_initial/';
req_domains = [1,2,3,4,7,14,15,16,18,19,21,22,27,28,35,37,39,41,42,44,46,49,65,71,72,79,84,87,96,108,112];
req_domains = [3,4];
doing_extra_files = 1;

% copy input files
max_subjs = 21;
for subj_no = 1 : max_subjs
    
    fprintf(['\n ' num2str(subj_no)])
    
    % check that this subject exists
    curr_filepath = [source_dir, 'sim_' num2str(subj_no) '.in'];
    if ~exist(curr_filepath, 'file')
        break
    end
    
    % copy input files
    if ~doing_extra_files
        file_exts = {['sim_' num2str(subj_no) '.mat'];
            ['sim_' num2str(subj_no) '.in'];
            ['sim_' num2str(subj_no) '_IN_1.bcs']};
        for file_no = 1 : length(file_exts)
            in_file = [source_dir, file_exts{file_no}];
            out_file = [PATHS.InputFiles, file_exts{file_no}];
            copyfile(in_file, out_file);
        end
    end
    
    % copy output files for each required domain
    for dom_no = 1 : length(req_domains)
        curr_dom = req_domains(dom_no);
        file_exts = {['sim_' num2str(subj_no) '_' num2str(curr_dom) '.his'], ...
            ['sim_' num2str(subj_no) '_out_' num2str(curr_dom) '.lum']};
        for file_no = 1 : length(file_exts)
            in_file = [source_dir, file_exts{file_no}];
            out_file = [PATHS.OutputFiles, file_exts{file_no}];
            if exist(in_file, 'file')
                method = 'rsync'; % or 'copyfile'
                do_copy(in_file, out_file, method);
            end
        end
    end
    
    % copy remaining output files
    if ~doing_extra_files
        file_exts = {['sim_' num2str(subj_no) '_period.tex']};
        for file_no = 1 : length(file_exts)
            in_file = [source_dir, file_exts{file_no}];
            out_file = [PATHS.OutputFiles, file_exts{file_no}];
            copyfile(in_file, out_file);
        end
    end
    
end

end

function do_copy(in_file, out_file, method)

if strcmp(method, 'rsync')
    command = ['rsync ' in_file ' ' , out_file];
    [status,cmdout] = system(command,'-echo');
elseif strcmp(method, 'copyfile')
    copyfile(in_file, out_file);
end

end