function value = pipeline_config_value(key)
% pipeline_config_value  One path from config.yaml, via pipeline_config.py
%   value = pipeline_config_value('pvs_root')
% Calls python3, which must be on the system path.

    repo_dir = fileparts(mfilename('fullpath'));
    helper = fullfile(repo_dir, 'pipeline_config.py');
    [status, out] = system(sprintf('python3 "%s" %s 2>/dev/null', helper, key));
    if status ~= 0
        error('pipeline_config.py could not return "%s" (status %d). Run it in a shell to see why.', key, status);
    end
    value = strtrim(out);
end
