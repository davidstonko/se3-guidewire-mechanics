function run_logged(fn, outfile)
%RUN_LOGGED  Run fn() and keep the command-window output, and the report of
%   any error, in the text file outfile. For long runs started remotely.
if exist(outfile, 'file'), fid = fopen(outfile, 'w'); fclose(fid); end
diary(outfile);  stop = onCleanup(@() diary('off'));
try
    fn();
catch err
    disp(getReport(err, 'extended', 'hyperlinks', 'off'));
end
fprintf('run_logged finished %s\n', char(datetime('now')));
end
