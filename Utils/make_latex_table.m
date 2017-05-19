function [ varargout ] = make_latex_table( T, varargin )
%MAKE_LATEX_TABLE Create a Latex table environment
%   TABLE_STR = MAKE_LATEX_TABLE( T ) Returns a string representation of
%   T that includes \begin{table} and \begin{tabular} (and corresponding
%   \end commands) environments. T may be an array, a table, or a string.
%   If given as a string, it must contain newlines processed by
%   sprintf or fprintf to transform the "\n" marker to the actual newline
%   character. If given as a table or array, it is passed to mat2latex to
%   generate the required string. If given as a table, the table
%   VariableNames are used as the column headers and the RowNames as the
%   row headers.
%
%   This function has a large number of parameter arguments that alter its
%   function:
%
%       'file' - allows the string to be directly output to a file.
%       Normally it just writes the string alone to that file, however
%       'insert' changes this behavior (see below). By default this will
%       ask before overwriting an existing file.
%
%       'insert' - instead of writing to a new file, this will instead try
%       to insert the table string into an existing file. It looks for
%       lines '%BEGIN <mark>' and '%END <mark>' and puts the table in
%       between those lines. <mark> is given by the 'marker' parameter (see
%       below). The %BEGIN or %END must be at the beginning of the line and
%       followed by at least one space before the <mark>. Everything
%       currently between the %BEGIN and %END lines is lost. This will back
%       up the existing file by copying it to <file>.bckp - but be wary of
%       running this function multiple times, as each time it will
%       overwrite the previous backup.
%
%       'marker' - the identifying string in the Latex file that it places
%       the table between in 'insert' mode. That is, if the value for this
%       is 'TABLE1', the table is placed between the lines "%BEGIN TABLE1"
%       and "%END TABLE1".
%
%       'overwrite' - boolean, if true, immediately overwrites the output
%       file without asking IF not operating in "insert" mode.
%
%       'caption' - the string to be placed in the \caption{} command for
%       the table environment. If not given, or if an empty string, the
%       \caption{} command will not be inserted at all.
%
%       'label' - the string to be place in the \label{} command for the
%       table environment. If not given, or if an empty string, the
%       \caption{} command will not be inserted at all.
%
%       'rownames' - cell array of strings that will be used as the row
%       names; overrides those given in T if it is a table.
%
%       'colnames' - cell array of strings that will be used as the column
%       names; overrides the VariableNames given in T if it is a table.
%
%       'm2l' - a cell array of options that mat2latex understands. These
%       will be passed to mat2latex if mat2latex is called to convert T to
%       a string.
%
%       'lines' - a 1-by-3 cell array of strings that have Latex commands
%       for the top, middle, and bottom horizontal lines. Defaults to
%       '\hline' for all three.

E = JLLErrors;

p = inputParser;
p.addParameter('file', '', @ischar);
p.addParameter('caption', '');
p.addParameter('label','');
p.addParameter('rownames', {});
p.addParameter('colnames', {});
p.addParameter('m2l', {});
p.addParameter('insert',false);
p.addParameter('marker','');
p.addParameter('overwrite',false);
p.addParameter('lines', {'\hline', '\hline', '\hline'});

p.parse(varargin{:});
pout = p.Results;

file_out = pout.file;
caption = pout.caption;
label = pout.label;
rownames = pout.rownames;
colnames = pout.colnames;
m2l_opts = pout.m2l;
do_insert = pout.insert;
insert_mark = pout.marker;
overwrite = pout.overwrite;
hlines = pout.lines;

if istable(T) && isempty(rownames)
    rownames = T.Properties.RowNames;
end
if istable(T) && isempty(colnames)
    colnames = T.Properties.VariableNames;
end

% Handle any % signs in the row/column names so that they aren't lost in
% the various sprintf calls.
rownames = strrep(rownames, '%', '%%');
colnames = strrep(colnames, '%', '%%');

% Check the input variables are of the right type
if ~ischar(T) && ~isnumeric(T) && ~istable(T)
    E.badinput('T must be a string, array, or table.')
end

if ~ischar(file_out)
    E.badinput('The parameter "file_out" must be a string')
end

if ~ischar(caption)
    E.badinput('The parameter "caption" must be a string')
end

if ~ischar(label)
    E.badinput('The parameter "label" must be a string')
end

if ~iscellstr(rownames)
    E.badinput('The parameter "rownames" must be a cell array of strings')
end

if ~iscellstr(colnames)
    E.badinput('The parameter "colnames" must be a cell array of strings')
end

if ~iscell(m2l_opts)
    E.badinput('The parameter "m2l" must be a cell array')
end

if ~isscalar(do_insert) || (~isnumeric(do_insert) && ~islogical(do_insert))
    E.badinput('The parameter "insert" must be a scalar number or boolean value')
end

if ~ischar(insert_mark)
    E.badinput('The parameter "marker" must be a string')
end

if ~isscalar(overwrite) || (~isnumeric(overwrite) && ~islogical(overwrite))
    E.badinput('The parameter "overwrite" must be a scalar number or boolean value')
end

if ~iscellstr(hlines) || numel(hlines) ~= 3
    E.badinput('The parameter "lines" must be a 3 element cell array of strings')
end

% Finally check the interrelationships
if do_insert && isempty(insert_mark)
    E.badinput('To use ''insert'' == true, you must specify a marker string')
end


%%%%% MAIN FUNCTION %%%%%
if ischar(T)
    latex_table_body = T;
elseif isnumeric(T) || istable(T)
    if istable(T)
        T = table2array(T);
    end
    latex_table_body = mat2latex(T, m2l_opts{:});
end

latex_table_body = strsplit(latex_table_body, '\n');
latex_table_body = latex_table_body(~iscellcontents(latex_table_body,'isempty'));

% If no row names, we don't need to reserve a column for them in the
% header, nor print them (obviously)
is_rownames = ~isempty(rownames);
n_columns = numel(strfind(latex_table_body{1},'&')) + 1 + is_rownames;

% Allow the column names to include one for the row names or not.
if ~isempty(colnames)
    if ~is_rownames
        % If no row names, the number of column names must be the number of
        % columns
        if numel(colnames) ~= n_columns
            E.badinput('Wrong number of column names: The parameter ''colnames'' must include an entry for each data column, if given')
        end
        
        header = [strjoin(colnames, ' & '), ' \\'];
    else
        % If row names, there need not be a column name for that one
        if numel(colnames) == n_columns - 1
            header = ['& ', strjoin(colnames, ' & '), ' \\'];
        elseif numel(colnames) == n_columns
            header = [strjoin(colnames, ' & '), ' \\ ', hlines{2}];
        else
            E.badinput('Wrong number of column names: The parameter ''colnames'' must include an entry for each data column, and may include one for the column of row names, if given');
        end
    end
    tabular_body = sprintf('%s %s\n ', hlines{1}, header);
else
    tabular_body = '';
end

for a=1:numel(latex_table_body)
    if a == numel(latex_table_body)
        sep = sprintf('%s \\n ', hlines{3});
    else
        sep = '\n ';
    end
    
    if is_rownames
        tabular_body = [tabular_body, rownames{a}, ' & ', latex_table_body{a}, sep];
    else
        tabular_body = [tabular_body, latex_table_body{a}, sep];
    end
end
tabular_body = sprintf(tex_in_printf(tabular_body));

% Now figure out how tabular should be formatted
if is_rownames
    fmt_str = sprintf('l%s', repmat('c', 1, n_columns-1));
else
    fmt_str = sprintf('%s', repmat('c', 1, n_columns));
end

% And whether caption and label should be included
if ~isempty(caption)
    cap_str = sprintf('\\caption{%s}\n ', caption);
else
    cap_str = '';
end

if ~isempty(label)
    label_str = sprintf('\\label{%s}\n ', label);
else
    label_str = '';
end

% Put it all together
table_str = [' \begin{table}\n',...
             ' \begin{tabular}{%1$s}\n'...
             ' %2$s\n'...
             ' \end{tabular}\n'...
             ' %3$s%4$s'...
             '\end{table}\n'];
table_str = tex_in_printf(table_str);
    
if isempty(file_out)
    varargout{1} = sprintf(table_str, fmt_str, tabular_body, cap_str, label_str);
elseif ~do_insert    
    if exist(file_out, 'file') && ~overwrite
        user_ans = ask_yn(sprintf('File %s exists. Overwrite?', file_out));
        if ~user_ans
            return
        end
    end
    fid = fopen(file_out,'w');
    if fid < 0
        E.callError('io_error', 'Could not open %s for writing', file_out);
    end
    fprintf(fid, table_str, fmt_str, tabular_body, cap_str, label_str);
    fclose(fid);
else
    insert_into_file(sprintf(table_str, fmt_str, tabular_body, cap_str, label_str), file_out, insert_mark);
end

end

function insert_into_file(table_str, file_name, mark)
E = JLLErrors;
table_str = strrep(table_str, '%', '%%');

bckp_name = sprintf('%s.bckp', file_name);
out_name = sprintf('%s.matlab', file_name);

if ~exist(file_name, 'file')
    E.filenotfound(file_name);
elseif exist(file_name, 'dir')
    E.callError('file_is_dir', 'Given file (%s) is actually a directory', msg);
end

[stat, msg] = copyfile(file_name, bckp_name);
if ~stat
    E.callError('could_not_copy','Could not make backup copy: %s', msg);
end

[curr_fid, msg] = fopen(file_name, 'r');
if curr_fid < 0
    E.callError('could_not_open','Could not open %s for reading: %s', file_name, msg);
end

[new_fid, msg] = fopen(out_name, 'w');
if new_fid < 0
    E.callError('could_not_open','Could not open %s for writing: %s', file_out, msg);
end

found_mark = false;
in_mark = false;
mark_regex_start = sprintf('^%%BEGIN\\s+%s',mark);
mark_regex_end = sprintf('^%%END\\s+%s',mark);
mark_regex = mark_regex_start;
tline = fgets(curr_fid);
while ischar(tline)
    nline = strrep(tline, '%', '%%');
    if ~in_mark
        fprintf(new_fid, tex_in_printf(nline));
    end
    if ~isempty(regexp(tline, mark_regex, 'once'))
        if ~in_mark
            if found_mark
                warning('Mark %s multiply defined', mark);
            end
            found_mark = true;
            in_mark = true;
            mark_regex = mark_regex_end;
        elseif in_mark
            fprintf(new_fid, tex_in_printf(table_str));
            fprintf(new_fid, tex_in_printf(nline));
            in_mark = false;
            mark_regex = mark_regex_start;
        end 
    end
    
    tline = fgets(curr_fid);
end
fclose(curr_fid);
fclose(new_fid);

if in_mark
    E.callError('bad_mark', 'Could not find ending mark line for %s', mark);
elseif ~found_mark
    E.callError('no_mark', 'Could not find mark %s', mark);
end

movefile(out_name, file_name);

end
