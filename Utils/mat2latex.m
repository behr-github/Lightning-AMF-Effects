function mat2latex(M,format_spec)
% MAT2LATEX Prints a matrix, M in Latex-table format
%   Takes a numeric matrix or a cell array  and prints out to the command window that matrix
%   formatted such that it can be copied into a Latex table. Optional second
%   argument is a format spec understood by fprintf to control how the numbers
%   will be formatted. Defaults to %g.

if nargin < 2
    format_spec = '%g';
end

E = JLLErrors;
if ~iscell(M) && ~isnumeric(M)
    E.badinput('Expecting a numeric matrix or a cell array')
elseif ndims(M) > 2
    E.badinput('Higher dimension arrays don''t make sense as a table')
end

if ~iscell(M)
    M = mat2cell(M, ones(1, size(M,1)), ones(1, size(M, 2)));
end

% The regular print string
fstr1 = sprintf('%s & ', format_spec);
fstr1b = '%s & ';
% The end-of-line print string
fstr2 = sprintf(tex_in_printf('%s \ \n',3), format_spec);
fstr2b = tex_in_printf('%s \ \n',3);


sz = size(M);
for a=1:sz(1)
    for b=1:sz(2)
        if isnumeric(M{a,b})
            if b < sz(2)
                fstr=sprintf(fstr1, M{a,b});
            else
                fstr=sprintf(fstr2, M{a,b});
            end
            % Replace e notation with x10^, remove + sign and leading 0s
            fstr = regexprep(fstr,'e',' \\times 10^{');
            fstr = regexprep(fstr,'(?<={-?)+?0*(?=\d+)','');
            %fstr = regexprep(fstr,'(?<={)\+?0*(?=\d+)','');
            fstr = regexprep(fstr,'(?<={-?\d+) (?=&|\\)','}$ ');
            if ~isempty(strfind(fstr,'{')) && isempty(strfind(fstr,'}'))
                fstr = strcat(fstr,'}$'); % handles the final line
            end
            % Numbers without an exponent will not have the ending $ added
            % yet
            if isempty(strfind(fstr,'$'))
                fstr = regexprep(fstr,' (?=&|\\)','$ ');
            end
            fprintf('$%s',fstr);
        else
            if b < sz(2)
                fprintf(fstr1b, M{a,b});
            else
                fprintf(fstr2b, M{a,b});
            end
        end
    end
end

fprintf('\n')
