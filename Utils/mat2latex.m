function mat2latex(M, format_spec, uncertainty_dim)
% MAT2LATEX( M ) Prints a matrix, M in Latex-table format
%   Takes a numeric matrix or a cell array  and prints out to the command window that matrix
%   formatted such that it can be copied into a Latex table. 
%
%   MAT2LATEX( M, FORMAT_SPEC ) Optional second argument is a format spec
%   string understood by fprintf to control how the numbers will be
%   formatted. Defaults to %g.
%
%   MAT2LATEX( M, UNCERTAINTY_DIM ) Allows you to include uncertainty
%   values that will be printed as v \pm u (x 10^e). Uncertainties must
%   alternate with values along the dimension given. For example, if this
%   is 1, then the first row of M should contain values, the second row the
%   uncertainties in the first row, the third row the next set of values,
%   the fourth row the uncertainties for the third row, and so on.
%
%   MAT2LATEX( M, 'uncertainty', UNCERTAINTY_DIM ) will format any numbers
%   such that the last digit in the value is the first digit in the
%   uncertainty.
%
%   MAT2LATEX( M, 'u', UNCERTAINTY_DIM ) is a shorthand for the last
%   syntax.
%
%   MAT2LATEX( M, FORMAT_SPEC, UNCERTAINTY_DIM ) combines the previous two
%   syntaxes.

if nargin < 2
    format_spec = '%g';
    uncertainty_dim = 0; % do not include uncertainty.
elseif nargin == 2
    if isnumeric(format_spec)
        uncertainty_dim = format_spec;
        format_spec = '%g';
    else
        uncertainty_dim = 0;
    end
end


E = JLLErrors;
if ~iscell(M) && ~isnumeric(M)
    E.badinput('Expecting a numeric matrix or a cell array')
elseif ~ismatrix(M)
    E.badinput('Higher dimension arrays don''t make sense as a table')
end

if ~ischar(format_spec)
    E.badinput('FORMAT_SPEC must be a string')
end
if ~isnumeric(uncertainty_dim) || uncertainty_dim < 0 || uncertainty_dim > 2
    E.badinput('UNCERTAINTY_DIM must be 0 (no uncertainties given), 1, or 2')
elseif uncertainty_dim > 0 && mod(size(M,uncertainty_dim),2) ~= 0
    E.badinput('M must have an even number of entries along the UNCERTAINTY_DIM')
end

if ~iscell(M)
    M = num2cell(M);
end

warn_neg_uncert = false;

% The regular print string
fstr1 = sprintf('%s & ', format_spec);
fstr1b = '%s & ';
% The end-of-line print string
fstr2 = sprintf(tex_in_printf('%s \ \n',3), format_spec);
fstr2b = tex_in_printf('%s \ \n',3);



sz = size(M);
if uncertainty_dim == 1
    astep = 2;
else
    astep = 1;
end

if uncertainty_dim == 2
    bstep = 2;
else
    bstep = 1;
end
for a=1:astep:sz(1)
    for b=1:bstep:sz(2)
        if isnumeric(M{a,b})
            if strcmpi(format_spec,'uncertainty') || strcmpi(format_spec, 'u')
                this_fstr = format_uncert;
            else
                this_fstr = format_normal;
            end
            
            fprintf('$%s',this_fstr);
        else
            if b < (sz(2) - bstep+1)
                fprintf(fstr1b, M{a,b});
            else
                fprintf(fstr2b, M{a,b});
            end
        end
    end
end

fprintf('\n')

if warn_neg_uncert
    warning('Negative values of uncertainty have been replaced with positive ones');
end

% Nested functions to parse numbers
    function fstr = format_uncert
        % Get the value and uncertainty following the usual rules, but
        % round them to the first non-zero place in the uncertainty first.
        v = M{a,b};
        if uncertainty_dim == 1
            u = M{a+1,b};
        elseif uncertainty_dim == 2
            u = M{a,b+1};
        else
            E.badinput('Cannot format the values by rounding to the first place of the uncertainty if no uncertainty given (uncertainty_dim must be >0).');
        end
        if u < 0
            warn_neg_uncert = true;
        end
        u = abs(u);
        place = 10^(floor(log10(u)));
        v = round(v/place)*place;
        u = round(u/place)*place;
        
        % Calculate the number of significant figures to be left in v
        % Testing if log10 is infinite rather than just u == 0 or v == 0
        % because it is if log10 = -Inf that directly causes the problem,
        % so I want to avoid any little floating point errors
        if isinf(log10(v)) && isinf(log10(u))
            fstr = '0 \pm 0$';
            if b < (sz(2) - bstep+1)
                fstr=sprintf(fstr1b, fstr);
            else
                fstr=sprintf(fstr2b, fstr);
            end
        else
            if ~isinf(log10(v))
                nsig = floor(log10(v)) - floor(log10(u)) + 1;
                
                % Format v accordingly and insert uncertainty
                fstr = sprintf('%#.*g',nsig,v);
            elseif ~isinf(log10(u));
                % If the value has become 0, then we will need to use u to set
                % the format, as long as it isn't 0!
                ustr = sprintf('%.1g',u);
                fstr = regexprep(ustr,'\d','0');
            else
                %fstr = '0 \pm 0
            end
            if b < (sz(2) - bstep+1)
                fstr=sprintf(fstr1b, fstr);
            else
                fstr=sprintf(fstr2b, fstr);
            end
            fstr = format_exponent(fstr);
            fstr = insert_uncertainty(u, fstr);
            fstr = strrep(fstr,'\\\\','\\');
        end
    end

    function fstr = format_normal
        if b < (sz(2) - bstep+1)
            fstr=sprintf(fstr1, M{a,b});
        else
            fstr=sprintf(fstr2, M{a,b});
        end
        
        fstr = format_exponent(fstr);
        
        % Add in uncertainty, if desired
        if uncertainty_dim > 0
            if uncertainty_dim == 1
                u = M{a+1,b};
            elseif uncertainty_dim == 2
                u = M{a,b+1};
            end
            if u < 0
                warn_neg_uncert = true;
                u = abs(u);
            end
            fstr = insert_uncertainty(u, fstr);
        end
    end

end

function fstr = format_exponent(fstr)
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
end

function fstr = insert_uncertainty(u, fstr)
% convert u to the same power of 10 as the value. first
% find the exponent in the string and convert it to a
% number.
[es, ee] = regexp(fstr,'(?<=10\^{)-?\d*(?=})');
if isempty(es)
    val_exp = 0;
else
    val_exp = str2double(fstr(es:ee));
end
u = u * 10^(-val_exp);
% figure out how many figures after the decimal point there
% are in the value
[ds,de] = regexp(fstr,'\.\d*');
if ~isempty(ds)
    ndec = de - ds;
else
    ndec = 0;
    de = regexp(fstr, '(\$|\o{40}\\times)')-1;
end
ustr = sprintf('%.*f', ndec, u);
% insert the uncertainty
fstr = sprintf('%s \\pm %s%s',fstr(1:de),ustr,fstr(de+1:end));
end

function [str, dec_ind] = round_str(str, index)
% Rounds the number in the string to the given index
dec_ind = regexp(str,'\.');
if isempty(dec_ind)
    dec_ind = length(str)+1;
elseif numel(dec_ind) > 1
    error('Multiple decimal points found');
end

% Handle rounding. Round numerically to the requested precision, then we
% will reformat the string afterwards.
num = str2double(str);
if index < dec_ind
    place = 10^(index - dec_in + 1);
else
    place = 10^(index - dec_in);
end
num = round(num * place)/place;

% Now either make the following numbers 0s or remove them, as appropriate
if index < dec_ind
    str(index+1:dec_ind-1)='0';
    str(dec_ind:end)=[];
else
    str(index+1:end)=[];
end
end
