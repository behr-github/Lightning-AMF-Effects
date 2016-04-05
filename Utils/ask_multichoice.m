function [ user_ans ] = ask_multichoice( prompt, allowed_options, varargin )
%ASK_MULTICHOICE Will ask the user to select from the allowed options
%   Asks the user the question given in PROMPT to choose from the cell
%   array ALLOWED_OPTIONS. Two parameters exist:
%
%   'default' - choose a default value for the response.
%
%   'softquit' - boolean, if false (default), will exit if the user enters
%   'q' at any point by throwing an error. If true, this will cause the
%   function to return a 0 (the number not the string).
%
%   This function always returns a string (unless softquit is true) so
%   you'll need to parse the string back into a number if you want it so.
%   Answers are also always returned in lower case.
%
%   Josh Laughner <joshlaugh5@gmail.com> 26 Jan 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
if ~ischar(prompt)
    E.badinput('prompt must be a string')
end
if ~iscellstr(allowed_options)
    E.badinput('allowed_options must be a cell array of strings')
end

p = inputParser;
p.addParameter('default',0);
p.addParameter('softquit',false,@(x) (islogical(x) && isscalar(x)));
p.parse(varargin{:});
pout = p.Results;

default = pout.default;
if default == 0
    use_default = false;
else
    if ~ischar(default) || ~ismember(default, allowed_options)
        E.badinput('default must be a string that is one of those in allowed_options')
    end
    use_default = true;
end
softquit = pout.softquit;

if use_default
    fprintf('%s (%s - %s is default): ', prompt, strjoin(allowed_options, ', '), default);
else
    fprintf('%s (%s): ', prompt, strjoin(allowed_options, ', '));
end

while true
    user_ans = lower(input('', 's'));
    if use_default && isempty(user_ans)
        user_ans = default;
        return
    elseif ismember(user_ans, allowed_options)
        return
    elseif strcmpi(user_ans, 'q')
        if softquit
            user_ans = 0;
        else
            E.userCancel()
        end
    else
        fprintf('You must choose one of %s, or enter q to quit: ', strjoin(allowed_options, ', '));
    end
    
end

end

