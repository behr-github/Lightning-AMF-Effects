classdef JLLErrors<handle
    %JLLErrors A class allowing more consistent custom error handling
    %   This class has a number of methods that can be called from within
    %   an error statement to return an appropriate, customized error
    %   structure.  These methods will construct the error identifier of
    %   the form 'calling_function_name:tag', which is the general form I
    %   have been using in my custom errors.
    %
    %   There are several generic errors that allow you to throw an
    %   exception with minimal customization; as you find exceptions that
    %   are common to many script, feel free to add them as methods to this
    %   class.
    %
    %   Three special methods are the addCustomError, callCustomError, and
    %   callError methods.  callError allows you to construct an error with
    %   the function name in the ID and everything else set as the
    %   arguments.  addCustomError and callCustomError work together to let
    %   you define reusable errors that will be available within that
    %   instance of this class.
    
    properties(SetAccess=private)
        callingfxn % The function that creates this class, set by constructor
    end
    
    properties(Access=private, Hidden=true)
        % The error stack needs to be set
        %errstack
        
        % Default tags and error messages for a number of common error
        % types
        scalar_tag = 'not_scalar';
        scalar_msg = 'The variable ''%s'' must be a scalar value';
        fnf_tag = 'file_not_found';
        fnf_msg = 'Could not load the file: %s';
        invalidinput_tag = 'invalid_input';
        invalidinput_msg = 'An input to this function is not valid; see documentation';
        invalidvar_tag = 'invalid_var';
        invalidvar_msg = 'The variable ''%s'' is not valid, see documentation';
        invalidvartype_tag = 'bad_var_type';
        invalidvartype_msg = 'The variable ''%s'' is a %s and must be an %s';
        numArgs_tag = 'wrong_number_arguments';
        numArgs_msg = 'The number of arguments must be between %d and %d.';
        numelMismatch_tag = 'numel_mismatch';
        numelMismatch_msg = 'The variables %s must have the same number of elements';
        dimMismatch_tag = 'dim_mismatch';
        dimMismatch_msg = 'The variables %s must have the same dimensions';
        badgeo_tag = 'latlon_data_mismatch';
        badgeo_msg = 'The data matrix input is not consistent with the size of the lat/lon inputs';
        tmf_tag = 'non_unique_file';
        tmf_msg = 'More than one %s file found meeting given criteria - try a more specific file prefix, if applicable';
        usercancel_tag = 'user_cancel';
        usercancel_msg = 'User cancelled run.';
        
        % A list of custom identifiers and messages and reference to those
        % entries
        custom_ids = {};
        custom_msgs = {};
        custom_refs = {};
    end
    
    methods
        function obj = JLLErrors()
            %Will identify the calling function, or default to 'base' if the calling function cannot be identified (i.e. if this is created from the command line).  This name will be used as the first part of the error identifier.
            s = dbstack(1);
            if numel(s)>0
                obj.callingfxn = s(1).name;
            else
                obj.callingfxn = 'base';
            end
        end
        
        function errstruct = notScalar(obj, varname)
            % Takes a single variable name.  Returns an error message that said variable must be a scalar.
            warning('JLLErrors.notScalar is deprecated - use "badvartype" instead');
            msg = sprintf(obj.scalar_msg,varname);
            errstruct = obj.makeErrStruct(obj.scalar_tag,msg);
        end
        
        function errstruct = badinput(obj,varargin)
            %Returns an error structure for a bad input.  With no arguments, the default message will be displayed; if one argument is passed, that specifies the message.
            if numel(varargin) == 0
                msg = obj.invalidinput_msg;
            elseif numel(varargin)==1
                msg = varargin{1};
            end
            
            errstruct = obj.makeErrStruct(obj.invalidinput_tag, msg);
        end
        
        function errstruct = badvar(obj, varname)
            % Takes a variable name as its only argument, will produce an error with a message indicating that that variable is not valid with no additional information.
            msg = sprintf(obj.invalidvar_msg,varname);
            
            errstruct = obj.makeErrStruct(obj.invalidvar_tag,msg);
        end
        
        function errstruct = badvartype(obj, var, rightclass)
            % Takes a variable of any type and a string with the correct class that that variable should be.  Will output a message indicating that the variable var is of type class(var) and should be of type rightclass.
            varname = inputname(2);
            wrongclass = class(var);
            msg = sprintf(obj.invalidvartype_msg, varname, wrongclass, rightclass);
            errstruct = obj.makeErrStruct(obj.invalidvartype_tag, msg);
        end
        
        function errstruct = badgeo(obj)
            errstruct = obj.makeErrStruct(obj.badgeo_tag, obj.badgeo_msg);
        end
        
        function errstruct = filenotfound(obj, filename)
            % Error when a file could not be found to be loaded. Takes one argument which describes the file that couldn't be loaded.
            msg = sprintf(obj.fnf_msg,filename);
            errstruct = obj.makeErrStruct(obj.fnf_tag, msg);
        end
        
        function errstruct = toomanyfiles(obj, varargin)
            % Error for use when finding file names to load using wildcard characters.  Takes one or no arguments, if one is given, it will describe what kind of file is trying to be loaded.
            if numel(varargin)>0
                msg = sprintf(obj.tmf_msg,varargin{1});
            else
                msg = sprintf(obj.tmf_msg,'');
            end
            
            errstruct = obj.makeErrStruct(obj.tmf_tag, msg);
        end
        
        function errstruct = numelMismatch(obj, varargin)
            % Error for when an arbitrary number of variables do not have the same number of elements. Takes at least two arguments (variable names as strings) up to an unlimited number of arguments.
            if numel(varargin)<2; error('JLLErrors:numelMismatch:too_few_inputs','JLLErrors.numelMismatch needs at least 2 variable names'); end
            varnamespec = [repmat('''%s'', ',1,numel(varargin)-1),'''%s'''];
            msgspec = sprintf(obj.numelMismatch_msg, varnamespec);
            msg = sprintf(msgspec, varargin{:});
            errstruct = obj.makeErrStruct(obj.numelMismatch_tag, msg);
        end
        
        function errstruct = dimMismatch(obj, varargin)
            % Error for when an arbitrary number of variables do not have the same number of elements. Takes at least two arguments (variable names as strings) up to an unlimited number of arguments.
            if numel(varargin)<2; error('JLLErrors:dimMismatch:too_few_inputs','JLLErrors.numelMismatch needs at least 2 variable names'); end
            varnamespec = [repmat('''%s'', ',1,numel(varargin)-1),'''%s'''];
            msgspec = sprintf(obj.dimMismatch_msg, varnamespec);
            msg = sprintf(msgspec, varargin{:});
            errstruct = obj.makeErrStruct(obj.dimMismatch_tag, msg);
        end
        
        function errstruct = userCancel(obj)
            % Error for when the user cancels out of a dialogue box. Takes no arguments.
            errstruct = obj.makeErrStruct(obj.usercancel_tag,obj.usercancel_msg);
        end
        
        function errstruct = numberArguments(obj,nmin,nmax)
            % Error for inputting the wrong number of arguments to a function. 
            narginchk(3,3);
            msg = sprintf(obj.numArgs_msg,nmin,nmax);
            errstruct = obj.makeErrStruct(obj.numArgs_tag, msg);
        end
        
        function errstruct = callError(obj, tag, msg)
            % A very simple method that creates an error with a custom message and id tag (second and first arguments respectively). The resulting error will have the identifier 'callingfxn:tag' and the specified message.
            errstruct = obj.makeErrStruct(tag,msg);
        end
        
        function obj = addCustomError(obj, varargin)
            % Allows the user to define custom error messages; this is
            % useful if there are errors you expect to use several times in
            % a particular function or script that are not common enough to
            % code into this class.  Takes 2 or 3 arguments: in the two arg
            % form, you pass the second half of the error identifer (i.e.
            % the ID will become 'calling_fxn_name:your_identifier'), and
            % the message to use.  In the three arg form, you pass a
            % reference string as the first argument; this is the string to
            % use when calling the error using the callError method.  If
            % only two arguments are given, the reference is set by default
            % to be the same as the user defined part of the identifier.
            if numel(varargin) == 2
                newid = varargin{1};
                newref = varargin{1};
                newmsg = varargin{2};
            elseif numel(varargin) == 3
                newref = varargin{1};
                newid = varargin{2};
                newmsg = varargin{3};
            end
            obj.custom_ids{end+1} = newid;
            obj.custom_msgs{end+1} = newmsg;
            obj.custom_refs{end+1} = newref;
        end
        
        function errstruct = callCustomError(obj,ref)
            % Call this function with the reference string set using
            % addCustomError to return a user defined error.
            xx = strcmp(obj.custom_refs,ref);
            if sum(xx)==0; error('JLLErrors:error_reference','Error reference was invalid or not set'); end
            tag = obj.custom_ids{xx};
            msg = obj.custom_msgs{xx};
            
            errstruct = obj.makeErrStruct(tag,msg);
        end
    end
    
    methods(Access=private, Hidden=true)
        function errstruct = makeErrStruct(obj, tag, msg)
            errid = sprintf('%s:%s',obj.callingfxn,tag);
            stack = dbstack(2,'-completenames');
            errstruct = struct('identifier',errid,'message',msg,'stack',stack);
        end
    end
    
end

