%%  MATLAB class to capture MOSFET properties.
%   Date of creation:   09-04-2019
%   Last Modified:      09-04-2019

classdef mosfet

    %%  Parameters and data
    properties
        ciss;
        coss;
        crss;
        param;
    end
    
    %%  Methods
    methods
        
        %%  Constructor
        function obj = mosfet(ciss, coss, crss, param)
            obj.ciss = ciss;
            obj.coss = coss;
            obj.crss = crss;
            obj.param = param;
        end
        
    end

end
