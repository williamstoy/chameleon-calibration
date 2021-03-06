classdef Chameleon
    properties
       verbose
       ser
       current_wavelength
    end
    
    methods
        %########################################################################%
        %# INITIALIZATION FUNCTIONS
        %########################################################################%
        function obj = Chameleon(com_port, verbose)
            obj.verbose = verbose;
            obj.ser = serialport(com_port, 19200, 'Timeout', 1);  % open serial port
            configureTerminator(obj.ser, 'CR');
            % initialize the connection by writing a carriage return
            obj.sendCmdGetResponse(' ');
            pause(0.05);

            obj.current_wavelength = 1000;
        end


        %########################################################################%
        %# COMMAND / RESPONSE PRIMITIVE FUNCTIONS
        %########################################################################%	
        function response = getResponse(obj, cmd)
            response = readline(obj.ser);
            %format the command so that it can be taken out with a regular expression
            cmd = regexprep(cmd, '[\?]', '\?');
            cmd = regexprep(cmd, '[\=]', '\=');
            %strip out \r and \n, CHAMELEON> and the command
            response = regexprep(response, ['[\\r|\\n|(CHAMELEON\>)|(', cmd, ')|\s]'], '');
        end

        function response = sendCmdGetResponse(obj, cmd)
            writeline(obj.ser, cmd);
            pause(0.05);
            response = obj.getResponse(cmd);
        end


        %########################################################################%
        %# WAVELENGTH FUNCTIONS
        %########################################################################%
        function setWavelength(obj, n)
            obj.current_wavelength = n;
            if(obj.verbose)
                disp(['Tuning the laser to ', n, 'nm']);
            end
            obj.sendCmdGetResponse(['vw=', num2str(n)]);
        end    

        function response = queryTunedStatus(obj)
            response = obj.sendCmdGetResponse('?ts');
        end

        function setWavelengthBlocking(obj, n)
            obj.setWavelength(n);
            while (obj.queryTunedStatus() ~= '0')
                if(obj.verbose)
                    disp('Waiting for tuning to finish...');
                end
                pause(0.25);
            end
            if(obj.verbose)
                disp('Laser is Tuned');
            end
        end

        %########################################################################%
        %# SHUTTER FUNCTIONS
        %########################################################################%
        function openShutter(obj)
            if(obj.verbose)
                disp('Opening the Shutter');
            end
            obj.sendCmdGetResponse('s=1');
        end

        function closeShutter(obj)
            if(obj.verbose)
                disp('Closing the Shutter');
            end
            obj.sendCmdGetResponse('s=0');
        end

        function response = queryShutterStatus(obj)
            response = obj.sendCmdGetResponse('?s');
        end

        function openShutterBlocking(obj)
            obj.openShutter()
            while (obj.queryShutterStatus() ~= '1')
                if(obj.verbose)
                    disp('Waiting for the shutter to open..');
                end
                pause(0.25);
            end

            if(obj.verbose)
                disp('Shutter is OPEN');
            end
        end

        function closeShutterBlocking(obj)
            obj.closeShutter();
            while (obj.queryShutterStatus() ~= '0')
                if(obj.verbose)
                    disp('Waiting for the shutter to close..');
                end
                pause(0.25);
            end

            if(obj.verbose)
                disp('Shutter is CLOSED');
            end
        end

        function queryRelativeHumidity(obj)
            if(obj.verbose)
                disp('Getting Relative Humidity');
            end
            obj.sendCmdGetResponse('?rh');
        end

        %########################################################################%
        %# SERIAL FUNCTIONS
        %########################################################################%
        function delete(obj)
            if obj.verbose
                disp('Deleting Object and Closing Serial Port');
            end
            
            obj.closeShutterBlocking();
            
            clear obj.ser;
        end
    end
end