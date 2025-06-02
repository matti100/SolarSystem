function [Ephemeris, x0, x_end, Trajectory, Velocity] = extractEphemeris(startTime, stopTime, unit, dt, frame, onlyIC_flag)

% Measurement units definition
if unit == 1
    % 1: day        Km/day
    um = 'KM-D';
    step_size = strcat( num2str(dt), 'd' );
elseif unit == 0
    % 0: second     Km/s
    um = 'KM-S';
    step_size = strcat( num2str(dt), 's' );
else
    error('unit flag not valid');
end

% Reference frame definition
if frame == 1
    % 1: Solar system baricenter
    refFrame = '500@0';
elseif frame == 0
    % 0: Sun
    refFrame = '500@10';
else
    error('frame flag not valid');
end

% Initialize variables
x0 = zeros(54,1);
x_end = zeros(54,1);

% Define API-url
url = 'https://ssd.jpl.nasa.gov/api/horizons.api';

% Define options
options = weboptions('Timeout', 30); % set a 30sec Timeout
    
% Defining planets IDs  
bodyCode = struct(...
    'Sun', '10', ...   
    'Mercury', '199', ...  
    'Venus', '299', ...      
    'Earth', '399', ...  
    'Mars', '499', ...    
    'Jupiter', '599', ...      
    'Saturn', '699', ...   
    'Uranus', '799', ...           
    'Neptune', '899');

dataLength = 21; % 21 characters used for each quantity

% Extract fieldnames from the planetsCode struct
fieldNames = fieldnames(bodyCode);

disp("--- Extracting Ephemeris ---");
fprintf("Start Date: %s", startTime);
fprintf("\nStop Date: %s", stopTime);
fprintf("\nURL: %s\n\n", url);

% If you want to extract only the IC, then lower the stopTime
if onlyIC_flag
    % Set the stopTime as the day after the startTime
    stopTime = strcat(startTime(1:end-1), num2str( str2num(startTime(end)) + 1) );
    fprintf("Ephem Stop Date: %s\n", stopTime);
end

for i = 1:length(fieldNames)
    
    % # planets
    j = 54*(i-1)/9;

    body = fieldNames{i};
    code = bodyCode.(body);
    
    fprintf("Body: %s", body);

    % Get response
    response = webread(...
        url, ...
        'format', 'text', ...
        'COMMAND', code, ...
        'EPHEM_TYPE', 'VECTORS', ...
        'CENTER', refFrame, ...
        'REF_PLANE', 'E', ...
        'VEC_TABLE', '2', ...
        'OUT_UNITS', um, ...
        'CAL_FORMAT', 'CAL', ...
        'START_TIME', startTime, ...
        'STOP_TIME', stopTime, ...
        'STEP_SIZE', step_size, ...
        options...
        );


    % Check if ephemeris are available for that time interval
    if strcmp(response(end-74: end-74+11), "No ephemeris")
        fprintf("\n\n !! Warning : %s\n\n", response(end-74:end));
        
        if onlyIC_flag == 0
            error('No ephemeris found');
        end
    end
            
    Ephemeris.(body) = response;

    %% Extract useful data
    ref1 = '$$SOE';
    ref2 = '$$EOE';
    index1 = strfind(response, ref1);
    index2 = strfind(response, ref2);

    data = response(index1:index2);
    
    refX = 'X';
    indexX = strfind(data, refX);
    refY = 'Y';
    indexY = strfind(data, refY);
    refZ = 'Z';
    indexZ = strfind(data, refZ);
    refVx = 'VX';
    indexVx = strfind(data, refVx);
    refVy = 'VY';
    indexVy = strfind(data, refVy);
    refVz = 'VZ';
    indexVz = strfind(data, refVz);

    stop_flag = true;
    n = 1;      % counter for indexX,Y,Z
    m = 1;      % counter fo indexVx,Vy,Vz

    while stop_flag
        
        % Extract positions and velocities of all planets
        i_x = indexX(n);
        i_y = indexY(n);
        i_z = indexZ(n);
        i_vx = indexVx(m);
        i_vy = indexVy(m);
        i_vz = indexVz(m);

        X = str2num( data(i_x+3:i_x+3+dataLength) );
        Y = str2num( data(i_y+3:i_y+3+dataLength) );
        Z = str2num( data(i_z+3:i_z+3+dataLength) );
        Vx = str2num( data(i_vx+3:i_vx+3+dataLength) );
        Vy = str2num( data(i_vy+3:i_vy+3+dataLength) );
        Vz = str2num( data(i_vz+3:i_vz+3+dataLength) );
    
        % X = str2num( data(i_x+3:i_y-2) );
        % Y = str2num( data(i_y+3:i_z-2) );
        % Z = str2num( data(i_z+3:i_vx-3) );
        % Vx = str2num( data(i_vx+3:i_vy-2) );
        % Vy = str2num( data(i_vy+3:i_vz-2) );
        % Vz = str2num( data(i_vz+3:i_vz+dataLength+3) );

        if (n == 1)
            Trajectory.(body) = [X; Y; Z];
            Velocity.(body) = [Vx; Vy; Vz];
        else
            Trajectory.(body) = horzcat(Trajectory.(body), [X; Y; Z]);
            Velocity.(body) = horzcat(Velocity.(body), [Vx; Vy; Vz]);
        end

        
        if (n == 1)
            % Extract initial position and velocity
            x0(j+1) = X;
            x0(j+2) = Y;
            x0(j+3) = Z;
            x0(j+4) = Vx;
            x0(j+5) = Vy;
            x0(j+6) = Vz;

        elseif (n == indexX(end))
            % Extract final position and velocity
            x_end(j+1) = X;
            x_end(j+2) = Y;
            x_end(j+3) = Z;
            x_end(j+4) = Vx;
            x_end(j+5) = Vy;
            x_end(j+6) = Vz;
        end
        
        % Update counters
        n = n + 2;               % ( X, Y, Z are double then Vx, Vy, Vz )
        m = m + 1;

        % Stop logic
        if (n == length(indexX)-1)
            stop_flag = false;
        end
        
        if onlyIC_flag == 1
            stop_flag = false;
        end
        
    end

    fprintf("\t...complete\n");

end

disp("--- Extraction successful ---");

end

