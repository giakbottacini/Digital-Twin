function DATA = CSM_read_DataFile(data_file, dim)
%CSM_READ_DATAFILE data_file parser
%
%   DATA = CSM_READ_DATAFILE(DATA_FILE) read the file specified by the string
%   DATA_FILE and put the fields values into the struct DATA

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>

if nargin < 2 || isempty(dim)
    dim = 2;
end

    %permette di definire una certa distribuzione dei parametri di rigidezza
    %all'interno della trave
spatial=1;

%% Set Default values

    %i flag verranno utilizzati per indicare quali siano i lati di contorno
    %in cui vengono definite le condizioni al contorno considerate. Questo deve
    %essere fatto per ogni direzione di interesse (dipenderà dalla dimensione
    %del problema)
for d = 1 : dim
    DATA.flag_dirichlet{d} = [];
    DATA.flag_neumann{d}   = [];
    DATA.flag_pressure{d}  = [];
    DATA.flag_robin{d}     = [];
    DATA.flag_clamp_points{d} = [];
end

DATA.flag_dirichletNormal = [];
        
switch dim
    
    case 2
        DATA.bcDir          = @(x,y,t,param)(0.*x.*y);
        DATA.bcNeu          = @(x,y,t,param)(0.*x.*y);
        DATA.bcPrex         = @(x,y,t,param)(0.*x.*y);
        DATA.force          = @(x,y,t,param)(0.*x.*y);
                            
    case 3
        
        DATA.bcDir          = @(x,y,z,t,param)(0.*x.*y);
        DATA.bcNeu          = @(x,y,z,t,param)(0.*x.*y);
        DATA.bcPrex         = @(x,y,z,t,param)(0.*x.*y);
        DATA.force          = @(x,y,z,t,param)(0.*x.*y);
end

%andamento temporale dei parametri meccanici
if spatial
    DATA.YoungSpatial = @(x,y,z,t,param)(0.*x.*y);
end

DATA.Output.ComputeVonMisesStress = false;
DATA.Output.ComputeAllStresses = false;

%% Read data_file and put problem-data into the DATA struct

%L: eval() permette di aprire un file .mat (in questo caso 'datafile.mat')
%e di eseguirlo mantenendo, nel workspace della funzione invocante, le
%variabili definite all'interno di esso
eval(data_file);
data_fields_name = fieldnames(data); %L: 'fieldnames' permette di riscrivere
             %in forma di vettore le stringhe contenenti le etichette utilizzate
             %all'interno della struttura 'data'

for i = 1 : length(data_fields_name)
    %L: i dati sono già stati raccolti nella data structure 'data'. Essi,
    %vengono però riscritti in una seconda struttura dati 'DATA'. Questo
    %perchè, nella successiva funzione dataParser ("analizzatore di dati"),
    %questi subiranno alcune modifiche.
    eval(['DATA.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
    
end


[ DATA ] = dataParser( DATA ); %L: dataParser opera una riscrittura dei dati
                    %(guidata dall'autore), in modo da tradurre alcuni dei
                    %comandi in ingresso in qualcosa di più facilmente
                    %interpretabile da MATLAB. Questo rende la gestione dei
                    %dati del problema più user-friendly.

end