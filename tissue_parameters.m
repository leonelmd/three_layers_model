%%% Tissue parameters to predict dielectric properties
%%% S Gabriel, R Lau, and C Gabriel, "The dielectric properties of
%%% biological tissues: III. Parametric models for the dielectric spectrum
%%% of tissues", Ohys Med Biol, vol. 41, pp. 2271-2293, 1996.
%%% (see Table 1, p2291)
%%%
%%% Leo Medina, Apr 2013

%Fat (not infiltrated)
fat = struct(...
    'e_inf',2.5,...    
    'de_1',3,...
    'tau_1',7.96e-12,...
    'alpha_1',0.2,...
    'de_2',15,...
    'tau_2',15.92e-9,...
    'alpha_2',0.1,...
    'de_3',3.3e4,...
    'tau_3',159.15e-6,...
    'alpha_3',.05,...
    'de_4',1e7,...
    'tau_4',7.958e-3,...
    'alpha_4',.01,...
    'sigma',.01 ...
);

%muscle (across)
muscle = struct(...
    'e_inf',4,...    
    'de_1',50,...
    'tau_1',7.23e-12,...
    'alpha_1',.1,...
    'de_2',7000,...
    'tau_2',353.68e-9,...
    'alpha_2',.1,...
    'de_3',1.2e6,...
    'tau_3',318.31e-6,...
    'alpha_3',.1,...
    'de_4',2.5e7,...
    'tau_4',2.274e-3,...
    'alpha_4',0,...
    'sigma',.2 ...
);

%muscle (along) %very rough estimation by Leo Medina
muscle_along = struct(...
    'e_inf',4,...    
    'de_1',70,... 
    'tau_1',7.23e-12,...
    'alpha_1',.1,...
    'de_2',500,... %different to across
    'tau_2',353.68e-9,...
    'alpha_2',.1,...
    'de_3',2e6,... 
    'tau_3',4e-4,... %different to across
    'alpha_3',.1,...
    'de_4',9e7,... %different to across
    'tau_4',1.8e-3,...
    'alpha_4',0,...
    'sigma',.2 ...
);

%skin (wet)
skin = struct(...
    'e_inf',4,...    
    'de_1',39,...
    'tau_1',7.96e-12,...
    'alpha_1',.1,...
    'de_2',280,...
    'tau_2',79.58e-9,...
    'alpha_2',0,...
    'de_3',3e4,...
    'tau_3',1.59e-6,...
    'alpha_3',.16,...
    'de_4',3e4,...
    'tau_4',1.592e-3,...
    'alpha_4',0.2,...
    'sigma',.0004 ...
);

%bone (cortical)
bone_cortical = struct(...
    'e_inf',2.5,...    
    'de_1',10,...
    'tau_1',13.6e-12,...
    'alpha_1',.2,...
    'de_2',180,...
    'tau_2',79.58e-9,...
    'alpha_2',.2,...
    'de_3',5e3,...
    'tau_3',159.15e-6,...
    'alpha_3',.2,...
    'de_4',1e5,...
    'tau_4',15.915e-3,...
    'alpha_4',0,...
    'sigma',.02 ...
);

%{
= struct(...
    'e_inf',,...    
    'de_1',,...
    'tau_1',,...
    'alpha_1',,...
    'de_2',,...
    'tau_2',,...
    'alpha_2',,...
    'de_3',,...
    'tau_3',,...
    'alpha_3',,...
    'de_4',,...
    'tau_4',,...
    'alpha_4',,...
    'sigma', ...
);
%}