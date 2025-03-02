%GENERALIZEDALPHA_TIMEADVANCE_DAMP class to Handle the time advancing scheme based 
%on Generalized Alpha scheme

%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri at epfl.ch>
%           Luca Rosafalco <luca.rosafalco at polimi.it>

classdef GeneralizedAlpha_TimeAdvance_Damped < handle
    
    properties (GetAccess = public, SetAccess = protected)
        M_U;
        M_dU;
        M_d2U;
        M_rhs_M;
        M_rhs_C;
        M_rhs_K;
        M_Csi;
        M_beta;
        M_gamma;
        M_alpha_f;
        M_alpha_m;
        M_timestep;
        M_stateSize;
    end
    
    methods (Access = public)
        
        %% Constructor
        function obj = GeneralizedAlpha_TimeAdvance_Damped( beta, gamma,...
                       alpha_m, alpha_f, timestep )
            obj.M_beta      = beta;
            obj.M_gamma     = gamma;
            obj.M_alpha_m   = alpha_m;
            obj.M_alpha_f   = alpha_f;
            obj.M_timestep  = timestep;
            obj.PrintCoefficients();
        end
        
        %% Initialize
        function obj = Initialize( obj, U, dU, d2U )
            if size(U,2) > 1
                error('Newmark_TimeAdvance.Initialize: incorrectSize of U')
            end
            obj.M_stateSize = size(U, 1);
            obj.M_U   = U;
            
            if size(dU,1) ~= obj.M_stateSize || size(dU,2) > 1
                error('Newmark_TimeAdvance.Initialize: incorrectSize of dU')
            end
            obj.M_dU  = dU;
            
            if size(d2U,1) ~= obj.M_stateSize || size(d2U,2) > 1
                error('Newmark_TimeAdvance.Initialize: incorrectSize of d2U')
            end
            
            obj.M_d2U = d2U;
            
            obj.M_Csi = 1 / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;
            
            obj.M_rhs_M = (1 - obj.M_alpha_m) / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - obj.M_alpha_m - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;
            
            obj.M_rhs_C = (obj.M_gamma * (1 - obj.M_alpha_f)) / (obj.M_beta * obj.M_timestep) * obj.M_U ...
                - ((obj.M_beta - obj.M_gamma*(1 - obj.M_alpha_f)) / obj.M_beta) * obj.M_dU ...
                - (1 - (obj.M_gamma / (2 * obj.M_beta))) * (1-obj.M_alpha_f) * obj.M_timestep * obj.M_d2U;
            
            obj.M_rhs_K = - obj.M_alpha_f * obj.M_U;            
        end
        
        %% Update
        function obj = Update( obj, U_np1 )
            
            d2U_np1        = 1 / ( obj.M_beta * obj.M_timestep^2) * U_np1 - obj.M_Csi;
            obj.M_dU       = obj.M_dU + obj.M_timestep * (obj.M_gamma * d2U_np1 + (1 - obj.M_gamma) * obj.M_d2U) ;
            obj.M_d2U      = d2U_np1;
            obj.M_U        = U_np1;
            
            obj.M_Csi = 1 / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;
            
            obj.M_rhs_M = (1 - obj.M_alpha_m) / (obj.M_beta * obj.M_timestep^2) * (obj.M_U + obj.M_timestep * obj.M_dU) ...
                + (1 - obj.M_alpha_m - 2*obj.M_beta)/(2*obj.M_beta) * obj.M_d2U;
            
            obj.M_rhs_C = (obj.M_gamma * (1 - obj.M_alpha_f)) / (obj.M_beta * obj.M_timestep) * obj.M_U ...
                - ((obj.M_beta - obj.M_gamma*(1 - obj.M_alpha_f)) / obj.M_beta) * obj.M_dU ...
                - (1 - (obj.M_gamma / (2 * obj.M_beta))) * (1-obj.M_alpha_f) * obj.M_timestep * obj.M_d2U;
            
            obj.M_rhs_K = - obj.M_alpha_f * obj.M_U;
        end
        
        %% Velocity and acceleration
        function [V_n, A_n] = VelAcc( obj )
            
            V_n = obj.M_dU;
            A_n = obj.M_d2U;
            
        end
        
        %% RhsContribute
        function Rhs_M = RhsContribute_M( obj )
            
            Rhs_M = obj.M_rhs_M;
            
        end
        
        %% RhsContribute damping
        function Rhs_C = RhsContribute_C( obj )
            
            Rhs_C = obj.M_rhs_C;
            
        end
        
        %% RhsContribute displacement
        function Rhs_K = RhsContribute_K( obj )
            
            Rhs_K = obj.M_rhs_K;
            
        end
        
        %% Mass coefficient
        function c = MassCoefficient( obj )
            
            c = (1 - obj.M_alpha_m) / ( obj.M_beta * obj.M_timestep^2) ;
            
        end
        
        %% Damping coefficient
        function c = DampingCoefficient( obj )
            
            c =  (obj.M_gamma * (1 - obj.M_alpha_f)) / (obj.M_beta * obj.M_timestep);
            
        end
        
        %% Stiffness coefficient
        function c = StiffnessCoefficient( obj )
            
            c =  1 - obj.M_alpha_f;
            
        end        
        
        %% PrintCoefficients
        function obj = PrintCoefficients( obj )
            
%             fprintf('\n=======================================');
%             fprintf('\nGeneralizedAlpha Time Advance Class');
%             fprintf('\n  beta     = %1.3f', obj.M_beta);
%             fprintf('\n  gamma    = %2.3f', obj.M_gamma);
%             fprintf('\n  alpha_m  = %1.3f', obj.M_alpha_m);
%             fprintf('\n  alpha_f  = %1.3f', obj.M_alpha_f);
%             fprintf('\n  timestep = %1.3e', obj.M_timestep);
%             fprintf('\n=======================================\n');
            
        end
    end
    
end