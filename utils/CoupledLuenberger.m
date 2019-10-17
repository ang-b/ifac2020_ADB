classdef CoupledLuenberger < handle
    %COUPLEDLUENBERGER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        A
        B
        C
        L
        Aij
        F
        xhat 
    end
    
    properties (Access = public)
        tSpan 
        x0
    end
    
        
    methods
        function self = CoupledLuenberger(A,B,C,Aij,L,tSpan,x0)
             self.A = A;
             self.B = B;
             self.C = C;
             self.Aij = Aij;
             self.L = L;
             self.F = A - L*C;
             
             if nargin == 7 && ~isempty(x0)
                 self.xhat = x0; 
                 self.x0 = x0;
             else
                 self.xhat = zeros(size(A,1), 1);
             end
             if nargin > 5  
                 self.tSpan = tSpan; 
             end

        end
        
        % 
        
        function self = setInitialCondition(self, x0)
            self.xhat = x0; 
            self.x0 = x0;
        end
        
        function y = yhat(self)
            y = self.C*self.xhat;
        end
        
        function self = estimate(self,u,y,xjd)
%             [~,intz] = ode45(@(t,x) self.stateEq(t,x,u,y,xjd,exo), [0 self.tSpan], self.xhat);
%             self.xhat = intz(end,:).';
              self.xhat = self.xhat + self.tSpan*self.stateEq([], self.xhat, u, y, xjd);
        end
        
    end
    
    methods (Access = private)
        function dz = stateEq(self,~,z,u,y,xjd)
            dz = self.F*z + self.B*u + self.L*y + self.Aij*xjd;
        end
    end
    
end

