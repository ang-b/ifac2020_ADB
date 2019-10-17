classdef UIO < handle
    %UIO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        A
        A1
        B
        C
        E
        F
        H
        K
        K1
        K2
        T
        feasibility
        xhat
    end
    
    properties (Access = public)
        tSpan
        z0
    end
    
    properties (Access = private)
        z    
    end
    
    methods
        % constructor
        function self = UIO(A,B,C,Aij,K1,tSpan)
            self.A = A;
            self.B = B;
            self.C = C;
            self.E = Aij;
            
            [self.feasibility, lrE] = checkExistence(self);
            if self.feasibility
                H = lrE * ((C*lrE)'*C*lrE)^(-1) * (C*lrE)';
                
                self.A1 = A - H*C*A;
                self.H = H;  
                self.E = lrE;
                
                if nargin > 4 && ~isempty(K1) 
                    self.setGainK1(K1); 
                end
                if nargin > 5 && ~isempty(tSpan) 
                    self.tSpan = tSpan; 
                end
                
                self.z = zeros(size(A,1), 1);
                self.z0 = self.z;
                
            else
                warning('Cannot design UIO for given system');
            end
        end
        
        % member methods      
      
        function self = setGainK1(self, K1)
            self.K1 = K1; 
            F_uio = self.A1 - K1 * self.C;

            HC  = self.H*self.C;
            self.T = eye(size(HC)) - HC;

            self.K2 = F_uio * self.H;
            self.F = F_uio;
            self.K = K1 + self.K2;
        end
        
        function self = setInitialCondition(self, z)
            % TODO: add dimensional checks etc
            self.z = z;
            self.z0 = z;
        end

        function self = assignFPoles(self, poles) % UNTESTED
            n = size(self.A, 1);
            poles = poles(:);
            if (length(poles) > n)
                warning('Too many poles specified, considering the first %d', n);
                poles = poles(1:n);
            elseif (length(poles) < n)
                error('Need to specify at least %d poles', n);   
            end
            
            K1 = place(self.A1.', self.C.', poles).';
            self.setGainK1(K1);
        end
        
        function self = estimate(self,u,y)
%             [~, intz] = ode45(@(t,x) self.stateEq(t,x,u,y),[0 self.tSpan], self.z);
%             self.z = intz(end,:).';
            self.z = self.z + self.tSpan*self.stateEq([], self.z, u, y);
            self.xhat  = self.z + self.H*y;
        end
    end

    methods (Access = private)
        function [cond, lrE] = checkExistence(self) 
            % perform a low rank factorisation first
            [Usv,Ssv,~] = svd(self.E);
            rankE = sum(diag(Ssv) ~= 0);
            if (rankE < size(self.A,1))
                lrE = Usv(:,1:rankE)*sqrtm(Ssv(1:rankE, 1:rankE));
            else
                lrE = self.E;
            end
            C = self.C;
            
            % rank condition check
            if (rank(C*lrE) ~= rank(lrE))
                cond = false;
            else
                cond = true;
            end
        end
        
        function dz = stateEq(self,~,z,u,y)
            dz = self.F*z + self.T*self.B*u + self.K*y;
        end
    end
    
end

