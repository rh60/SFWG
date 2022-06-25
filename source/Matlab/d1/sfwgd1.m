classdef sfwgd1
    % Stabilizer Free Weak Galerkin Finite Element Method
    % solver class dim=1

    properties
        D
        domain
        x
        loc
        K, w        
        sol
        h, k
        monoms
        f
    end

    methods
        function obj = sfwgd1(problem,polyorder)
            [M,obj.monoms] = sfwgd1.MassMatrix(polyorder+1,1);
            obj.D = M^(-1);
            obj.h = 1;
            obj.k = polyorder;
            obj.domain = problem.domain;
            obj.f = problem.f;
            obj.w = zeros(polyorder+2, polyorder+3);
            obj.K = zeros(polyorder+3);
        end
        
        % u0 is polynomial
        function w = weakdiff(obj,u1,u0,u2)
            n = size(obj.D,1);
            b = zeros(n,1);
            q = 1;
            for i=1:n
                b(i) = - sfwgd1.L2product(u0,polyder(q),obj.h) + ...
                    u2*polyval(q,obj.h)-u1*polyval(q,0);
                q = [q 0];
            end
            w = obj.D*b;
        end    

        function obj = solve(obj,nele)
            obj.x = linspace( obj.domain(1),obj.domain(2),nele+1);
            obj.w(:,1) = obj.weakdiff(1,0,0);
            p=1;
            for i = 2:obj.k+2
                obj.w(:,i) = obj.weakdiff(0,p,0);
                p = [p 0];
            end
            obj.w(:,end) = obj.weakdiff(0,0,1);
            % local stiffness matrix
            n = size(obj.K,1);
            for i = 1:n
                for j = 1:i
                    obj.K(i,j) = sfwgd1.L2product(obj.w(:,i).',obj.w(:,j).',obj.h);
                end
            end
            obj.K=obj.K+tril(obj.K,-1)';
            N = nele*(obj.k+1) + nele-1;
            G = sparse(N,N);            
            G(1:n-1,1:n-1) = obj.K(2:n,2:n);
            ins = n-1;
            for i = 2:nele-1
                r=ins:ins+n-1;
                G(r,r) = G(r,r) + obj.K;
                ins = ins+n-1;
            end
            G(ins:end,ins:end) = G(ins:end,ins:end) + obj.K(1:n-1,1:n-1);
            rhs = zeros(N,1);
            ins = 1;
            for i=1:nele
                for j=1:obj.k+1
                    rhs(ins) = sfwgd1.L2product(obj.f,obj.monoms{j},obj.h);
                    ins = ins+1;
                end
                ins=ins+1;
            end
            obj.sol = G\rhs;
        end
    end

    methods (Static)
        function [M,mon]=MassMatrix(k,h)
            M = zeros(k+1);
            mon = cell(1,k+1); 
            mon{1}=1;
            for i=2:k+1
                mon{i} = [mon{i-1} 0];
            end
            for i=1:k+1               
                for j=1:i                    
                    M(i,j)=sfwgd1.L2product(mon{i},mon{j},h);                    
                end                
            end
            M=M+tril(M,-1)';
            M=flip(M,2);
        end

        function r=L2product(p,q,h)
            pint=polyint(conv(p,q));
            r=polyval(pint,h)-polyval(pint,0);
        end
    end
end

