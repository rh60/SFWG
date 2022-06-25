classdef sfwgd1
    % Stabilizer Free Weak Galerkin Finite Element Method
    % solver class dim=1

    properties
        M, D, k, h, monoms, domain, f        
        x, K, w
        sol, lift                
    end

    methods
        function obj = sfwgd1(problem,polyorder)
            [obj.M,obj.monoms] = sfwgd1.MassMatrix(polyorder+1,1);
            obj.D = obj.M^(-1);
            obj.h = 1;
            obj.k = polyorder;
            obj.domain = problem.domain;
            obj.f = problem.f;
            obj.w = zeros(polyorder+2, polyorder+3);
            obj.K = zeros(polyorder+3);
        end

        % u0 is polynomial
        function w = weakdiff(obj,u1,u0,u2)
            L2 = @sfwgd1.L2product;
            n = size(obj.D,1);
            b = zeros(n,1);            
            for i=1:n
                q = obj.monoms{i};
                b(i) = - L2(u0,polyder(q),obj.h) + ...
                    u2*polyval(q,obj.h)-u1*polyval(q,0);                
            end
            w = obj.D*b;
        end

        function obj = solve(obj,nele)
            L2 = @sfwgd1.L2product;
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
                    obj.K(i,j) = L2(obj.w(:,i).',obj.w(:,j).',obj.h);
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
                    rhs(ins) = L2(obj.f,obj.monoms{j+1},obj.h); % !!!!!!!!
                    ins = ins+1;
                end
                ins=ins+1;
            end
            obj.sol = [0; G\rhs; 0];
            obj = Lift(obj);
        end

        function obj = Lift(obj)
            L2 = @sfwgd1.L2product;
            nele = numel(obj.x)-1;
            m = obj.monoms{1};
            mon = [ [m 0] obj.monoms ];
            Dmon = cellfun(@polyder, mon(1:end-1), 'UniformOutput',false);
            n = numel(Dmon); 
            A = zeros(n); b = zeros(n,1);
            for i=1:n
                for j=1:n
                    A(i,j) = L2(Dmon{i},Dmon{j},obj.h);
                end
            end
            obj.lift = zeros(nele,obj.k+3);
            for i=1:nele
                u = obj.sol((i-1)*n+1:i*n+1);
                w=obj.w*u;
                % w=obj.weakdiff(obj.sol(1),obj.sol(2),obj.sol(3))
                for j=1:n
                    b(j) = L2(w',Dmon{j},obj.h);
                end
                uh = A\b;
                uh = [uh' 0];
                u0 = u(2:end-1);
                uh(end) = (L2(u0,1,obj.h)-L2(uh,1,obj.h))/L2(1,1,obj.h);
                obj.lift(i,:)=uh;
            end            
        end
    end

    methods (Static)
        function [M,mon] = MassMatrix(k,h)
            L2 = @sfwgd1.L2product;
            M = zeros(k+1);
            mon = cell(1,k+1);
            mon{1}=1;
            for i=2:k+1
                mon{i} = [mon{i-1} 0];
            end
            mon=flip(mon);
            for i=1:k+1
                for j=1:i
                    M(i,j) = L2(mon{i},mon{j},h);
                end
            end
            M=M+tril(M,-1)';
            %M=flip(M,2);
        end

        function r=L2product(p,q,h)
            pint=polyint(conv(p,q));
            r=polyval(pint,h)-polyval(pint,0);
        end
    end
end

