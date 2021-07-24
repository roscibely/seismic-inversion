function [refl,dp, gnorm, dx] = gradd(d,w,alpha,maxiter, tol)
%Exemplo: [refl,dp, gnorm] = gradd(dado_sismico,w,0.08,100, 10^(-3));
%         [refl3,dp3, gnorm] = gradd(dado_sismico,w,0.1,200, 10^(-7));    
%[r20,d20, norm200] = gradd(dado_sismico,w,0.009,190, 10^(-7));
    E = 0.0001;
    nw = length(w);
    [nt,ntraces] = size(d);
    d = [d;zeros(nw-1,ntraces)];
    W = convmtx(w,nt); %retorna a matriz de convolução da wavelet w
    % Tolerancia
    tol = 1e-6;
    % Numero maximo de interação 
    %maxiter = 1000;
    %perturbação mínima permitida
    dxmin = 1e-6;
    %tamanho do passo 
    %alpha = 0.01;
    % inicializaçao da norma do gradiente, reflectividade, contador niter, e
    % pertubação
    gnorm = inf; r = zeros(nt,1); niter = 0; dx = inf; 
    X3 = zeros(nt,nt); 



    for itrace = 1:ntraces;
        for i=1:nt
            X3(nt,nt)= (((E./(abs(r(i))'*abs(r(i))+E)^(3/2)))); %X3(r)
        end
        s = d(:,itrace);
        psi = 0.5; 
        a(itrace) = 1*psi^(itrace-1); %Parametro de regularização
 
        H = W'*W + a(itrace)*X3; %Matriz Hessian : H = W'W + alpha*X3(r)
        %Ho = inv(H); %Matriz Hessiana Inversa
 
        % Gradiente Descendente:
        while(gnorm>=tol & (niter <= maxiter & dx >= dxmin)) %Mudar o operado and por or
            % Calcular o gradiente:
            F = (r./sqrt((abs(r)'*abs(r)+E)));
            g = W'*W*r - W'*s + a(itrace)*F; %Add a função da regra geometrica
            gnorm = norm(g);
            %Atualizar o valor de r:
            B1 = 0.001; 
            B2 = 0.003; 
            h= B1*(g'*g)/(g'*(H*g)) + B2*(g'*(H*g))/(g'*(H'*H*g)); %parâmetros critérios de Rayleigh onde B1 e B2 os constantes atribuida pelo usuario

            rnew = r - h*g;
            %rnew = r - alpha*g;
            % checa o valor de r
            if ~isfinite(rnew)
                display(['Numero de interação: ' num2str(niter)])
                error('r é inf ou NaN')
            end
            % atualizar
            niter = niter + 1;
            dx = norm(rnew-r);
            r = rnew;   
        end

    refl(:,itrace) = r;
    niter =0;
    end


    n2 = floor(nw/2);
    temp = refl;
    dp = conv2(refl,w);
    dp = dp(1:nt,:);
    refl = [zeros(n2,ntraces);temp(1:nt-n2,:)];

 return; 
