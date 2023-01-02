function [CL3,CM3,xcp] = vortex_lattice(N)
    V=1 %V no infinito, em m/s
    c=1 %corda no perfil, em m
    rho=1.2 %Densidade do ar ao nível do mar, em kg/m^3
    alfa=-14:14 %Define alfa nos valores de -14° a 14° 
    for l=1:length(alfa) %Define para qual alfa vai ser rodada a rotina
        alfa_r=deg2rad(alfa(l)) %Transforma o alfa de graus para radianos
    
        for i=1:N
            xj(i) = [(1+4*(i-1))/4]*(1/N) %Cria a matriz dos pontos de aplicação
            xk(i) = [(3+4*(i-1))/4]*(1/N) %Cria a matriz dos pontos de controle
        end
        
        for k=1:N
            for j=1:N
                A(k,j)=-(1/(2*pi))*(1/(xk(k)-xj(j))) %Cria a matriz dos coeficientes de influência
            end
        end
    
        for i=1:N
            if xk(i) <= 0.4
                z(i)=0.25*(0.8-2*xk(i)) %Cria a matriz das derivadas do câmber, caso xk <= 0.4
            else
                z(i)=0.111*(0.8-2*xk(i)) %Cria a matriz das derivadas do câmber, caso xk > 0.4
            end
        end
        
        B=transpose(z-alfa_r) %Calcula a matriz B
        circulacao = linsolve(A,B) %Resolve o sistema de equações para encontrar a circulação
        
        coefl(l)=2*sum(circulacao) %Calcula o coeficiente de sustentação para cada alfa
        
        for i=1:N
            M(i)=circulacao(i)*(1/4-xj(i)) %Calcula o momento de cada circulação em relação ao ponto a 1/4 da corda
        end
        
        coefm(l)=2*sum(M) %Calcula o coeficiente de momento para cada alfa

    end

    plot(alfa,coefl) %Plota a curva Cl alfa
    CL3=coefl(18) %Define o coeficiente de sustentação para alfa=3°
    CM3=coefm(18) %Define o coeficiente de momento a 1/4 da corda para alfa=3°
    xcp=1/4-(CM3/CL3) %Calcula a posição do centro de pressão para alfa=3°

end
