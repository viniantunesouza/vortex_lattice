function [gama] = wessinger(N)
    V=1 %V no infinito, em m/s
    AR=6 %Alongamento
    alfa=1 %Define alfa como 1 rad 
    l=AR/(2*N) %Distância entre a linha de simetria e o lado da ferradura
    x=-0.5 %Distância constante entre o ponto de controle e de aplicação, em relação a qualquer filamento
    z=0 %Valor de z para a asa plana
    
    for i=1:N
        for j=1:N
            y(i,j)=2*l*(i-j) %Calcula as posições de influência
        end
    end
    
    for i=1:N
        for j=1:N
            w(i,j)=((-x/(4*pi*(x^2+z^2)))*(((y(i,j)+l)/sqrt(x^2+(y(i,j)+l)^2+z^2))-((y(i,j)-l)/sqrt(x^2+(y(i,j)-l)^2+z^2))))-(((y(i,j)-l)/(4*pi*(z^2+(y(i,j)-l)^2)))*(1-(x/sqrt(x^2+(y(i,j)-l)^2+z^2))))+(((y(i,j)+l)/(4*pi*(z^2+(y(i,j)+l)^2)))*(1-(x/sqrt(x^2+(y(i,j)+l)^2+z^2)))) %Calcula a matriz de velocidade induzida
        end
    end
    
    B=V*(-alfa*ones(N,1)) %Calcula a matriz B
    circulacao = linsolve(w,B) %Resolve o sistema de equações para encontrar a circulação
    
    for i=1:(N/2) %os coeficientes são dados para metade da asa
        gama(i)=-circulacao(i+(N/2))/(AR*V) %calcula o gama
    end
    
    for i=1:(N/2)
        eta(i)=(l+(2*l*(i-1)))/(AR/2)%calcula o eta
    end
    
    figure()
    scatter(eta,gama)%plota um gráfico de pontos
    grid on
    title('Gráfico de pontos da distribuição de circulação')
    xlabel('\eta')
    ylabel('\gamma')
    figure()
    plot(eta,gama)%plota um gráfico de linhas
    grid on
    title('Gráfico de linhas da distribuição de circulação')
    xlabel('\eta')
    ylabel('\gamma')
    
end
