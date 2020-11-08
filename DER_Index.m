%This function computes the (Duclos, Esteban & Ray, 2004) polarization index.
%An example of computed results for Mexico, Central America, Panama and Dominican Republic is available on the following IDB Technical Note from (Eggers & López-Marmolejo, 2020)(page 16):
%https://publications.iadb.org/es/polarizacion-instituciones-y-conflicto-una-aplicacion-mexico-el-istmo-centroamericano-y-republica
%In depth mathematical explanaitions and proofs are directly available on (Duclos, Esteban & Ray, 2004)

function [DER] = DER_Index(data,cY,alpha) %The inputs to get the index are the dataset (matrix), the column with the income variable (or whatever the variable you want to compute the index for) and the alpha parameter, where 0 is equivalent to the Gini index. 

    if alpha == 0 %If we use, alpha = 0, the index will be the same as the Gini index, which is computed using the function defined below.
        DER = GINI_Index(data,cY);
    else
        data = data(:,cY)/mean(data(:,cY)); %Normalizing the data so the scale has no effect on the results.
        data = sort(data);
        N = length(data);
        sigma = sqrt(var(data));
        h = 4.7*(1/sqrt(N))*sigma*alpha^(0.1); %Gaussian kernel optinal bandwidth
        h2 = h^2;
        Y = data;
        F = 1;

        for i = 1:N
            Y = data - data(i)*ones(N,1);
            Y = Y/h;
            Y = Y.^2;
            K = mean((1/(h*sqrt(2*pi)))*exp(-0.5*Y)); %Gaussian kernel.
            F = [F, K];
        end

        F = F(1,2:end);
        F = F.^alpha;

        mu = mean(data);
        A = zeros(1,N);

        for i = 1:N
            A(1,i) = mu + data(i)*((1/N)*(2*i - 1) - 1) - (1/N)*(2*sum(data(1:i-1)) + data(i));
        end

        DER = F*A';
        DER = DER/(2*N);
    end
    
end

function [Gini] = GINI_Index(data,cY)
    
    data = data(:,cY)/mean(data(:,cY));
    data = sort(data);
    N = length(data);
    sigma = sqrt(var(data));
    
    F = ones(1,N);
    
    mu = mean(data);
    A = zeros(1,N);
    
    for i = 1:N
        A(1,i) = mu + data(i)*((1/N)*(2*i - 1) - 1) - (1/N)*(2*sum(data(1:i-1)) + data(i));
    end
    
    Gini = F*A';
    Gini = Gini/(2*N);
    
end