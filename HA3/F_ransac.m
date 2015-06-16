function [Korrespondenzen_robust]=F_ransac(Korrespondenzen,varargin)
 %% Input parser
    P = inputParser;
    
    P.addOptional('epsilon', 0.1, @(x) isnumeric(x));
    % (default) = 0.1; ca. 10 Prozent falsch
    
    P.addOptional('p', 0.95, @(x) isnumeric(x) && x>0 && x<=1);
    % (default) = 0.95; zu 95% sicher dass kein Ausreißer in set

    P.addOptional('tolerance', 5, @(x) isnumeric(x));
    % (default) = 0;
    % Toleranz, dass Punkt noch zum Consensus-set gerechnet wird

    P.parse(varargin{:});

    %% Extrahiere die Variablen aus dem Input-Parser
    epsilon = P.Results.epsilon; 
    p = P.Results.p; 
    tolerance = P.Results.tolerance;
    
    %% Perform RanSaC
    %Anzahl der Daten k
    k=8;
    S=ceil(log(1-p)/log(1-(1-epsilon)^k));
    rand_index=zeros(S,k);
    consensus_index=zeros(S,size(Korrespondenzen,2));
    dist=zeros(S,size(Korrespondenzen,2));
    for s=1:S
        rand_index(s,:)=randi([1,size(Korrespondenzen,2)],[1,k]);
        rand_set=Korrespondenzen(:,rand_index(s,:));
        F=achtpunktalgorithmus(rand_set);
        dist(s,:)=compute_Sampson(F,Korrespondenzen);
        consensus_index(s,:)=dist(s,:)<=tolerance;
    end
    num_con=max(sum(consensus_index,2));
    maxIndex = find(sum(consensus_index,2) == num_con);
    if size(maxIndex,1) > 1
        sum_dist=sum(dist,2);
        [~,min_ind]=min(sum_dist(maxIndex));
        maxIndex=maxIndex(min_ind);
    end
    Korrespondenzen_robust=Korrespondenzen(:,consensus_index(maxIndex,:)==1);
end
function [dist] = compute_Sampson(F,Korrespondenzen)
    % Convert to homogenous coordinates
    korr_hom1=Korrespondenzen(1:2,:);
    korr_hom1(3,:)=ones(length(Korrespondenzen),1);
    korr_hom2=Korrespondenzen(3:4,:);
    korr_hom2(3,:)=ones(length(Korrespondenzen),1);
    % Compute parts of Sampson distance calculation
    e3_hat=apply_dachoperator([0;0;1]);
    norm1=e3_hat*F*korr_hom1;
    norm1=sum(norm1.*norm1,1);
    norm2=korr_hom2'*F*e3_hat;
    norm2=sum(norm2.*norm2,2)';
    norm=norm1+norm2;
    zaehler=sum((korr_hom2'*F)'.*korr_hom1,1).^2;
    dist=zaehler./norm;   
end

function [x_dach] = apply_dachoperator(x)
    x_dach=[0,-x(3),x(2);
            x(3),0,-x(1);
            -x(2),x(1),0];
end