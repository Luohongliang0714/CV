function [] = zeige_korrespondenzen( I1, I2, Korrespondenzen )
%Zeige Korrespondenzpunktpaare an
%% Zeige die Robusten Korrespondenzpunktpaare
    figure('name', 'Nach RanSac');
    imshow(uint8(I1))
    hold on
    plot(Korrespondenzen(1,:),Korrespondenzen(2,:),'r*')
    imshow(uint8(I2))
    alpha(0.5);
    hold on
    plot(Korrespondenzen(3,:),Korrespondenzen(4,:),'g*')
    for i=1:size(Korrespondenzen,2)
        hold on
        x_1 = [Korrespondenzen(1,i), Korrespondenzen(3,i)];
        x_2 = [Korrespondenzen(2,i), Korrespondenzen(4,i)];
        line(x_1,x_2);
    end
    hold off
end

