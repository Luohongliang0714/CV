function [] = show_back_proj( I, x2, x2_back)
%% Zeige Rückprojektion
    figure('name', 'Backprojection to Second camera frame');
    imshow(uint8(I))
    hold on
    plot(x2(1,:),x2(2,:),'r*')
    plot(x2_back(1,:),x2_back(2,:),'g*')
    for i=1:size(x2,2)
        hold on
        x = [x2(1,i), x2_back(1,i)];
        y = [x2(2,i), x2_back(2,i)];
        line(x,y);
    end
    hold off
end

