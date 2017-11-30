function [msd, lags] = polymer_msd(data, box_size, frame_rate)
%input data is an Nx3 position matrix

N = length(data);
lags = 1:N/3;

%need to account for periodic boundaries!!!
for i = 1:length(lags)
    dispx = zeros(N-lags(i), 1);
    dispy = zeros(N-lags(i), 1);
    dispz = zeros(N-lags(i), 1);
    for j = 1:N-lags(i)
        dispx(j) = data(j+lags(i),1) - data(j,1);
        dispy(j) = data(j+lags(i),2) - data(j,2);
        dispz(j) = data(j+lags(i),3) - data(j,3);
        if dispx(j) >= box_size/2
            dispx(j) = dispx(j) - box_size;
        elseif dispx(j) <= -box_size/2
            dispx(j) = dispx(j) + box_size;
        end
        if dispy(j) >= box_size/2
            dispy(j) = dispy(j) - box_size;
        elseif dispy(j) <= -box_size/2
            dispy(j) = dispy(j) + box_size;
        end
        if dispz(j) >= box_size/2
            dispz(j) = dispz(j) - box_size;
        elseif dispz(j) <= -box_size/2
            dispz(j) = dispz(j) + box_size;
        end
    end
    
    lags(i) = lags(i)/frame_rate;
    squared_disp = dispx.^2 + dispy.^2 + dispz.^2;
    msd(i) = mean(squared_disp);

end
