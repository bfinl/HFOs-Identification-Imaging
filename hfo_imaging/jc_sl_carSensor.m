%% Source imaging
% This function performs common average sreference for the input signals.
% This input signals are assumed to have dimensions as nChannel x nTime.
%
%--------------------------------------------------------------------
% Zhengxiang Cai
% 2020.08.21
% Document and commit for repository.


function phi_car = jc_sl_carSensor(phi)

% centralize the sensor data
% do common average reference

[nChan,~] = size(phi);

car = mean(phi,1);
phi_car = phi - repmat(car,nChan,1);

end