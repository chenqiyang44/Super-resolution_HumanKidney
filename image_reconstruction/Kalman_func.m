function [Q_loc_estimate] = Kalman_func(Q_loc_meas)
    duration = size(Q_loc_meas,1);  %how long the Quail flies
    %but we'll assume he's just repeatedly sampling over time at a fixed interval

    %% Define update equations (Coefficent matrices): A physics based model for where we expect the Quail to be [state transition (state + velocity)] + [input control (acceleration)]
    A = [1 0 1 0;0 1 0 1;0 0 1 0;0 0 0 1]; % state transition matrix:  expected flight of the Quail (state prediction)
    %B = [dt^2/2; dt]; %input control matrix:  expected effect of the input accceleration on the state.
    C = [1 0 0 0;0 1 0 0]; % measurement matrix: the expected measurement given the predicted state (likelihood)
    %since we are only measuring position (too hard for the ninja to calculate speed), we set the velocity variable to
    %zero.

    %% define main variables
    Q_pos_z = Q_loc_meas(2,1);%initized state--it has two components: [position; velocity] of the Quail
    Q_pos_x = Q_loc_meas(2,2);
    Q_vel_z = Q_loc_meas(2,1)-Q_loc_meas(1,1); 
    Q_vel_x = Q_loc_meas(2,2)-Q_loc_meas(1,2);
    Q_estimate = [Q_pos_z Q_pos_x Q_vel_z Q_vel_x]';  %x_estimate of initial location estimation of where the Quail is (what we are updating)
    %Observ_var = round(40/(1540/15.625e6/8*1e6));  %measurement noise: How mask-blinded is the Ninja (stdv of location, in meters)
    %Predict_var = round(40/(1540/15.625e6/8*1e6));
    Ez = [16 0;0 16];% Ez convert the measurement noise (stdv) into covariance matrix
    Ex = 1 * eye(4);% Ex convert the process noise (stdv) into covariance matrix
    P = Ex; % estimate of initial Quail position variance (covariance matrix)

    %% Do kalman filtering
    %initize estimation variables
    %Q_loc_estimate = [];
    Q_loc_estimate = [Q_loc_meas(1,1) Q_loc_meas(1,2)]; %  Quail position estimate
    P_estimate = P;
    P_mag_estimate = [];
    predic_state = [];
    predic_var = [];
    for t = 2:size(Q_loc_meas,1)
        % Predict next state of the quail with the last state and predicted motion.
        Q_estimate = A * Q_estimate;
        predic_state = [predic_state; Q_estimate(1:2)'] ;
        %predict next covariance
        P = A * P * A' + Ex;
        predic_var = [predic_var; P] ;
        % predicted Ninja measurement covariance
        % Kalman Gain
        K = P*C'*inv(C*P*C'+Ez);
        % Update the state estimate.
        Q_estimate = Q_estimate + K * (Q_loc_meas(t,:)' - C * Q_estimate);
        % update covariance estimation.
        P =  (eye(4)-K*C)*P;
        %Store for plotting
        Q_loc_estimate = [Q_loc_estimate; Q_estimate(1:2)'];
        %vel_estimate = [vel_estimate; Q_estimate(2)];
        %P_mag_estimate = [P_mag_estimate; P(1)];

    end

    % Plot the results
%     figure(1);clf
%     plot(Q_loc_meas(:,2), Q_loc_meas(:, 1), 'Color', 'g')
%     hold on
%     plot(Q_loc_estimate(:,2), Q_loc_estimate(:, 1), 'Color', 'r')
end
