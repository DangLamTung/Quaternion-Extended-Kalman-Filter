log_file_full = load('banxoay_480euler_500Hz/final/log_25.TXT');
starting_row = 1;
end_row = 100000;
log_file_test  = log_file_full(starting_row:end_row,:);

N = max(size(log_file_test)); % Size of tested sample
T = 0.002; % Sampling time
time = 0:T:(N-1)*T; % Time array

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change unit of sensor values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gyro = log_file_test(:,4:6)*pi/180000;        % unit: rad/s
acc  = log_file_test(:,7:9)/10000;       % unit: m/s^2
mag  = log_file_test(:,10:12)/10000;          % unit: gauss


t = 0.002
x =               [1                       % gyrox
                  0                       % gyroy
                  0                       % gyroz
                  0];                     %
P = eye(4)
Q = [   0.00001  0  0  0  ;
        0  0.00001  0  0  ;
        0  0  0.00001  0  ;
        0  0  0  0.00001  ;
         ];
R = [  0.0000178 ,0 ,0  ;
       0 ,0.0000134 ,0  ;
       0 ,0, 0.0000192  ]
R_full = [  0.0000178  0  0 0 0 0 ;
       0  0.0000134  0 0 0 0 ;
       0  0  0.0000192 0 0 0 ;
       0  0  0  0.00000372  0  0;
       0  0  0  0  0.00000231  0;
       0  0  0  0  0  0.00000155];
roll_estimated = zeros(N,1);
pitch_estimated = zeros(N,1);
yaw_estimated = zeros(N,1);
Gx_estimated = zeros(N,1);
Gy_estimated = zeros(N,1);
Gz_estimated = zeros(N,1);

k = 1

while k~=N
    
w = [gyro(k, 1), gyro(k, 2), gyro(k, 3)];
y =  [acc(k, 1), acc(k, 2), acc(k, 3)];
mag_norm = sqrt(mag(k,1)^2 + mag(k,2)^2 + mag(k,3)^2);

y_full = [acc(k, 1), acc(k, 2), acc(k, 3),mag(k,1)/mag_norm,mag(k,2)/mag_norm,mag(k,3)/mag_norm];

    

  

s_trans = [[-x(2),-x(3),-x(4)];
           [x(1),-x(4),x(3)];
           [x(4), x(1), -x(2)];
           [-x(3), x(2),x(1)]];
s_up = -s_trans*0.5*t;

Ca = [[-x(3) x(4) -x(1) x(2)];
      [x(2) x(1) x(4) x(3)];
      [x(1) -x(2) -x(3) x(4)]];
 Cm = [[x(4) x(3) x(2) x(1)];
     [x(1) -x(2) x(3) -x(4)];
     [-x(2) -x(1) x(4) x(3)]];
   H =[ [-2*Ca]];
   H_full =[ [-2*Ca ]; [2*Cm]];

    A = [[eye(4)]];
        
      B = [s_trans*0.5*t];
%predict route

     x = A*x + B*transpose(w);
     P = A*P*transpose(A) +Q ;
% % %quaternion rotation matrix
% % 
    K = P*transpose(H)/(H*P*transpose(H) + R);
    K_full = P*transpose(H_full)/(H_full*P*transpose(H_full) + R_full);
%       y;
    
    s = (transpose(y) - H*x);
    s_full = (transpose(y_full) - H_full*x);
    m = [ transpose(y) H*x];
     x = x + K_full*s_full;
     P = (eye(4) - K_full*H_full)*P ;
    norm = sqrt(x(1)^2 + x(2)^2 + x(3)^2 + x(4)^2);
    x(1) = x(1)/norm;
    x(2) = x(2)/norm;
    x(3) = x(3)/norm;
    x(4) = x(4)/norm;
    [roll_estimated(k), pitch_estimated(k), yaw_estimated(k)] = IMU_Quaternion2Euler(x(2), x(3), x(4), x(1));
       if(roll_estimated(k) > 0)
           roll_estimated(k) = roll_estimated(k) - 180;
       
       else
         roll_estimated(k) = roll_estimated(k) + 180;
       end
       pitch_estimated(k) = -pitch_estimated(k);
       yaw_estimated(k) = yaw_estimated(k);
     
%     if(yaw_estimated(k) > 180)
%         roll_estimated(k) = roll_estimated(k) - 180;
%     end
%     if(yaw_estimated(k) < -180)
%         roll_estimated(k) = roll_estimated(k) + 180;
%     end
    k = k+1;   
    Gx_estimated(k)     = Gx_estimated(k-1) +gyro(k, 1)*t;
    Gy_estimated(k)     = Gy_estimated(k-1) +gyro(k, 2)*t;
    Gz_estimated(k)     = Gz_estimated(k-1) +gyro(k, 3)*t;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot data and error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read and extract reference data from commercial IMU and encoder
IMU = log_file_test(:,1:3)/1000;
ENC = (log_file_test(:,16:18) - 32000) * 360 / 10000;
ENC(:,2) = -ENC(:,2);
% ENC(:,3) = -ENC(:,3);

offset = IMU(1,:) - ENC(1,:);
error  = IMU - ENC;
error(:,1) = error(:,1) - offset(1); 
error(:,2) = error(:,2) - offset(2);
error(:,3) = error(:,3) - offset(3);

% RMS_error = [norm(error(:,1))/sqrt(N);
%              norm(error(:,2))/sqrt(N);
%              norm(error(:,3))/sqrt(N)]

gain = 100;
figure;
subplot(3,1,1);
hold on;grid on;
 plot(time,ENC(:,1),'r');
 plot(time,IMU(:,1),'b');
 plot(time,error(:,1)*gain,'g');
 plot(time,log_file_test(:,4)/10000,'black');
 legend('enc','imu','err*100','rate/10');
 title('Roll');
% % 
 subplot(3,1,2);
 hold on;grid on;
 plot(time,ENC(:,2),'r');
 plot(time,IMU(:,2),'b');
 plot(time,error(:,2)*gain,'g');
 legend('enc','imu','err*100','rate/10');
 title('Pitch');

subplot(3,1,3);
 hold on;grid on;
 plot(time,ENC(:,3),'r');
 plot(time,IMU(:,3),'b');
 plot(time,error(:,3)*gain,'g');
 plot(time,log_file_test(:,6)/10000,'black');
 legend('enc','imu','err*100','rate/10');
 title('Yaw');
 
% % % Plot estimated angular velocity and measured
 figure;
 subplot(3,1,1);
hold on;grid on;
plot(time,gyro(:,1),'r');
  plot(time,Gx_estimated,'b');
  legend('gyro measured', 'estimated');
  title('x-axis');
  % % % 
  subplot(3,1,2);
  hold on;grid on;
  plot(time,gyro(:,2),'r');
  plot(time,Gy_estimated,'b');
 legend('measured', 'estimated');
  title('y-axis');
  
  subplot(3,1,3);
  hold on;grid on;
 plot(time,gyro(:,3),'r');
  plot(time,Gz_estimated,'b');
  legend('measured', 'estimated');
  title('z-axis');

% Plot estimated orientation and error
error = zeros(N,3);
error = ENC - [roll_estimated, pitch_estimated, yaw_estimated];
figure;
subplot(3,1,1);
hold on;grid on;
plot(time,ENC(:,1),'r');
 plot(time,roll_estimated,'b');
%  plot(time,error(:,1),'g');
legend('enc', 'estimated', 'error');
title('Roll');

subplot(3,1,2);
hold on;
grid on;
plot(time,ENC(:,2),'r');
 plot(time,pitch_estimated,'b');
%  plot(time,error(:,2),'g');
legend('enc', 'estimated', 'error');
title('Pitch');

subplot(3,1,3);
hold on;grid on;
plot(time,ENC(:,3),'r');
 plot(time, yaw_estimated,'b');
%  plot(time,error(:,3),'g');
legend('enc', 'estimated', 'error');
title('Yaw');

