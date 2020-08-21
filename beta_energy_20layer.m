clear all
close all

path(pathdef)

CALC_SLOPE = 1;
CALC_SLOPE_BY_WAVENUMBER = 0;


% Growth rates at the "kink"
% 0.0    2.40e-8
% 0.1    2.09e-8
% 0.2    1.92e-8
% 0.3    1.60e-8
% 0.4    1.11e-8


directory = '/project/oce/bachman/MITgcm/dedalus/two_layer_beta/';
output_dir = '21';
files = dir(fullfile(directory, strcat(output_dir, '_1*')));
cd(directory)

deltaU = 0.0013;
Ld = 9e4;
nx = 512;
dx = Ld / 8; %1/(gridratio / Ld);
layers = 20;
beta1 = (str2num(output_dir)-21) * 3.08641975308642e-12 / 2;
beta = [(beta1 + deltaU /(4.5e3^2)) beta1*ones(1,17) (beta1 - deltaU /(4.5e3^2))]


snaps = strcat(directory, files(1).name, strcat('/snapshots/snapshots.h5'));

%h5info(snaps, '/') 
time = h5read(snaps,'/scales/sim_time');
y = h5read(snaps,'/scales/y/1.0');
x = h5read(snaps,'/scales/x/1.0');

for k = 1:layers
  n = num2str(k);
  eval([strcat('psi', n, ' = [];')])
  eval([strcat('v', n, ' = [];')])
  eval([strcat('q', n, ' = [];')])
end

h5info(snaps, '/') 
time = h5read(snaps,'/scales/sim_time');
y = h5read(snaps,'/scales/y/1.0');
x = h5read(snaps,'/scales/x/1.0');

for k = 1:layers
    n = num2str(k);
 eval([strcat('psi', n, ' = cat(3, psi', n, ', squeeze(h5read(snaps, ''/tasks/psi', n, ''')));')]);
% eval([strcat('v', n, ' = cat(3, v', n, ', squeeze(h5read(snaps, ''/tasks/v', n, ''')));')]);
 eval([strcat('q', n, ' = cat(3, q', n, ', squeeze(h5read(snaps, ''/tasks/q', n, ''')));')]); 
end

KE = fft2(psi10) .* conj(fft2(psi10));
tmp = KE(1,:,end);
yay = find(tmp > 0.001 * max(tmp));

  for k = 1:layers
    n = num2str(k);
    eval([strcat('psie', n, ' = cat(1, psi', n, '(end,:,:), psi', n, ', psi', n, '(1,:,:));')])
    eval([strcat('psie', n, ' = cat(2, psie', n, '(:,end,:), psie', n, ', psie', n, '(:,1,:));')])
  end
  'Finished making psi.'
 
  for k = 1:layers
    n = num2str(k);
    eval([strcat('v', n, '= 0.5*(psie', n, '(2:end-1,3:end,:) - psie', n, '(2:end-1,1:end-2,:)) / dx;')])
    eval([strcat('u', n, '= -0.5*(psie', n, '(3:end,2:end-1,:) - psie', n, '(1:end-2,2:end-1,:)) / dx;')])
  end    

'Done reading output.'

%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SLOPE %%%%%%%%%%%%%

if CALC_SLOPE == 1
for k = 1:layers-1
  n = num2str(k);
  eval([strcat('psidiff', n, ' = [];')])
  eval([strcat('vmean', n, ' = [];')])
  eval([strcat('qmean', n, ' = [];')])  
end

for k = 1:layers-1
 n = num2str(k);
 np1 = num2str(k+1);
 eval([strcat('psidiff', n, ' = psi', n, ' - psi', np1, ';')]);
 eval([strcat('vmean', n, ' = 0.5*(v', n, ' + v', np1, ');')]);
 eval([strcat('qmean', n, ' = 0.5*(q', n, ' + q', np1, ');')]); 
end


for k = 1:layers-1
  n = num2str(k);
  eval([strcat('mean_psidiff2_', n, ' = [];')])
  eval([strcat('mean_vpsidiff', n, ' = [];')])
  eval([strcat('mean_qpsidiff', n, ' = [];')])  
  eval([strcat('mean_vq', n, ' = [];')])  
end
dt = diff(time);
for k = 1:layers-1
 n = num2str(k);
 eval([strcat('mean_psidiff2_', n, ' = squeeze(mean(squeeze(mean(psidiff', n, '.^2 ,1)),1));')]);
 eval([strcat('mean_vpsidiff', n, ' = squeeze(mean(squeeze(mean(vmean', n, ' .* psidiff', n, ',1)),1));')]);  % v'b'
 eval([strcat('mean_qpsidiff', n, ' = squeeze(mean(squeeze(mean(qmean', n, ' .* psidiff', n, ',1)),1));')]);  % b'q'
 eval([strcat('mean_vq', n, ' = squeeze(mean(squeeze(mean(vmean', n, ' .* qmean', n, ',1)),1));')]); % v'q'
end
yayayaya
for k = 1:layers-1
  n = num2str(k);
  eval([strcat('num', n, ' = [];')])
  eval([strcat('denom', n, ' = [];')])
  eval([strcat('num21', n, ' = [];')])
  eval([strcat('num22', n, ' = [];')])  
  eval([strcat('denom2', n, ' = [];')])  
end
for k = 1:layers-1
 n = num2str(k);
 eval([strcat('num', n, ' = diff(mean_psidiff2_', n, ') ./ transpose(dt);')]);
 eval([strcat('num', n, ' = [0 num', n, '];')]);
 eval([strcat('denom', n, ' = mean_vpsidiff', n, ';')]); 
 eval([strcat('num21', n, ' = diff(mean_qpsidiff', n, ') ./ transpose(dt);')]);  
 eval([strcat('num21', n, ' = [0 num21', n, '];')]); 
 eval([strcat('num22', n, ' = beta(', n,') * mean_vpsidiff', n, ';')]); 
 eval([strcat('denom2', n, ' = mean_vq', n, ';')]);  
end

for k = 1:layers-1
  n = num2str(k);
  eval([strcat('C', n, ' = [];')])
  eval([strcat('C2', n, ' = [];')])  
  eval([strcat('C3', n, ' = [];')])   
  eval([strcat('C4', n, ' = [];')])      
end
for k = 1:layers-1
 n = num2str(k);    
 eval([strcat('C', n, ' = 1 - (1 ./ (2*deltaU)) * num', n, './ denom', n, ';')])
 eval([strcat('C2', n, ' = 1 - (1 ./ (deltaU)) * (num21', n, ' + num22', n, ')./ denom2', n, ';')]) 
 eval([strcat('C3', n, ' = -(1 ./ (deltaU)) *  num22', n, './ denom2', n, ';')])   
 eval([strcat('C4', n, ' = -(1 ./ (deltaU)) *  num21', n, './ denom2', n, ';')])  
 eval([strcat('vpbp', n, ' = denom', n, ';')])
end

figure
C_sum = 0;
C2_sum = 0;
C3_sum = 0;
C4_sum = 0;
count = 0;
for k = 2:layers-2
    count = count + 1;
    n = num2str(k);
    nm1 = num2str(k-1);
    subplot(2,2,1)
    eval([strcat('plot(C', n, ')')])
    eval([strcat('C_sum = C_sum + C', n, ';')])
    eval([strcat('C_array(', nm1,') = mean(C', n, '(100:end));')])
    hold on
    xlim([20 150])
    ylim([0 1])
    
    subplot(2,2,2)
    eval([strcat('plot(C2', n, ')')])
    eval([strcat('C2_sum = C2_sum + C2', n, ';')])
    eval([strcat('C2_array(', nm1,') = mean(C2', n, '(100:end));')])
    hold on
    xlim([20 100])
    ylim([0 5]) 
    
    % Pos def
    subplot(2,2,3)
    eval([strcat('plot(C3', n, ')')])
    eval([strcat('C3_sum = C3_sum + C3', n, ';')])
    eval([strcat('C3_array(', nm1,') = mean(C3', n, '(100:end));')])
    hold on
    xlim([20 150])
    ylim([0 2])
    
    % Neg def
    subplot(2,2,4)
    eval([strcat('plot(C4', n, ')')])
    eval([strcat('C4_sum = C4_sum + C4', n, ';')])
    eval([strcat('C4_array(', nm1,') = mean(C4', n, '(100:end));')])
    hold on
    xlim([20 100])
    ylim([-2 0])     
    
    eval([strcat('vpbp_array(', nm1,') = mean(vpbp', n, '(100:end));')])
    
end
    C_avg = C_sum / count;
    C2_avg = C2_sum / count;
    C3_avg = C3_sum / count;
    C4_avg = C4_sum / count;    
    subplot(2,2,1)
    plot(C_avg, 'Color', 'k', 'LineWidth', 2)
    subplot(2,2,2)
    plot(C2_avg, 'Color', 'k', 'LineWidth', 2)
    subplot(2,2,3)
    plot(C3_avg, 'Color', 'k', 'LineWidth', 2)
    subplot(2,2,4)
    plot(C4_avg, 'Color', 'k', 'LineWidth', 2)
    
%save(strcat('C_array', output_dir, '.mat'), 'C_array')
%save(strcat('C2_array', output_dir, '.mat'), 'C2_array')
%save(strcat('C3_array', output_dir, '.mat'), 'C3_array')
%save(strcat('C4_array', output_dir, '.mat'), 'C4_array')

save(strcat('vpbp_array', output_dir, '.mat'), 'vpbp_array')
end

%%%%%%%% SLOPE BY WAVENUMBER %%%%%%%%%%%%%
if CALC_SLOPE_BY_WAVENUMBER == 1
    wavenumbers = [0:256 -255:1:-1];
   
for k = 1:layers-1
  n = num2str(k);
  eval([strcat('psidiff', n, ' = [];')])
  eval([strcat('vmean', n, ' = [];')])
  eval([strcat('qmean', n, ' = [];')])  
end

for k = 1:layers-1
 n = num2str(k);
 np1 = num2str(k+1);
 eval([strcat('psidiff', n, ' = psi', n, ' - psi', np1, ';')]);
 eval([strcat('vmean', n, ' = 0.5*(v', n, ' + v', np1, ');')]);
 eval([strcat('qmean', n, ' = 0.5*(q', n, ' + q', np1, ');')]);
end
'Done making b, v, and q.'

for k = 1:layers-1
 n = num2str(k);
 np1 = num2str(k+1);
 eval([strcat('clear psi', n, ';')]);
 eval([strcat('clear psie', n, ';')]);
 eval([strcat('clear u', n, ';')]);
 eval([strcat('clear v', n, ';')]);
 eval([strcat('clear q', n, ';')]);
end
'Cleared old variables.'

dt = diff(time);
for k = 1:layers-1
 n = num2str(k);
 for t = 2:length(time)  
 eval([strcat('dbdt', n, '(:,:,t) = (psidiff', n, '(:,:,t) - psidiff', n, '(:,:,t-1)) ./ dt(t-1);')]);
 eval([strcat('dqdt', n, '(:,:,t) = (qmean', n, '(:,:,t) - qmean', n, '(:,:,t-1)) ./ dt(t-1);')]);
 end
 eval([strcat('dbdt', n, '(:,:,1) = 0;')]);
 eval([strcat('dqdt', n, '(:,:,1) = 0;')]); 
end
'Done making dbdt.'

for k = 1:layers-1
 n = num2str(k)
 for t = 1:size(dbdt1,3)
     eval([strcat('b = psidiff', n, '(:,:,t);')])
     eval([strcat('dbdt = dbdt', n, '(:,:,t);')])
     eval([strcat('v = vmean', n, '(:,:,t);')])
     eval([strcat('q = qmean', n, '(:,:,t);')])
     eval([strcat('dqdt = dqdt', n, '(:,:,t);')])

     eval([strcat('num', n, '(:,:,t) = fft2(b) .* conj(fft2(dbdt)) + conj(fft2(b)) .* fft2(dbdt);')]);
     eval([strcat('denom', n, '(:,:,t) = fft2(b) .* conj(fft2(v)) + conj(fft2(b)) .* fft2(v);')]);
     eval([strcat('num2', n, '(:,:,t) = fft2(b) .* conj(fft2(v)) + conj(fft2(b)) .* fft2(v);')]);
     eval([strcat('denom2', n, '(:,:,t) = fft2(v) .* conj(fft2(q)) + conj(fft2(v)) .* fft2(q);')]);
     eval([strcat('num3', n, '(:,:,t) = fft2(q) .* conj(fft2(dbdt)) + conj(fft2(q)) .* fft2(dbdt) + conj(fft2(b)) .* fft2(dqdt)  + fft2(b) .* conj(fft2(dqdt));')]); 
%     eval([strcat('num3', n, '(:,:,t) = fft2(q) .* conj(fft2(b)) + conj(fft2(b)) .* fft2(dqdt);')]);      
 end
end
'Done making num and denom.'

% 11:   15:25

for k = 1:layers-1
 n = num2str(k)
 eval([strcat('tmp1 = squeeze(num', n, '(1,:,:));')])
 eval([strcat('nummean(', n, ',:) = squeeze( mean( tmp1(:,100:end), 2));')]);  % 100:end
 eval([strcat('tmp2 = squeeze(denom', n, '(1,:,:));')])
 eval([strcat('denommean(', n, ',:) = squeeze( mean( tmp2(:,100:end), 2));')]);
 eval([strcat('tmp = tmp1 ./ tmp2;')])
 eval([strcat('tmp3 = squeeze(num2', n, '(1,:,:));')])
 eval([strcat('nummean2(', n, ',:) = squeeze( mean( tmp3(:,100:end), 2));')]);  % 100:end
 eval([strcat('tmp4 = squeeze(denom2', n, '(1,:,:));')])
 eval([strcat('denommean2(', n, ',:) = squeeze( mean( tmp4(:,100:end), 2));')]);
 eval([strcat('tmp5 = squeeze(num3', n, '(1,:,:));')])
 eval([strcat('nummean3(', n, ',:) = squeeze( mean( tmp5(:,100:end), 2));')]);  % 100:end
 %eval([strcat('ratio_max(', n, ',:) = squeeze( max( tmp1(:,25:end) ./ tmp2(:,25:end), [], 2));')]);
 %eval([strcat('ratio_min(', n, ',:) = squeeze( min( tmp1(:,2:end) ./ tmp2(:,2:end), [], 2));')]);
end
'Done time-averaging.'

 pos_finder = zeros(1,512);
% finder = find(nummean(10,:) > 0.000001 * max(nummean(10,:)));
% pos_finder(finder) = 1;

load(strcat('growth_rate_lim', output_dir, '.mat'))
pos_finder(lim) = 1;

figure
C_wavenumber = (1 - (1 ./ deltaU) * (nummean ./ denommean)) .* pos_finder;
%C_wavenumber = 1 - (1 ./ deltaU) * ratio_min;
plot(C_wavenumber(10,:))
ylim([-1 2])
xlim([0 50])

save(strcat('C_b_wavenumbers', output_dir, '.mat'), 'C_wavenumber')

% n = 20
% plot(1 - squeeze(num10(1,n,:)) ./ squeeze(deltaU*denom10(1,n,:)) )
%plot(squeeze(num10(1,n,:))  )
%plot(squeeze(deltaU*denom10(1,n,:)) )

term3 = nummean2 *beta1 ./ (denommean2 * (-deltaU));
term2 = nummean3 ./ (denommean2 * (-deltaU));
residual = C_wavenumber - 1 - term2 - term3;

save(strcat('term3_wavenumbers', output_dir, '.mat'), 'term3')
save(strcat('term2_wavenumbers', output_dir, '.mat'), 'term2')
save(strcat('residual_wavenumbers', output_dir, '.mat'), 'residual')


% plot(nummean2(10,:) *beta1 ./ (denommean2(10,:) * - deltaU) )
% hold on
% plot(nummean3(10,:) ./ (denommean2(10,:) * - deltaU) )
% xlim([0 30])
% ylim([-3 3])

end
