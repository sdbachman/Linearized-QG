%clear all
 
path(pathdef)


layer = 10;

% Growth rates at the "kink"
% 0.0    2.40e-8
% 0.1    2.09e-8
% 0.2    1.92e-8
% 0.3    1.60e-8
% 0.4    1.11e-8

hfig = figure;
set(hfig, 'Position', [100 100 1000 600])
set(hfig, 'Color',[1 1 1]);
hold on

for outdir = [21]%:14]
directory = '/project/oce/bachman/MITgcm/dedalus/two_layer_beta/';
output_dir = num2str(outdir)
files = dir(fullfile(directory, strcat(output_dir, '_1*')));
cd(directory)

deltaU = 0.0013;
U = 0.025;
Ld = 9e4;
nx = 512;
dx = Ld / 8; %1/(gridratio / Ld);
L = nx * dx;

wavelengths =  L ./ [Inf 1:256 -255:1:-1];
k = 2*pi ./ wavelengths;
mu = k * Ld;


snaps = strcat(directory, files(1).name, strcat('/snapshots/snapshots.h5'));

%h5info(snaps, '/') 
time = h5read(snaps,'/scales/sim_time');
y = h5read(snaps,'/scales/y/1.0');
x = h5read(snaps,'/scales/x/1.0');

for k = layer
  n = num2str(k);
  eval([strcat('psi', n, ' = [];')])
end

for k = layer
  n = num2str(k);
  eval([strcat('psi', n, ' = cat(3, psi', n, ', squeeze(h5read(snaps, ''/tasks/psi', n, ''')));')]);
end

  for k = layer
    n = num2str(k);
    eval([strcat('psi', n, ' = cat(1, psi', n, '(end,:,:), psi', n, ', psi', n, '(1,:,:));')])
    eval([strcat('psi', n, ' = cat(2, psi', n, '(:,end,:), psi', n, ', psi', n, '(:,1,:));')])
  end
  'Finished making psi.'
 
  for k = layer
    n = num2str(k);
    eval([strcat('v', n, '= 0.5*(psi', n, '(2:end-1,3:end,:) - psi', n, '(2:end-1,1:end-2,:)) / dx;')])
    eval([strcat('u', n, '= -0.5*(psi', n, '(3:end,2:end-1,:) - psi', n, '(1:end-2,2:end-1,:)) / dx;')])
  end    
  'Finished making u and v.'
  

for k = layer
  for t = 1:length(time)
    n = num2str(k);
    eval([strcat('u_hat = fft2(u', n, '(:,:,t));')])
    u_hats = conj(u_hat);
    eval([strcat('v_hat = fft2(v', n, '(:,:,t));')])
    v_hats = conj(v_hat);
    KE_hat(:,:,t) = 0.5 * (u_hat .* u_hats + v_hat .* v_hats);
  end
end


dt = diff(time);


for k = 1 %:20%1:512
  for t = 2:size(dt)
      %growth_rate2(k,t,:) = (1 ./ (2*dt(t))) .* log10(KE_hat(k,:,t+1) ./ KE_hat(k,:,t));

      growth_rate(t,:) = (1 ./ (2*dt(t))) .* log10(KE_hat(k,:,t+1) ./ KE_hat(k,:,t));
  end
end  
     % growth_rate = squeeze(mean(growth_rate2,1));

% 
% for l = 20:50%1:size(growth_rate,2)
%     plot(growth_rate(:,l)) 
%     pause(0.5)
%     title(num2str(l))
%     gt(l) = getframe(gcf);
% end
% movie(gt)

% 11 30

var = std(growth_rate(50:end,:),1);
mn = mean(growth_rate(50:end,:),1);

%11  1e-10
%12  1e-11
%13  2e-12
gr = zeros(1,512)
if str2num(output_dir) == 21
    lim = [1:26];
elseif str2num(output_dir) == 22
    lim = [10:28];
elseif str2num(output_dir) == 23
    lim = [13:33];
elseif str2num(output_dir) == 24
    lim = [16:38];  %50
elseif str2num(output_dir) == 25
    lim = [19:50];  %50     
end
gr(lim) = 1;

%threshold = 1e-9

gr = gr .* mn;
gr(lim(1)) = 0;

%plot(mu(lim), 2*gr(lim) * Ld / U, 'LineWidth', 2)

% plot([1:512], 2*gr * Ld / U)
% xlim([1 50])

save(strcat('growth_rate', num2str(outdir), '.mat'), 'gr')
save(strcat('growth_rate_lim', num2str(outdir), '.mat'), 'lim')


end

% xlim([0 4])
% set(gca, 'FontSize', 12)
