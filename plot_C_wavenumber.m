clear all
close all

path(pathdef)

deltaU = 0.0013;
U = 0.025;
Ld = 9e4;
nx = 512;
dx = Ld / 8; %1/(gridratio / Ld);
L = nx * dx;

ind = [1:512];
wavelengths =  L ./ [Inf 1:256 -255:1:-1];
k = 2*pi ./ wavelengths;
mu = k * Ld;


hfig1 = figure(1);
set(hfig1, 'Position', [100 100 800 300])
set(hfig1, 'Color',[1 1 1]);
ylim([0 1])
xlim([0 3.5])
hold on
box on
grid on
set(gca, 'YTick', [0 0.25 0.5 0.75 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

hfig1 = figure(2);
set(hfig1, 'Position', [100 100 800 300])
set(hfig1, 'Color',[1 1 1]);
ylim([0 1e-7])
xlim([0 3.5])
hold on
box on
grid on
set(gca, 'YTick', [0 0.25 0.5 0.75 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

hfig2 = figure(3);
set(hfig2, 'Position', [100 100 800 300])
set(hfig2, 'Color',[1 1 1]);
ylim([-2 1.5])
xlim([0 3.5])
hold on
box on 
grid on
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])



C = linspecer(5)
count = 0;
for outdir = [21:25]%:25]%:26]%:14]
  count = count + 1;  
   
  load(strcat('../C_b_wavenumbers', num2str(outdir), '.mat'))
  load(strcat('../growth_rate', num2str(outdir), '.mat'))
  load(strcat('../growth_rate_lim', num2str(outdir), '.mat'))

  if outdir ~= 21
  load(strcat('../term2_wavenumbers', num2str(outdir), '.mat'))
  load(strcat('../term3_wavenumbers', num2str(outdir), '.mat'))
  load(strcat('../residual_wavenumbers', num2str(outdir), '.mat')) 
  end
  
%   if outdir == 21
%       tmp = lim(end);
%       C_wavenumber(10,tmp) = 0;
%   end
  lim_ext = [lim(2:end-1)];
  Cn = (ind > 20) .* (C_wavenumber(10,:) < 0);
  if outdir > 22
  finder = find(Cn == 1);
  finder2 = finder(1);
  finder3 = find(lim_ext == finder2);
  lim_ext(finder3:end) = [];  
  end
  
  
  figure(1)
  plot(mu(lim_ext), C_wavenumber(10,lim_ext), 'LineWidth', 2, 'Color', C(count,:))
  g = plot([0], mean(C_wavenumber(10,lim_ext)), 'o', 'Color', C(count,:), 'MarkerFaceColor', C(count,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8)

  figure(2)
  tmp = C_wavenumber(10,lim_ext) .* gr(lim_ext);
  plot(mu(lim_ext), tmp, 'LineWidth', 2, 'Color', C(count,:))
  g = plot([0], mean(tmp), 'o', 'Color', C(count,:), 'MarkerFaceColor', C(count,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
  
  
  if outdir ~= 21
  figure(3)
  plot(mu(lim_ext), term2(10,lim_ext), 'LineWidth', 2, 'Color', C(count,:), 'LineStyle', '--')
  plot(mu(lim_ext), term3(10,lim_ext), 'LineWidth', 2, 'Color', C(count,:)) 
  plot(mu(lim_ext), residual(10,lim_ext), 'LineWidth', 2, 'Color', C(count,:), 'LineStyle', ':')
  end
end

figure(1)
export_fig C_wavenumber.eps -q101 -m2.5

figure(2)
ylim([0 max(tmp)])
set(gca, 'YTick', [0 .25*max(tmp) .5*max(tmp) .75*max(tmp) max(tmp)])
export_fig C_wavenumber_normalized.eps -q101 -m2.5

figure(3)
export_fig C_b_terms_wavenumber.eps -q101 -m2.5