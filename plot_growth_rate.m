clear all
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


hfig = figure;
set(hfig, 'Position', [100 100 800 300])
set(hfig, 'Color',[1 1 1]);
hold on

C = linspecer(5);
count = 0;
for outdir = [21:25]%:26]%:14]
  count = count + 1;
  
  load(strcat('../growth_rate', num2str(outdir), '.mat'))
  load(strcat('../growth_rate_lim', num2str(outdir), '.mat'))

  load(strcat('../C_b_wavenumbers', num2str(outdir), '.mat'))
  lim_ext = lim;% [lim(1:end-1)];

  if outdir == 21
      tmp = lim(end);
      gr(tmp) = 0;
  end
  Cn = (ind > 20) .* (C_wavenumber(10,:) < 0);
  if outdir > 22
  finder = find(Cn == 1);
  finder2 = finder(1);
  finder3 = find(lim_ext == finder2);
  lim_ext(finder3:end) = [];
  end
  
  plot(mu(lim_ext), 2*gr(lim_ext) * Ld / U, 'LineWidth', 2, 'Color', C(count,:))
 
end

xlim([0 3.5])
set(gca, 'FontSize', 12)
grid on
%set(gca, 'XTickLabel', [])
%set(gca, 'YTickLabel', [])
box on
%export_fig growth_rates.eps -q101 -m2.5
