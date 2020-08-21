clear all

path(pathdef)

deltaU = 0.0013;
U = 0.025;
Ld = 9e4;
nx = 512;
dx = Ld / 8; %1/(gridratio / Ld);
L = nx * dx;

wavelengths =  L ./ [Inf 1:256 -255:1:-1];
k = 2*pi ./ wavelengths;
mu = k * Ld;


hfig1 = figure(1);
set(hfig1, 'Position', [100 100 200 600])
set(hfig1, 'Color',[1 1 1]);
ylim([0 2])
xlim([1 17])
hold on
box on
grid on
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])


hfig2 = figure(2);
set(hfig2, 'Position', [100 100 200 600])
set(hfig2, 'Color',[1 1 1]);
ylim([0 2])
xlim([1 17])
hold on
box on 
grid on
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

hfig3 = figure(3);
set(hfig3, 'Position', [100 100 200 600])
set(hfig3, 'Color',[1 1 1]);
ylim([-2 2])
xlim([1 17])
hold on
box on 
grid on
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

C = linspecer(5)
count = 0;
for outdir = [21:25]%:25]%:26]%:14]
  count = count + 1;  
    
  load(strcat('../C_array', num2str(outdir), '.mat'))
  load(strcat('../C2_array', num2str(outdir), '.mat'))
  load(strcat('../C3_array', num2str(outdir), '.mat'))
  load(strcat('../C4_array', num2str(outdir), '.mat'))  
  
  figure(1)
  g = plot(C_array, 'LineWidth', 2, 'Color', C(count,:));
  clr = get(g, 'Color');
  %plot([1], mean(C_array), 'o', 'Color', clr, 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
  
  
  figure(2)
  if outdir ~= 21
    g2 = plot(C2_array, 'Color', clr, 'LineWidth', 2)
    %plot([1], mean(C2_array), 'o', 'Color', clr, 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
  end
  
  
  residual = C_array - 1 - C3_array - C4_array;
  figure(3)
  if outdir ~= 21
    % C3 is the pos def term
    g3 = plot(C3_array, 'Color', clr, 'LineWidth', 2)
    %plot([1], mean(C3_array), 'o', 'Color', clr, 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
    
    % C4 is the d/dt term
    g4 = plot(C4_array, 'Color', clr, 'LineWidth', 2, 'LineStyle', '--')
    %plot([1], mean(C4_array), 's', 'Color', clr, 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
    
    g5 = plot(residual, 'Color', clr, 'LineWidth', 2, 'LineStyle', ':')
    %plot([1], mean(residual), 'd', 'Color', clr, 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerSize', 8)

  end
  
end


count = 0;
for outdir = [21:25]%:25]%:26]%:14]
  count = count + 1; 

  load(strcat('../C_array', num2str(outdir), '.mat'))
  load(strcat('../C2_array', num2str(outdir), '.mat'))
  load(strcat('../C3_array', num2str(outdir), '.mat'))
  load(strcat('../C4_array', num2str(outdir), '.mat'))  
  
  figure(1)
  g = plot([1], mean(C_array), 'o', 'Color', C(count,:), 'MarkerFaceColor', C(count,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
  clr = get(g, 'Color');
  
  figure(2)
  if outdir ~= 21
    plot([1], mean(C2_array), 'o', 'Color', clr, 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
  end
  
  
  residual = C_array - 1 - C3_array - C4_array;
  figure(3)
  if outdir ~= 21
    % C3 is the pos def term
    plot([1], mean(C3_array), 'o', 'Color', clr, 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
    
    % C4 is the d/dt term
    plot([1], mean(C4_array), 'o', 'Color', clr, 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerSize', 8)
    
    plot([1], mean(residual), 'o', 'Color', clr, 'MarkerFaceColor', clr, 'MarkerEdgeColor', 'k', 'MarkerSize', 8)

  end
  
end  

% figure(1)
% export_fig C_b.eps -q101 -m2.5
% 
% figure(2)
% export_fig C_q.eps -q101 -m2.5
% 
% figure(3)
% export_fig C_b_terms.eps -q101 -m2.5
