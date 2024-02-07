function graphics_setup(plot_type)

switch plot_type
case '1by1'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 9, 8]);
set(gca,'fontsize',8)

case '1by1.5'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 17, 9.5]);
set(gca,'fontsize',12)

case '1by2'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 18, 8]);
set(gca,'fontsize',8)

case '1by3'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 24, 7]);
set(gca,'fontsize',8)

case '1by3.5'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 24, 7.5]);
set(gca,'fontsize',8)

case '1by3.5_2'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 20, 6.5]);
set(gca,'fontsize',8)

case '2by2'
set(gcf, 'units', 'centimeters', 'Position',  [3, 3, 18, 14]);
set(gca,'fontsize',8)

case '3by2'
set(gcf, 'units', 'centimeters', 'Position',  [0, -2, 18, 22]);
set(gca,'fontsize',8)
end

return