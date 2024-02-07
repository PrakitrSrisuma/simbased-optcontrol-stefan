function export_figures(figure,filename)

filename1 = fullfile('Figures',  [filename,'.jpg']);
exportgraphics(figure, filename1,'Resolution',1000)
filename2 = fullfile('Figures',  [filename,'.pdf']);
exportgraphics(figure, filename2,'Resolution',1000)

return