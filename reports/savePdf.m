function savePdf(h,name)
    set(gca(), 'LooseInset', get(gca(), 'TightInset'));
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','manual','PaperUnits','Inches', ...
        'PaperPosition',[0.02*pos(3),0.02*pos(4),1.02*pos(3),1.02*pos(4)], ...
        'PaperSize',1.05*[pos(3), pos(4)])
    print(h,name,'-dpdf','-r0')
end