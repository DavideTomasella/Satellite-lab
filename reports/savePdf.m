function savePdf(h,name,nomod)
    if nargin<3
        nomod = false;
    end
    set(gca(), 'LooseInset', get(gca(), 'TightInset'));
    set(h,'Units','Inches');
    pos = get(h,'Position');
    if nomod
        title("")
        set(h,'PaperPositionMode','manual','PaperUnits','Inches', ...
            'PaperPosition',[0.02*pos(3),0.02*pos(4),0.83*pos(3),1.02*pos(4)], ...
            'PaperSize',[0.95*pos(3), 1.05*pos(4)])
    else
        set(h,'PaperPositionMode','manual','PaperUnits','Inches', ...
            'PaperPosition',[0.02*pos(3),0.02*pos(4),1.02*pos(3),1.02*pos(4)], ...
            'PaperSize',1.05*[pos(3), pos(4)])
    end
    print(h,name,'-dpdf','-r0')
end