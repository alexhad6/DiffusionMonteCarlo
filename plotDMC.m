function plotDMC(p0,R,bins,limits,guided,plotfunc)
%PLOTDMC Plots the DMC distribution against the true wavefunction p0.
if guided
    h = histogram(R,bins,'Visible','off','Normalization','pdf');
    centers = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
    scaledValues = h.Values;
else
    h = histogram(R,bins,'Visible','off','Normalization','pdf');
    centers = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
    scaledValues = h.Values./sqrt(h.BinWidth*sum(h.Values.^2));
end
xvals = (limits(1):range(limits)/1000:limits(2));
hold on
p1 = bar(centers,scaledValues,1,'FaceColor',[0.4,0.6,0.8]);
if plotfunc
    p2 = plot(xvals,arrayfun(p0, xvals),'LineWidth',2);
    legend([p2, p1],'Actual','DMC')
else
    legend(p1,'DMC')
end
hold off
xlim(limits)
xlabel('Distance (bohr)')
if guided
    ylabel('Probability Density')
else
    ylabel('Wave Function')
end
end
