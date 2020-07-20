function plotDMC(psi0,R,bins,limits)
%PLOTDMC Summary of this function goes here
%   Detailed explanation goes here
h = histogram(R,bins,'Visible','off');
centers = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
scaledValues = h.Values./sqrt(h.BinWidth*sum(h.Values.^2));
xvals = (limits(1):range(limits)/1000:limits(2));
hold on
p1 = bar(centers,scaledValues,1,'FaceColor',[0.4,0.6,0.8]);
p2 = plot(xvals,arrayfun(psi0, xvals),'LineWidth',2);
hold off
xlim(limits)
legend([p2, p1],'Actual','DMC')
end
