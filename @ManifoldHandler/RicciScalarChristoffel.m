function S = RicciScalarChristoffel(ig, C, DC1)

    S = 0;
    for a = 1:2
        for b = 1:2
            t1 = sum(diag(squeeze(DC1(:,a,b,:))));
            t2 = sum(diag(squeeze(DC1(:,a,:,b))));
            t3 = 0;
            t4 = 0;
            for c = 1:2
                for d = 1:2
                    t3 = t3 + C(d,a,b)*C(c,c,d);
                    t4 = t4 + C(d,a,c)*C(c,b,d);
                end
            end
            S = S + ig(a,b) * (t1 - t2 + t3 - t4);
        end
    end
end
