function Cijk = Christoffel(ig, D1, i, j, k)

    Cijk = ig(i,:)*(D1(:,j,k)+D1(:,k,j)-squeeze(D1(j,k,:)))/2.0;
    
end