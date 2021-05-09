function ID = STAR_NAME(Data, Pair,ac_er_ang, ac_er_dis,ID_known)
    SizeD = size(Data);
    m=0;
    for i = ID_known+1:(SizeD(1)-1)
        for j = i+1:SizeD(1)
            m = m+1;
            K = STAR_ID(Data,Pair,ac_er_ang, ac_er_dis,ID_known+1,i,j);
            SK= size(K);
            TriStar(1:SK(1),m,:) = K(:,:);
        end
    end

    S = size(TriStar);
    k = 0;
    PivotStar = [0,0];
    for i = 1:S(2)
        PivotStar(:,3) = 1;
        for j = 1:S(1)
            k = k + 1;
            con    = double(sum(PivotStar(:,1)==TriStar(j,i))==0)*2-1;
            con1   = double((TriStar(j,i)==0)&(con<0));
            PivotStar(k,1) = (con)*((TriStar(j,i))+con1);
            pos    = double(uint8(PivotStar(:,1)==TriStar(j,i)));
            PivotStar(:,2) = pos.*PivotStar(:,3)+PivotStar(:,2);
            PivotStar(:,3) = 1-pos;
        end
    end
    PivotStar = PivotStar(PivotStar(:,1)>0,:);
    [~,index]= sort(PivotStar(:,2),'descend');
    PivotStar = PivotStar(index,:);
    if length(PivotStar)<1
        ID = [-99999;-99999];
    else
        ID = PivotStar(1,1:2);
    end
end