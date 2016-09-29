function logodds=calculateLogodds(planes,center,c_intensity,c_diameter,xymaximavals,xydetdiameters,anisotropy,xycoverage)
%calculate per plane vs per maxima a diameter and delta xy value and cut
%off if these go over thresholds

sp=size(planes);

if(max(sp)>0)
    logodds=zeros(1,sp(1));
    [features]=calc_disk_feature_vector(planes,c_intensity,c_diameter,center,xymaximavals,xydetdiameters,anisotropy)   ;

    global allbadlm
    global allgoodlm
    global allbadrm
    global allgoodrm

    global allbadlc
    global allbadrc

    global allgoodrc
    global allgoodlc

    pinlall=mvnpdf(features',allgoodlm,allgoodlc);
    polall=mvnpdf(features',allbadlm,allbadlc);
    pinrall=mvnpdf(features',allgoodrm,allgoodrc);
    porall=mvnpdf(features',allbadrm,allbadrc);

    pivot=center(3);
    topplane=min(min(planes(:,3)),pivot);
    bottomplane=max(max(planes(:,3)),pivot);
    counter=1;

    %loop ignores center planes leaving them to default of valid and
    %logodds zero
    for j=topplane:bottomplane-1

        leftdiskindexs=find(planes(:,3)==j);
        rightdiskindexs=find(planes(:,3)==j+1);
        leftdisks=planes(leftdiskindexs,:);
        rightdisks=planes(rightdiskindexs,:);
        sizes=size(leftdisks);
        sizes2=size(rightdisks);

        if(j<pivot)
            for h=1:sizes(1)
                pi=pinlall(counter);
                po=polall(counter);
                if(xycoverage(leftdisks(h,4))>13)
                    logodds(leftdiskindexs(h))=log(pi/po);
                end
                %else logodds is assigned value0
                counter=counter+1;
            end

        end
        if(j>=pivot)
            for h=1:sizes2(1)
                pi=pinrall(counter);
                po=porall(counter);
                if(xycoverage(rightdisks(h,4))>13)
                    logodds(rightdiskindexs(h))=log(pi/po);
                end
                counter=counter+1;
            end
        end
    end

end












