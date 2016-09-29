function [esequence] = processFPBifurcation( esequence,t,i,minsize,d1length )
%given that current bifurcation is FP fix it by deleting the shorter branch
%and making sure the longer is d1
d1cur=esequence{t}.suc(i,1);
d2cur=esequence{t}.suc(i,2);
tcur1=esequence{t}.suc_time(i,1);
tcur2=esequence{t}.suc_time(i,2);
if(d1length==minsize)%delete 1
    esequence{tcur1}.pred(d1cur)=-1;%unlink pred side of condemmed
     esequence{tcur1}.pred_time(d1cur)=-1;%unlink pred side of condemmed
    for j=1:minsize %mark short branch for deletion
        esequence{tcur1}.delete(d1cur)=1;
        d1curbk=d1cur;
        d1cur=esequence{tcur1}.suc(d1cur,1);
        tcur1=esequence{tcur1}.suc_time(d1curbk,1);
    end
    esequence{t}.suc_time(i,1)=esequence{t}.suc_time(i,2);
    esequence{t}.suc(i,1)=esequence{t}.suc(i,2);
else %delete 2
    esequence{tcur2}.pred(d2cur)=-1;
    esequence{tcur2}.pred_time(d2cur)=-1;
    for j=1:minsize %mark short branch for deletion
        esequence{tcur2}.delete(d2cur)=1;
        d2curbk=d2cur;
        d2cur=esequence{tcur2}.suc(d2cur,1);
        tcur2=esequence{tcur2}.suc_time(d2curbk,1);
    end
end
esequence{t}.suc_time(i,2)=-1;%unlink suc side of condemmed
esequence{t}.suc(i,2)=-1;

end

