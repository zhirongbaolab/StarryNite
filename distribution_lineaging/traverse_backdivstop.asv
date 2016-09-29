function [depth ] = traverse_backdivstop( esequence,t,i )
%traverses esequence sucessor links (always first sucessor if multiple) and
%returns depth before dead end or division
if(t==1||esequence{t}.pred(i)==-1||esequence{esequence{t}.pred_time(i)}.suc(esequence{t}.pred(i),2)~=-1)
    depth=1;
    return;
else
    depth=traverse_backdivstop(esequence,esequence{t}.pred_time(i),esequence{t}.pred(i));
    depth=depth+1;
end
end

