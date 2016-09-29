function [depth ] = traverse_back( esequence,t,i )
%traverses esequence sucessor links (always first sucessor if multiple) and returns depth before dead end
if(t==1||esequence{t}.pred(i)==-1)
    depth=1;
    return;
else
    depth=traverse_back(esequence,esequence{t}.pred_time(i),esequence{t}.pred(i));
    depth=depth+1;
end
end

