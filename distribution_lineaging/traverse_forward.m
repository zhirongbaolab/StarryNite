function [depth ] = traverse_foward( esequence,t,i )
%traverses esequence sucessor links (always first sucessor if multiple) and returns depth before dead end
if(esequence{t}.suc(i,1)==-1)
    depth=1;
    return;
else
    depth=traverse_forward(esequence,esequence{t}.suc_time(i,1),esequence{t}.suc(i,1));
    depth=depth+1;
end
end

