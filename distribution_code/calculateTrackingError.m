%calculate tracking error by testing for every matched correct nucleus at
%t+1  is its predecessor(if matched) matched to the predecessor of its
%match


function [correct,wrong]=calculateTrackingError(p1_sucessors,p1_test_sucessors,match1,match1r, match2, match2r)


correct=0;
wrong=0;
s=size(p1_sucessors);

%should make ordering of children not matter

for i=1:s(1)
    if(match1r(i)~=-1) %real nucleus at t was matched to detected result (otherwise det FN)
        if(p1_sucessors(i,1)~=-1&&p1_test_sucessors(match1r(i),1)~=-1) %has sucessor 1 (otherwise death or sucessor det FN)
            if(p1_sucessors(i,1)==match2(p1_test_sucessors(match1r(i),1)))
                correct=correct+1;
            else
                if (p1_test_sucessors(match1r(i),2)~=-1&&p1_sucessors(i,1)==match2(p1_test_sucessors(match1r(i),2)))
                    correct=correct+1;
                else
                    wrong=wrong+1;
                end
            end
        end
        %ditto test for sucessor two
        if(p1_sucessors(i,2)~=-1&&p1_test_sucessors(match1r(i),2)~=-1) %has sucessor 1 (otherwise death or sucessor det FN)
            if(p1_sucessors(i,2)==match2(p1_test_sucessors(match1r(i),2)))
                correct=correct+1;
            else
                if(p1_test_sucessors(match1r(i),1)~=-1&&p1_sucessors(i,2)==match2(p1_test_sucessors(match1r(i),1)))
                      correct=correct+1;
                else
                wrong=wrong+1;
                end
            end
        end
    end
end


