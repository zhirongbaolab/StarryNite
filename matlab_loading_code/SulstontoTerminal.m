function [ terminal ] = SulstontoTerminal(sulston,partlist,deathlist )
%UNTITLED2 Summary of this function goes here
%Detailed explanation goes here
terminal={};
for j=1:length(sulston)
    sulstonj=sulston{j};
    if(isempty(sulstonj))
        terminal{j}=[]; %final answer
    else
    terminal{j}=[]; %ensures some answer in first round blank if no match
    for i=1:length(partlist)
        if(strcmp(lower(sulstonj),lower(partlist{i,2})))
            terminal{j}=partlist{i,1};
         
        end
    end
    if(isempty(terminal{j}))
        %is death or non terminal
        for i=1:length(deathlist)
            if(strcmp(lower(sulstonj),lower(deathlist{i})))
                terminal{j}='Cell Death';
            
            end
        end
        if isempty(terminal{j})
           terminal{j}='Non-terminal Cell';
        end
        
    end
    
    end
end
terminal=terminal';

