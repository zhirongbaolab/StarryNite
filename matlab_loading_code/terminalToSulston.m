function [ sulston ] = terminalToSulston(terminal,partlist )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sulston=[];
for i=1:length(partlist)
    if(strcmp(lower(terminal),lower(partlist{i,1})))
        sulston=partlist{i,2};
    end
end

end

