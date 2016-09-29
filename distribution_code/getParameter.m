function param=getParameter(name,numcells,location)
global parameters;

%find correct index into parameter arrays
if(numcells>min(parameters.staging))
    timeindex=max(find(parameters.staging<numcells))+1;
else
    timeindex=1;
end
%default values
if(length(parameters.(name))==1)
    param=parameters.(name);
else
    param=parameters.(name)(timeindex);
end

%if regions are defined for this time and value overwrite value
%if isempty(parameters.regions{timeindex})
if(isfield(parameters,'regions')&&length(parameters.regions)>=timeindex) %regions defined
    if(isfield(parameters.regions{timeindex},name))%regional def for this variable
        areadef=parameters.regions{timeindex}.area;

        if(location(1)>areadef(1)&&location(1)<=areadef(2)&&location(2)>areadef(3)&&location(2)<=areadef(4)&&location(3)>areadef(5)&&location(3)<=areadef(6))
            param=parameters.regions{timeindex}.(name);

        end
    end
end
