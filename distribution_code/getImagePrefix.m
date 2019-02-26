function [ prefix ] = getImagePrefix(imageName )
%matlab port of java code for image loading
% distinguish between _t###.TIF and diSPIM conventions
%cut out  just the number
slicepattern='(.+?)\d+.{4,5}$';
%result=(regexp(imageName,slicepattern,'tokens','forceCellOutput'));
%prefix=(result{1}{1}{1});
result=(regexp(imageName,slicepattern,'tokens'));
prefix=(result{1}{1});
%{
%partially ported this but gave up because its equivalent to shorter above
    String tID_8bitConvention = "-t";
     String tID_16bitConvention = "_t";
     String TIF_ext = ".TIF";
     String tif_ext = ".tif";
     String planeStr = "-p";

if (imageName.indexOf(tID_16bitConvention) ~= -1)
    prefix= imageName(0: imageName.lastIndexOf(tID_16bitConvention) + tID_16bitConvention.length());
else
    %the most we'll do to check if it's diSPIM is see if the characters between the last _ or - and the extension
    % are numbers
    
    % see if it's a fused diSPIM image or a single view (restrict the search to only the filename, not the full path)
     timeAppendCharacterType = '0';
    %// fused
    if ((imageName.substring(imageName.lastIndexOf(FORWARDSLASH))).lastIndexOf(UNDERSCORE) ~= -1)
        timeAppendCharacterType ='_';
    else
        if ((imageName.substring(imageName.lastIndexOf(FORWARDSLASH))).lastIndexOf('-') ~= -1)  // single view
            timeAppendCharacterType = '-';
        else
            ['Couldnt extract image prefix from: " + ',imageName,' \nUnable to find character type before time.")'];
            prefix='';
        end
    end
    
    %now that we know whether it's fused or single view, use that correct character type to extract the prefix
    String potentialTimeStr = imageName.substring(imageName.lastIndexOf(timeAppendCharacterType)+1, imageName.indexOf(tif_ext));
    
    for  i = 0:potentialTimeStr.length()
        if (~Character.isDigit(potentialTimeStr.charAt(i)))
             ['Couldnt extract image prefix from: " + ',imageName,' \nUnable to find character type before time.")'];
            prefix='';
        end
    end
    
    % // if we've reached here, we know that the potentialTimeStr is all digits, so we'll treat this as the time and everything
    %     // before it will be considered the image prefix
    %System.out.println("Found image prefix: " + imageName.substring(0, imageName.lastIndexOf(timeAppendCharacterType)+1));
    prefix=imageName.substring(0, imageName.lastIndexOf(timeAppendCharacterType)+1);
    
end
%}
end

