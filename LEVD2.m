function [outputR,outputI,newDCR,newDCI,LmaxR,LminR,LmaxI,LminI] = LEVD2(inputR,inputI,preDCR,preDCI,premaxR,preminR,premaxI,preminI)
%对输入的矩阵进行LEVD处理
size1 = length(inputR);
LmaxR = max(inputR);
LminR = min(inputR);
LmaxI = max(inputI);
LminI = min(inputI);


%get real part variance
temp = inputR - ones(size1,1)*inputR(1);
dsum = sum(temp)/size1;
vsum = sum(temp.^2)/size1;

temp = inputI - ones(size1,1)*inputI(1);
dsum = dsum + sum(temp)/size1;
vsum = vsum + sum(temp.^2)/size1;

newDCR = preDCR;
newDCI = preDCI;

if vsum + dsum * dsum > 150  %15000 is power threshold
    if (LmaxR > premaxR || (LmaxR > preminR + 22 && (LmaxR-LminR) > 22 * 4)) %amplitude threshold
        
    else
        LmaxR = premaxR;
    end
    
    if (LminR < preminR || (LminR < premaxR - 22 && (LmaxR-LminR) > 22 * 4))
        
    else
        LminR = preminR;
    end
    
    if (LmaxI > premaxI || (LmaxI > preminI + 22 && (LmaxI-LminI) > 22 * 4)) %amplitude threshold
        
    else
        LmaxI = premaxI;
    end
    
    if (LminI < preminI || (LminI < premaxI - 22 && (LmaxI-LminI) > 22 * 4))
        
    else
        LminI = preminI;
    end
    
    if ((LmaxR - LminR) > 22 && (LmaxI - LminI) > 22)
        newDCR = (1-0.25)*preDCR + (LmaxR + LminR)/2*0.25;
        newDCI = (1-0.25)*preDCI + (LmaxI + LminI)/2*0.25;
    end
    

end

outputR = ones(size1,1)*newDCR;
outputI = ones(size1,1)*newDCI;

end