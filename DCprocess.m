function [ph1,newDCR,newDCI,freCount] = DCprocess(x,maxR,minR,preDCR,maxI,minI,preDCI,freCount)
%process the DC with the angle
size2 = length(x);
tempval = sum(abs(x));
newDCR = preDCR;
newDCI = preDCI;
if(tempval/120 > 150)
    ph = angle(x);
    ph1 = unwrap(ph);
    for i = 2:size2
        if abs(ph1(i)-ph1(i-1)) > pi/4
            newDCR = (1 - 0.5)*preDCR + (maxR + minR)/2*0.5;
            newDCI = (1 - 0.5)*preDCI + (maxI + minI)/2*0.5;
        
           
        end

    end
    freCount = 1;
else
    ph1 = angle(x);
end
end