function jd = julday(y,m,d,h)
  if m <= 2, y = y-1; m = m+12; end
     jd = floor(365.25*(y+4716))+floor(30.6001*(m+1))+d+h/24-1537.5; 
end
   
 

