function [H,status] = nls_h(x,y)
 H(1) = 2 * x(3) * y(1) ;
 H(2) = 0 ;
 H(3) = 2 * y(2) ;
 H(4) = 2 * x(1) * y(1) ;
 H(5) = 0 ;
 H(6) = 0 ;
 status = 0 ;
end
