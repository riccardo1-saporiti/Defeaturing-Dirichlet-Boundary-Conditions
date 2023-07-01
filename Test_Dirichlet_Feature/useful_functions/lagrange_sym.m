function [ y , dP_dx ] = lagrange_sym( X , V )

	syms x
	
	P = 0;
	
	for ii = 1 : length( X )
		PJ = V( ii );
		X_filt = X( X ~= X( ii ) );
		for j  = 1 : length( X ) - 1 
			PJ = PJ * ( x - X_filt( j ) ) / ( X( ii ) - X_filt( j ) ) ;
		end
		P = P + PJ ; 
	end
	
	y = matlabFunction( P ) ; 
    dP_dx = matlabFunction( diff( P , x ) );
    
return 