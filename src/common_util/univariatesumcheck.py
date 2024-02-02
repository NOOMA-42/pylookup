

""" 
 Compute D(X ) = M (X, α) = Σmj=1 µj (α)ˆ
τcol(j)(X ) and find R(X ), Q2(X ) such that
D(X )t(X ) − φ(α) = XR(X ) + zI (X )Q2(X )
– Set ˆR = XN−m+2 – Output π2 = [D]2 = [D(x)]2, [R]1 = [R(x)]1, [ ˆR]1 = [ ˆR(x)]1, [Q2]1 = [Q2(x)]1.
"""

mu_j
rho_hat_col_j
#D(X )t(X ) − φ(α) = XR(X ) + zI (X )Q2(X )

Q2X, XRX = (Dx * tx - phi(alpha)).polydiv(zI(x))  