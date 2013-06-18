repmat = function(X,m,n){
##R equivalent of repmat (matlab)
if (is.matrix(X)) {
mx = dim(X)[1]
nx = dim(X)[2]
} else {
   mx = 1
   nx = length(X)
}
matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}
