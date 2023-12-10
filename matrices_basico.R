A.1 <- c(0,1,1,0)
A.2 <- c(1,0,1,1)
A.3 <- c(1,1,0,0)
A.4 <- c(0,1,0,0)
A <- rbind(A.1,A.2,A.3,A.4)
At <- t(A)
M.1 <- c(1,1,0,0)
M.2 <- c(0,1,1,1)
M.3 <- c(1,0,1,0)
M.4 <- c(0,0,0,1)
M <- rbind(M.1,M.2,M.3,M.4)
Mt <- t(M)
M*Mt
Q <- A
degs <- rowSums(A)
for (i in (1:length(degs)))
  Q[i,i] <- degs[i]
print(Q)
print(M%*%Mt)
N.1 <- c(-1,-1,0,0)
N.2 <- c(0,1,1,-1)
N.3 <- c(1,0,-1,0)
N.4 <- c(0,0,0,1)
N <- rbind(N.1,N.2,N.3,N.4)
Nt <- t(N)
L <- - A
for (i in (1:length(degs)))
  L[i,i] <- degs[i]
print(L)
print(N%*%Nt)
u <- c(1,2,1,4)
print("u")
print(u)
print("A")
print(A)
# The adjacency matrix can be viewed as an operator g = Au = SUM fi sobre eij
g <- A%*%u    # Para cada fila, el resultado es la suma de los valores de u que corresponden a un enlace
print("Au")
print(g)
# It can also be viewed as a quadratic form u'Au = SUM fifj sobre (eij)
h <- u%*%A%*%u
print("uAu'")
print(h)

print("u")
print(u)
print("L")
print(L)
# La laplaciana es un operador Lu = 0.5 SUM (fi-fj)^2 sobre eij
ol <- L%*%u
print("L*u")
print(ol)
# Y su forma cuadrática es
fc <- u%*%L%*%u
print("u*L*u'")
print(fc)