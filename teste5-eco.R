#teste do IV

install.packages("AER")
install.packages("stargazer")
install.packages("robustbase")
install.packages("gmm")

library("AER")
library("stargazer")
library("robustbase")
library("gmm")

data("CigarettesSW")
summary(CigarettesSW)
data95 <- subset(CigarettesSW, year == "1995")
data95

#nossas variáveis 

y <- data95$packs #  quantidades de cigarro vendido
x <- (data95$price / data95$cpi) #  preço de cada maço
z1 <- ((data95$taxs - data95$tax) / data95$cpi ) # Variavel instrumental o
#imposto de venda
z2 <-  data95$income <- with(data95, income / population / cpi) #renda
z3 <- data95$cigtax <- with(data95, tax/cpi) #taxa sobre o cigarro

#criando matrizes

I <- matrix(rep(1,48),48,1) #Matriz 1x48 de 1
X <- as.matrix(cbind(I,x,z2)) #Matriz X
Y <- as.matrix(cbind(y)) #Matriz com a variÃ¡veil Y
Z <- as.matrix(cbind(I,z1,z2,z3)) #Matriz com a variÃ¡vel Z


#mqo para o erro
Transposta_X = t(X) #transposta de X
Transposta_Z = t(Z) #tranposta de Z

#Calcular erro do MQO
Inversa_Transposta_X_vezes_X = solve(Transposta_X %*% X)
P = X %*% Inversa_Transposta_X_vezes_X  %*% Transposta_X
M = diag(48) - P
e = M%*%Y
e_quadro = e^2


# calcule D
d = matrix(e_quadro, 48, 48)
D =  diag(diag(d))

#para W heterocedastico, temos o S
s_heter = (Transposta_Z %*% D %*% Z)/48

#w de hetero

W_heter = solve(s_heter)

#beta hetero 

beta_gmm_het = (solve(t(X) %*% Z %*% W_heter %*% t(Z) %*% X) %*% (t(X) %*% Z %*% W_heter %*% t(Z) %*% Y))


#para conferir

gmm(y ~ x + z2,Z)



#para W com homoced, temos que
sigmma_2 = (t(e)%*%e)/45
s_homo = as.numeric(sigmma_2)* ((t(Z) %*% Z)/48)
w_homo = solve(s_homo)

#beta homo
beta_gmm_homo = (solve(t(X) %*% Z %*% w_homo %*% t(Z) %*% X) %*% (t(X) %*% Z %*% w_homo %*% t(Z) %*% Y))


#para conferir
reg_iv <- ivreg(Y ~ X | Z)




