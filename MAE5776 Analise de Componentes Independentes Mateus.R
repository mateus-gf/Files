#=======================#
# Criptografia com ICA  #
#=======================#

#---------------------------------------------------------------------------------------
#
# É um problema de processamento de sinais
# Descriptografia de sinais não-gaussianos por meio de blind source separation
#
#
# Autor: Mateus Gonzalez de Freitas Pinto
# Instituto de Matemática e Estatística - Universidade de São Paulo
# Julho de 2022
#
#---------------------------------------------------------------------------------------

#====================#
# Pacotes utilizados #
#====================#

rm(list = ls())
cat("\014")

library(DescTools)
library(rootSolve)
library(stargazer)
library(wrassp)
library(tuneR)
library(audio)
library(ica)

#=========================#
#   Funções do programa   #
#=========================#

zero.crossings <- function(x){
  aprx <- approxfun(1:length(x), x)
  return(length(uniroot.all(aprx, interval = range(1:length(x)))))
}

relative.error <- function(signal, approximation){
  return(sqrt(sum((signal-approximation)^2))/sqrt(sum(signal^2)))
}

#========================#
# Parâmetros do programa #
#========================#

#Obs: na prática preciso de numeros veradadeiramente aleatórios. 
#Aqui só pego seed pra fins de reprodutibilidade da Key 
#Uma aplicação de criptografia precisa de fato de números verdadeiramente aleatórios

set.seed(1234)

beta <- 2

#======================#
#      Diretórios      #
#======================#

Slocal <- "D:\\audio_processing-master\\"

#Leitura dos arquivos

train_audio <-  tuneR::readWave(paste0(Slocal, "audio_file.wav"))

#O que o audio nos dá? 
str(train_audio)

#Comentários
#1. O aúdio dá 18593 pontos amostrais de 4000 taxa de amostragem. Isso dá 4.6 segundos = 18593/4000
#2. O áuido é de 16 bits ==> pressão do som é mapeada para números que vão de -2^15 a 2^15-1
#3. Como converter audio para float? 

x <- train_audio@left/2^(train_audio@bit-1)
length(x) /train_audio@samp.rate

#Represento como sinal plotando amplitude no tempo

timeArray <-  (0:(18593-1))/train_audio@samp.rate
som <- x[1:16384]
temp <- timeArray[1:16384]

#É um áudio de um coração batendo
plot(x = temp, y=som, main = "Aúdio - batidas de coração",
        col = 'darkblue', type='l', xlab='Tempo em segundos',
     ylab = "Amplitude")
grid()

#==========================#
#                          #
#  CRIPTOGRAFIA DOS DADOS  #
#                          #
#==========================#

#Separo em 2 frames
#E se não for par? Faz com padding ou qqr coisa do tipo 

frame1 <- som[1:8192]
frame2 <- som[8193:16384]

#Cada frame em 4 partições

s11    <- frame1[1:2048]
s12    <- frame1[2049:4096]
s13    <- frame1[4097:6144]
s14    <- frame1[6145:8192]

s21    <- frame2[1:2048]
s22    <- frame2[2049:4096]
s23    <- frame2[4097:6144]
s24    <- frame2[6145:8192]

#Chave conforme o artigo define
key <- runif(16384, min = -2, max = 2)

#Particiono a chave em 4 partes

key1 <- key[1:8192]
key2 <- key[8193:16384]

k11    <- key1[1:2048]
k12    <- key1[2049:4096]
k13    <- key1[4097:6144]
k14    <- key1[6145:8192]

k21    <- key2[1:2048]
k22    <- key2[2049:4096]
k23    <- key2[4097:6144]
k24    <- key2[6145:8192]

dev.off()

#================================#
#Levanto as waveform informations#
#================================#

#Faço interpolação linear pra contar os zeros

zc11 <- zero.crossings(scale(s11))
zc12 <- zero.crossings(scale(s12))
zc13 <- zero.crossings(scale(s13))
zc14 <- zero.crossings(scale(s14))

zc11;zc12;zc13;zc14

zc21 <- zero.crossings(scale(s21))
zc22 <- zero.crossings(scale(s22))
zc23 <- zero.crossings(scale(s23))
zc24 <- zero.crossings(scale(s24))

#SPMIN
spmin11 <- min(scale(s11))
spmin12 <- min(scale(s12))
spmin13 <- min(scale(s13))
spmin14 <- min(scale(s14))
spmin21 <- min(scale(s21))
spmin22 <- min(scale(s22))
spmin23 <- min(scale(s23))
spmin24 <- min(scale(s24))


#SPMAX

spmax11 <- max(scale(s11))
spmax12 <- max(scale(s12))
spmax13 <- max(scale(s13))
spmax14 <- max(scale(s14))
spmax21 <- max(scale(s21))
spmax22 <- max(scale(s22))
spmax23 <- max(scale(s23))
spmax24 <- max(scale(s24))

#Cosntruo como tabela

tab.1 <- rbind(c(zc11, spmin11, spmax11),
               c(zc12, spmin12, spmax12),
               c(zc13, spmin13, spmax13),
               c(zc14, spmin14, spmax14))

colnames(tab.1) <- c("NZC", "MIN", "MAX")
tab.1 <- as.data.frame(tab.1)

tab.2 <- rbind(c(zc21, spmin21, spmax21),
               c(zc22, spmin22, spmax22),
               c(zc23, spmin23, spmax23),
               c(zc24, spmin24, spmax24))

colnames(tab.2) <- c("NZC", "MIN", "MAX")
tab.2 <- as.data.frame(tab.2)

stargazer::stargazer(cbind(tab.1, tab.2),title = "Waveform - primeiro frame", summary = F, digits = 2, label = "waveform_frame1")

#Média e desvio padrão 

ms11 <- mean(s11) ; sds11 <- sd(s11)
ms12 <- mean(s12) ; sds12 <- sd(s12)
ms13 <- mean(s13) ; sds13 <- sd(s13)
ms14 <- mean(s14) ; sds14 <- sd(s14)
ms21 <- mean(s21) ; sds21 <- sd(s21)
ms22 <- mean(s22) ; sds22 <- sd(s22)
ms23 <- mean(s23) ; sds23 <- sd(s23)
ms24 <- mean(s24) ; sds24 <- sd(s24)

#Gero a matriz s2p a partir das chaves

s2p1 <- cbind(s11, s12, s13, s14, k11, k12, k13, k14)
s2p2 <- cbind(s21, s22, s23, s24, k21, k22, k23, k24)

#Gera a matriz de mixing

set.seed(1321)

B <- matrix(NA,nrow = 4, ncol = 4)

for(i in 1:4){
  B[i,] <- runif(4,-1,1)
}

Ae <- cbind(B, beta*B)

xp1 <- Ae %*% t(s2p1)
xp2 <- Ae %*% t(s2p2)

#Sinais originais

par(mfrow = c(2,4))

plot.ts(s11, axes = T, xy.labels = F)
plot.ts(s12, axes = T, xy.labels = F)
plot.ts(s13, axes = T, xy.labels = F)
plot.ts(s14, axes = T, xy.labels = F)
plot.ts(s21, axes = T, xy.labels = F)
plot.ts(s22, axes = T, xy.labels = F)
plot.ts(s23, axes = T, xy.labels = F)
plot.ts(s24, axes = T, xy.labels = F)
mtext("Sinais originais", side = 3, line = -2, outer = TRUE)

#Sinais criptografados

par(mfrow = c(2,4))

plot.ts(k11, axes = T, xy.labels = F)
plot.ts(k12, axes = T, xy.labels = F)
plot.ts(k13, axes = T, xy.labels = F)
plot.ts(k14, axes = T, xy.labels = F)
plot.ts(k21, axes = T, xy.labels = F)
plot.ts(k22, axes = T, xy.labels = F)
plot.ts(k23, axes = T, xy.labels = F)
plot.ts(k24, axes = T, xy.labels = F)
mtext("Sinais de chave", side = 3, line = -2, outer = TRUE)

dev.off()

par(mfrow = c(2,4))


for(i in 1:nrow(xp1)){
  plot.ts(xp1[i,], xy.labels = F, axes = T)
}

for(i in 1:nrow(xp2)){
  plot.ts(xp2[i,], xy.labels = F, axes = T)
}
mtext("Sinais criptografados", side = 3, line = -2, outer = TRUE)




dev.off()

#==========================#
#                          #
#     DESCRIPTOGRAFIA      #
#                          #
#==========================#


#Como funciona? 

#Na prática, eu preciso passar o sinal chave se for puramente aleatório pro receptor da mensagem 
#Se for pseudoaleatório, eu passo a semente. 
#Com isso, faço o append e fastICA pra obter o sinal de volta

vec1 <- (cbind(xp1[1,],xp1[2,], xp1[3,],xp1[4,], k11, k12, k13, k14))
vec2 <- (cbind(xp2[1,],xp2[2,], xp2[3,],xp2[4,], k21, k22, k23, k24))

v1 <- ica::icafast(vec1, nc = 8)
v2 <- ica::icafast(vec2, nc = 8)

plot.ts(v1$S, main="Sinal chave x Audio parte 1")
plot.ts(v2$S, main="Sinal chave x Audio parte 2")

#Vem embaralhado e fora de escala 

par(mfrow=c(2,4))
plot.ts(v1$S[,5])
plot.ts(v1$S[,6])
plot.ts(v1$S[,7])
plot.ts(v1$S[,8])
plot.ts(s11) 
plot.ts(s12)
plot.ts(s13)
plot.ts(s14)
mtext("Sinais fora de escala (cima) e originais (baixo)", side = 3, line = -2, outer = TRUE)
dev.off()

par(mfrow=c(2,4))
plot.ts(v2$S[,5])
plot.ts(v2$S[,6])
plot.ts(v2$S[,7])
plot.ts(v2$S[,8])
plot.ts(s21) 
plot.ts(s22)
plot.ts(s23)
plot.ts(s24)
mtext("Sinais fora de escala (cima) e originais (baixo)", side = 3, line = -2, outer = TRUE)

dev.off()

#Agora basta fazer o rescaling
#Depois calculo o signal to noise ratio

zcS15 <- zero.crossings(v1$S[,5])
zcS16 <- zero.crossings(v1$S[,6])
zcS17 <- zero.crossings(v1$S[,7])
zcS18 <- zero.crossings(v1$S[,8])

maxS15 <- max(v1$S[,5])
maxS16 <- max(v1$S[,6])
maxS17 <- max(v1$S[,7])
maxS18 <- max(v1$S[,8])

minS15 <- min(v1$S[,5])
minS16 <- min(v1$S[,6])
minS17 <- min(v1$S[,7])
minS18 <- min(v1$S[,8])


zcS25 <- zero.crossings(v2$S[,5])
zcS26 <- zero.crossings(v2$S[,6])
zcS27 <- zero.crossings(v2$S[,7])
zcS28 <- zero.crossings(v2$S[,8])

minS25 <- min(v2$S[,5])
minS26 <- min(v2$S[,6])
minS27 <- min(v2$S[,7])
minS28 <- min(v2$S[,8])

maxS25 <- max(v2$S[,5])
maxS26 <- max(v2$S[,6])
maxS27 <- max(v2$S[,7])
maxS28 <- max(v2$S[,8])


#Comparao NZC por ordenação

tabC1 <- cbind(c(5:8),
               c(zcS15,zcS16,zcS17,zcS18),
               c(minS15,minS16,minS17,minS18),
               c(maxS15,maxS16,maxS17,maxS18))


colnames(tabC1) <- c("SINAL","NZC", "MIN", "MAX")
tabC1 <- as.data.frame(tabC1)

tabC1
tab.1

comp1 <- cbind(tab.1, tabC1)

#Gero tabela em Latex para poder incluir no relatório
stargazer::stargazer(comp1, summary = F, label = "waveform1", digits = 2)

#Ordem final: 6, 5, 8,  7

#Faço o mesmo pro frame 2

tabC2 <- cbind(c(5:8),
               c(zcS25,zcS26,zcS27,zcS28),
               c(minS25,minS26,minS27,minS28),
               c(maxS25,maxS26,maxS27,maxS28))


colnames(tabC2) <- c("SINAL","NZC", "MIN", "MAX")
tabC2 <- as.data.frame(tabC2)

tabC2
tab.2

comp2 <- cbind(tab.2, tabC2)

stargazer::stargazer(comp2, summary = F, label = "tab:waveform2", digits = 2)


#Reconstrução

s11hat <- ms11 + (sds11*v1$S[,6])
s12hat <- ms12 + (sds12*v1$S[,5])
s13hat <- ms13 + (sds13*v1$S[,8])
s14hat <- ms14 + (sds14*v1$S[,7])

s1.hat <- c(s11hat, s12hat, s13hat, s14hat)
s1.true <- c(s11, s12,s13,s14)

#Faço pra s2

s21hat <- ms21 + (sds21*v2$S[,5])
s22hat <- ms22 + (sds22*v2$S[,6])
s23hat <- ms23 + (sds23*v2$S[,7])
s24hat <- ms24 + (sds24*v2$S[,8])

s2.hat <- c(s21hat, s22hat, s23hat, s24hat)
s2.true <- c(s21, s22,s23,s24)

#Comparação para o sinal inteiro 

descr <- c(s1.hat, s2.hat)

plot.ts(descr, main = "Comparação Sinal Original x Descriptografado", col='darkblue', 
        ylim = c(-1,1))
lines(som, col='red', lty = 'dashed')
legend(0, -0.6, c("Reconstruído", "Original"), col= c("blue","red"), lty=c(1,2)) 

dev.off()

relative.error(signal = som,approximation = descr)
Metrics::rmse(actual = som, predicted = descr)
Metrics::mae(som, descr)

#E se quiser ouvir o áudio novo? 
#Vou criar um arquivo .wav e salvar

#Parâmetros do áudio: taxa de amostragem, precisão em bits, 

sampling <-  4000
bit      <-  16
pcm      <-  TRUE

audio.desc  <- Wave(left = descr*(2^(bit-1)), bit=bit, pcm=pcm, samp.rate = sampling)
writeWave(object = audio.desc, filename = "D:\\audio descriptografado.wav")
