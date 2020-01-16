#R tiene un sin fin de usos, por ejm para hacer cálculos simples
#MATEMATICA BASICA
2+3+sqrt(4) #sqrt es para hallar la raíz cuadrada
2*6 #para multiplicar
2^3 #para el exponencial
#MATEMATICA BASICA II
x <-c(1,2,3) #Con esto creo un vector 
y <-c(12,9,8)
plot(x,y) #para plotear los puntos de los vectores x e y
#OTRA FORMA DE ESTUDIAR VECTORES
mivector <- c(8, 6, 9, 10, 11, 12, 1) 
mivector #para ver el contenido de las variables de "mivector"
mivector [4] #esto es para ver que número hay en la posición 4 de "mivector" 
#USANDO LA FUNCION LISTAR
nuevalista <- list(name="proteína", nombre="DNA", mivector)
nuevalista
nuevalista[[3]]
#También se puede nombrar 
nuevalista$nombre #esto me dará el valor de la función nombre
#Introducción al análisis de secuencias
library("seqinr") #llamemos a la librería seqinr
M7 <- read.fasta(file = "M7.fasta") #para leer el fasta llamado M7 previamente guardado en mi PC
M7seq <- M7[[1]] #Se puede acceder a un elemento de una lista de objetos con doble corchete
#de aquí que guardamos la secuencia en una nueva variable llamada M7seq, es decir pone la secuencia en un vector
#Si mi interés es la región que va desde la base 20 hasta 100
M7seq[20:100]
M7seq
#Para hallar el largo de la secuencia
length(M7seq)
#si quiero saber la composición de las bases nitrogenadas de la secuencia
table(M7seq)
#Otra forma de contar
count(M7seq, 1) #hace lo mismo que table, al poner 1 es porque lo cuenta de 1 en 1
count(M7seq, 2)
count(M7seq, 3)
M7table <- count(M7seq,1)
table(M7seq)
M7table[[3]] #como ordena en alfabético serían las guaninas 
M7table[["a"]] #para encontrar por elemento (por eso entre comillas) en este caso adenina 
#contenido de GC
(1254+1135)*100/(5946) #de una forma manual
GC(M7seq) #calcula el contenido de GC de toda la secuencia
GC(M7seq)*100
GC(M7seq[20:100])*100
#para un primer
GC(M7seq[100:120])*100
#Para contar el contenido de GC por cada 1000 bases
starts <- seq(1, length(M7seq)-1000, by = 1000) #esto lo guarda en el vector starts
starts #Esto es para ver desde donde empezará a contar el GC
n <- length(starts) #me crea un valor n con 5 largos
for (i in 1:n) {
  trozar <- M7seq[starts[i]:(starts[i]+999)]
  trozarGC <- GC(trozar)
  print (trozarGC)
}
#hago este loop ordenando que separe en trozos a las secuencias

#Ahora lo hago cada 500 bases
starts <- seq(1, length(M7seq)-500, by = 500)
n <- length(starts) # Como en la anterior para encontrar el largo 
trozarGCs <- numeric(n) # Esto hace un vector del mismo largo al vector "starts" 
for (i in 1:n) {
  trozar <- M7seq[starts[i]:(starts[i]+499)] 
  trozarGC <- GC(trozar)
  print(trozarGC)
  trozarGCs[i] <- trozarGC #agrego esta línea
}
#Para plotear el contenido de GC a lo largo de la secuencia
plot(starts,trozarGCs,type="b",xlab="Posición inicial de nucleótido",ylab="Contenido de GC")
#Llamar a una BD
choosebank()
choosebank("genbank")
