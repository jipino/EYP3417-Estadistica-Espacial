## Ayudantía 10 - Estadística Espacial

## Librerías
library(spatstat)
library(spatstat.geom)
library(spatstat.explore)

## Ejercicio 01

## Cargar datos
data(japanesepines)
data(cells)
data(redwood)

## Descripción breve de los conjuntos de datos utilizados

## japanesepines:
## Contiene las ubicaciones de 65 brotes de pino negro japonés (Pinus thunbergii) recolectados en un
## bosque natural. El área original (5.7 x 5.7 m) fue reescalada al cuadrado unitario. Es un ejemplo
## clásico de patrón compatible con aleatoriedad espacial completa (CSR).

## cells:
## Registra las posiciones de los centros de 42 células biológicas observadas por microscopía óptica.
## También reescalado al cuadrado unitario, este conjunto de datos es un ejemplo canónico de patrón
## inhibido o con repulsión entre puntos (ordenado).

## redwood:
## Muestra las posiciones de 62 brotes y plántulas de Sequoia gigante (Sequoiadendron giganteum)
## en una subregión cuadrada (~63.1 pies de lado), reescalada al cuadrado unitario. Es un caso clásico
## de patrón con agrupamiento o clustering espacial.

## Fuente común: Todos los datos han sido utilizados ampliamente en la literatura (Ripley, Diggle)
## para ilustrar distintas configuraciones espaciales: aleatoria (japanesepines), repulsiva (cells),
## y agrupada (redwood).


## Verificar clases, usaremos la función is.ppp(), la cual nos
## permite saber si un objeto es del tipo point pattern planar, el cual
## corresponde al formato estándar para representar patrones de puntos en 2D.

## Vamos a combinar la función anterior, con un is.ppp(), en el caso de que alguno
## sea FALSE, se tendrá un error, en caso contrario no ocurre nada.

stopifnot(is.ppp(japanesepines), is.ppp(cells), is.ppp(redwood))

## No ocurrió nada, luego podemos continuar con el análisis de nuestros conjuntos
## de datos.

## Visualización de los patrones
par(mfrow = c(1,3), mar = c(2,2,2,2))
plot(japanesepines, main = "Japanese Pines")
plot(cells, main = "Cells")
plot(redwood, main = "Redwood")
par(mfrow = c(1,1))

## Sobre el conjunto de datos de japanesepines, podemos apreciar que se tiene una
## distribución aparentemente aleatoria, en la que los puntos no se ven ni muy
## juntos ni muy separados, luego como hipótesis podemos establecer que es
## compatible con aleatoriedad espacial completa (CSR)

## En segundo lugar, si hablamos del conjunto de datos de cells, se observa una
## distribución aparentemente ordenada y regular. Los puntos están equiespaciados
## luego como hipótesis podemos establecer que existe inhibición o repulsión
## entre puntos.

## Finalmente, respecto al conjunto de datos redwood, se aprecia una clusterización
## de los puntos, es decir, hay zonas con alta concentración de puntos y otras vacías.
## Luego, como hipótesis, podemos decir que hay agrupamiento espacial.

## En primer lugar, realizaremos el test G, para ver si existe aleatoriedad 
## espacial completa (CSR)

## Parámetros: numero de simulaciones, nivel de significancia y valor
## del ranking usado para los percentiles en las bandas de confianza.

n_sim <- 100
alpha <- 0.05
rank_value <- alpha * (n_sim + 1)

## Test de aleatoriedad espacial completa con función G:

## Aplicamos la función G (distancia al vecino más cercano) sobre tres patrones:
## japanesepines, cells y redwood

## Para cada patrón tenemos:
## 1. Se genera una banda de confianza a partir de 100 simulaciones de CSR.
## 2. Se calcula la curva observada G(r).
## 3. Se compara con la banda de confianza construida a partir de las simulaciones.

## El parámetro 'rank' está ajustado a nivel de significancia α = 0.05.
## Si la curva observada se encuentra fuera de la banda de confianza, se rechaza CSR.
## Si está completamente dentro, entonces no hay evidencia para rechazar CSR.

set.seed(13052025)
jap_Genv <- envelope(japanesepines, fun = Gest, correction = "all", nsim = n_sim, rank = rank_value)
cells_Genv <- envelope(cells, fun = Gest, correction = "all", nsim = n_sim, rank = rank_value)
red_Genv <- envelope(redwood, fun = Gest, correction = "all", nsim = n_sim, rank = rank_value)

## También se graficará para complementar
par(mfrow = c(1,3))
plot(jap_Genv, main = "Función G - Japanese Pines")
plot(cells_Genv, main = "Función G - Cells")
plot(red_Genv, main = "Función G - Redwood")

## Para cada patrón de puntos se generaron 100 simulaciones bajo el modelo
## CSR y se construyeron bandas de confianza al 95%
## Significancia del test (Monte Carlo): 2 / (100 + 1) = 0.0198

## jap_Genv: El rango recomendado de r es [0,0.104], la curva observada G(r)
## cae completamente dentro de la banda de confianza, luego no se rechaza
## el CSR, por lo que el patrón es compatible con aleatoriedad espacial completa.

## cells_Genv: El rango recomendado de r es [0,0.15], la curva observada G(r) cae
## por debajo de la banda de confianza en buena parte del rango, luego se rechaza
## CSR y se indica repulsión entre puntos, pues la curva está bajo las bandas de
## confianza

## red_Genv: El rango recomendado de r es [0,0.0598], la curva observada G(r) cae
## se ubica sobre las bandas de confianza, luego se rechaza la CSR y se detecta
## agrupamiento (clustering espacial)

## Test función F 
set.seed(13052025)
jap_Fenv <- envelope(japanesepines, fun = Fest, correction = "all", nsim = n_sim, rank = rank_value)
cells_Fenv <- envelope(cells, fun = Fest, correction = "all", nsim = n_sim, rank = rank_value)
red_Fenv <- envelope(redwood, fun = Fest, correction = "all", nsim = n_sim, rank = rank_value)

## Visualizamos función F
plot(jap_Fenv, main = "Función F - Japanese Pines")
plot(cells_Fenv, main = "Función F - Cells")
plot(red_Fenv, main = "Función F - Redwood")

## El test F (espacio vacío) se usa para evaluar si el patrón de puntos observado
## presenta agrupamiento, inhibición o CSR. Para cada uno de los tres patrones
## de puntos se simularon 100 patrones CSR, se calcula la función F observada
## y sus bandas de confianza al 95%, la significancia es de 2 / (100 + 1) = 0.0198
## con una corrección de borde: Kaplan-Meier (km)

## jap_Fenv: Vemos que el rango recomendado de r es [0, 0.111], la curva observada
## F(r) está completamente contenida en las bandas de confianza, luego no se rechaza
## CSR, por lo que el patrón parece aleatorio respecto al espacio vacío.

## cells_Fenv: El rango recomendado de r es [0,0.086], la curva observada F(r) está
## por encima de las bandas de confianza desde valores bajos de r, luego se rechaza CSR
## de ello, hay evidencia de inhibición (los puntos ocupan el espacio uniformemente)

## red_Fenv: El rango recomendado es [0, 0.154], con una curva observada F(r) que está
## por debajo de las bandas de confianza, se rechaza CSR y se evidencia el agrupamiento (más
## espacio vacío del esperado por CSR.)

## Ejercicio 02

## Cargamos los datos desde github
url <- "https://raw.githubusercontent.com/mgimond/Spatial/main/Data/ppa.RData"
load(url(url, "rb"))

## Visualización del patrón sobre el estado de Massachusetts

## Primero asignamos el área de estudio (ma) al objeto starbucks, esto es relevante,
## pues starbucks es un ppp (patrón de puntos) y posteriormente grafiquemos la
## distribución de ellos en el estado de Massachusetts


par(mfrow = c(1,1),mar = c(3, 3, 3, 3))  
Window(starbucks) <- ma
plot(ma,main = "Distribución de tiendas Starbucks en Massachusetts",
     col = "white",border = "black",lwd = 1.2,axes = FALSE)         

points(starbucks, col = "gold2", pch = 19, cex = 0.7) 
legend("bottomleft", legend = "Starbucks", pch = 19, col = "gold2", bty = "n")

## El gráfico muestra la ubicación espacial de 171 tiendas de Starbucks dentro del
## estado de Massachusetts, se aprecia una fuerte concentración de puntos en la región
## noreste del estado. A medida que nos alejamos de esa zona, la densidad de tiendas
## disminuye visiblemente, con algunas áreas del estado mostrando vacíos espaciales.
## Lo anterior sugiere un patrón de agrupamiento, en el que las tiendas se concentran
## en zonas urbanas densamente pobladas.

## A partir de la visualización realizada, se puede anticipar que el patrón espacial
## no es compatible con aleatoriedad espacial completa (CSR). Se espera que las 
## funciones K y L evidencien agrupamiento, con curvas observadas por encima de las
## bandas de simulación CSR.

## Función K - cuenta acumulada de vecinos en círculos
set.seed(13052025)
star_Kenv <- envelope(starbucks, fun = Kest, correction = "all", nsim = n_sim, rank = rank_value)
plot(star_Kenv, main = "Función K - Starbucks")

## La función K evalúa la cantidad acumulada de puntos vecinos dentro de un radio
## r alrededor de cada punto del patrón, comparado con lo esperado bajo CSR.

## Del gráfico, tenemos: La curva negra es K observada del patrón de
## Starbucks. La curva roja punteada es K teórica bajo CSR, el área gris son
## las bandas de confianza de simulaciones CSR (100 reps.) y se aplica una
## corrección de borde isotrópica ("iso")

## Por otro lado, del test: La significancia de Monte Carlo es 2 / (100 + 1) = 0.0198
## y el rango recomendado de distancias es de [0, 38488] (en unidades UTM)

## Del gráfico, vemos que la curva observada está completamente por encima de las bandas
## de confianza generadas bajo CSR. Lo anterior indica un número de vecinos acumulados
## significativamente mayor a lo esperado bajo aleatoriedad.

## En conclusión, se rechaza la hipótesis de CSR, luego existe evidencia clara de 
## agrupamiento espacial, es decir, las tiendas Starbucks tienden a localizarse más
## cerca entre sí de lo que ocurriría si estuvieran distribuidas aleatoriamente en
## el espacio.



## Función L 
set.seed(13052025)
star_Lenv <- envelope(starbucks, fun = Lest, correction = "all", nsim = n_sim, rank = rank_value)
plot(star_Lenv, main = "Función L - Starbucks")


## La función L es una transformación de la función K que facilita la interpretación
## al linealizar el crecimiento esperado bajo CSR.

## Del gráfico tenemos: La curva negra, L observada del patrón de Starbucks, la
## curva roja punteada, L teórica bajo CSR, el área gris son las bandas de confianza
## generadas por 100 simulaciones CSR y se aplica una corrección de borde aplicada
## isotrópica ("iso")

## Del test: Nivel de significancia Monte Carlo: 2 / 101 = 0.0198, con un rango
## recomendado de r [0, 38488] (en unidades UTM)

## Del gráfico, se puede ver que la curva está consistentemente por encima de las
## bandas de confianza a lo largo de todo el rango, esto indica que las distancias
## acumuladas entre puntos son más cortas de lo esperado bajo CSR, lo que implica
## un exceso de cercanía entre puntos

## En conclusión, se rechaza la hipótesis de CSR, existe una fuerte evidencia de
## agrupamiento espacial en la distribución de tiendas Starbucks de Massachussetts.
## Este resultado es consistente con el test de la función K y con la visualización
## inicial del patrón espacial.
