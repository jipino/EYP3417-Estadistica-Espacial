## ============================================================
##  AYUDANTÍA 09 — EYP3417 Estadística Espacial
##  Tema: Procesos Puntuales — Intensidad, PPP, Funciones K/L y Modelos Log-lineales
## ============================================================

library(spatstat)
library(ggplot2)
library(patchwork)
library(viridisLite)

set.seed(2026)


## ============================================================
##  PROBLEMA 1: INTENSIDAD ESPACIAL Y EL PROCESO DE POISSON (PPP)
## ============================================================
## ------------------------------------------------------------
## 1(b) Simulación de PPP homogéneo y no homogéneo en [0,1]^2
## ------------------------------------------------------------

## Simulamos la ventana solicitada
ventana <- owin(xrange = c(0, 1), yrange = c(0, 1))

## PPP homogéneo con intensidad constante λ = 100
ppp_homogeneo <- rpoispp(lambda = 100, win = ventana)

## PPP no homogéneo con función de intensidad λ(x, y) = 200 · e^(-3x)
lambda_no_hom <- function(x, y) { 200 * exp(-3 * x) }
ppp_no_homogeneo <- rpoispp(lambda = lambda_no_hom, win = ventana)

cat(sprintf("PPP homogéneo:    n = %d puntos simulados (λ = 100 esperado)\n",
            ppp_homogeneo$n))
cat(sprintf("PPP no homogéneo: n = %d puntos simulados (λ(x,y) = 200·e^(-3x))\n\n",
            ppp_no_homogeneo$n))

## Podemos ver de la simulación, que ambos conteos son consistentes con lo esperado bajo el modelo,
## el número total de puntos de un PPP es Poisson con media μ(S) = ∫_S λ(u) du.

## Para el caso homogéneo tenemos μ([0,1]²) = λ · área = 100 · 1 = 100 →  n = 105, es decir, la media
## está dentro de la variabilidad típica de una Poisson(100), cuya desv. est. es 10.

## Para el caso no homogéneo μ([0,1]²) = ∫₀¹∫₀¹ 200·e^(-3x) dx dy = (200/3)(1 − e^(-3)) ≈ 66.3  →  n = 50
## un valor algo por debajo de la media, pero perfectamente plausible para una Poisson(66.3) (desv. est ≈ 8.1)

## Es completamente esperable que el proceso no homogéneo genere muchos menos puntos en total que el homogéneo,
## ya que su intensidad λ(x,y) = 200·e^(-3x) promedia apróx. 66 en toda la ventana, muy por debajo de la intensidad
## constante λ = 100 del proceso homogéneo

## Gráficos comparativos de ambas configuraciones
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
plot(ppp_homogeneo, main = "PPP homogéneo: λ = 100",
     pch = 21, cex = 0.9, cols = "#2c3e8c", bg = "#5b7fd4")
plot(ppp_no_homogeneo, main = "PPP no homogéneo: λ(x,y) = 200·e^(-3x)",
     pch = 21, cex = 0.9, cols = "#9c1f1f", bg = "#e0584f")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## Se aprecia del gráfico que:

## El PPP Homogéneo (izquierdo): los puntos cubren la ventana de manera
## aproximadamente uniforme, es decir, no se aprecian zonas con mayor o
## menor densidad de puntos más allá de la fluctuación aleatoria propia del
## muestreo. Esto es justamente lo que se espera de un λ(u) = cte.

## Para el caso de PPP no homogéneo (derecha): se observa una clara concentración
## de puntos hacia el borde izquierdo del cuadrado, y una progresiva dilución (puntos
## cada vez más escasos) a medida que x crece hacia 1. El gradiente de densidad de
## puntos de izquierda a derecha es visualmente evidente y replica fielmente la forma
## exponencial decrediente de λ en la dirección x.

## ------------------------------------------------------------
## 1(c) Estimación no paramétrica de la intensidad (kernel gaussiano)
## ------------------------------------------------------------

## Estimación con bandwidth razonable (seleccionado por regla de Diggle)
bw_optimo <- bw.diggle(ppp_no_homogeneo)
intensidad_kernel <- density.ppp(ppp_no_homogeneo, sigma = bw_optimo, kernel = "gaussian")

## Estimaciones con bandwidth muy pequeño y muy grande para comparar
intensidad_chica  <- density.ppp(ppp_no_homogeneo, sigma = bw_optimo / 5,  kernel = "gaussian")
intensidad_grande <- density.ppp(ppp_no_homogeneo, sigma = bw_optimo * 5,  kernel = "gaussian")

cat(sprintf("Bandwidth óptimo (regla de Diggle): sigma = %.4f\n\n", bw_optimo))

## La regla de Diggle nos entrega σ ≈ 0.0555, que es un bandwidth bastante pequeño en
## relación al tamño de la ventana [0,1]² (≈ 6% de su lado), lo que es razonable dado
## que solo se cuenta con n = 50 puntos para estimar la superficie de intensidad. Con pocos
## datos conviene un suavizado moderado que no diluya por completo la información disponible.
## Los paneles pequeño y grande usan σ/5 ≈ 0.0111 y σ·5 ≈ 0.2775 respectivamente, para exhibir
## los dos extremos del trade-off sesgo-varianza


paleta_intensidad <- viridis(128, option = "plasma")

par(mfrow = c(1, 3), mar = c(1, 1, 3, 3))
plot(intensidad_chica,  main = sprintf("Bandwidth pequeño (σ = %.3f)", bw_optimo / 5),
     col = paleta_intensidad, ribargs = list(las = 1))
points(ppp_no_homogeneo, pch = 21, cex = 0.5, col = "white", bg = "white")
plot(intensidad_kernel, main = sprintf("Bandwidth óptimo (σ = %.3f)", bw_optimo),
     col = paleta_intensidad, ribargs = list(las = 1))
points(ppp_no_homogeneo, pch = 21, cex = 0.5, col = "white", bg = "white")
plot(intensidad_grande, main = sprintf("Bandwidth grande (σ = %.3f)", bw_optimo * 5),
     col = paleta_intensidad, ribargs = list(las = 1))
points(ppp_no_homogeneo, pch = 21, cex = 0.5, col = "white", bg = "white")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

# INTERPRETACIÓN VISUAL Y EMPÍRICA DEL EFECTO DEL BANDWIDTH (con σ = 0.0597
# como referencia central, σ/5 ≈ 0.012 y σ·5 ≈ 0.299 como extremos):

## Para el Bandwidth muy pequeño:
## Tenemos que el estimador de intensidad queda dominado por el ruido muestral, es
## decir, aparecen muchos picos aislados y agudos (manchas amarillas puntuales) centrados
## casi exactamente en cada punto observado, rodeados de un fondo azul oscuro de intensidad
## prácticamente nula (sobreajuste). El mapa resultante es bastante irregular y no logra capturar
## el patrón suave subyacente λ(x, y) = 200·e^(-3x); en cambio retrata casi literalmente la nube
## de puntos individual. Alta varianza, bajo sesgo.

## Bandwidth muy grande: 
## El estimador sobresuaviza, pues la información se promedia sobre regiones tan amplias que el
## patrón de decaimiento exponencial en x se reduce a un degradado suave y casi lineal de izquierda
## a derecha, perdiendo el detalle de la verdadera forma exponencial y estirando la influencia de
## los puntos mucho más allá de su entorno real (alto sesgo, baja varianza)

## Bandwidth óptimo: Logra un balance sesgo-varianza adecuado, puesel mapa de intensidad estimado
## muestra manchas suaves de alta intensidad concentradas hacia el borde izquierdo (x ≈ 0) que decaen gradualmente hacia la derecha 
## recuperando visualmente  el patrón de decaimiento exponencial en la dirección x, sin verse dominado por fluctuaciones puntuales del 
## muestreo ni perder la estructura real por exceso de promediado.

## En conclusión la elección del bandwidth es crítica, pues un valor muy pequeño conduce a sobreajuste (ruido, imagen de los propios puntos)
## y uno muy grande conduce a un subajuste (pérdida de estructura, superficie casi constante). Ambos extremos distorsionan la lectura
## del fenómeno espacial subyacente. El valor entregado por bw.diggle() ofrece un compromiso razonable entre ambos extremos para este
## patrón.


## ============================================================
##  PROBLEMA 2: DEPENDENCIA DE SEGUNDO ORDEN — FUNCIONES K Y L
## ============================================================

## ------------------------------------------------------------
## 2(b) Estimación de L(h) - h para el conjunto de datos 'cells'
## ------------------------------------------------------------

data(cells)

cat(sprintf("Conjunto 'cells': n = %d puntos (centros celulares) en ventana unitaria\n\n",
            cells$n))

L_cells <- Lest(cells, correction = "Ripley")

plot(L_cells, . - r ~ r, main = "Función L(h) - h para el patrón 'cells'",
     legend = TRUE, xlab = "h (distancia)", ylab = "L(h) - h")

## Podemos ver que la curva empírica L^(h) - h (línea negra) es no monótona, teniendo:
## 1. h pequeño (h ≲ 0.15): Cae abruptamente a mín ~ -0.09 en h ≈ 0.10. 
## Indica REPULSIÓN fuerte (inhibición de corto alcance tipo hard-core, d_min ≈ 0.10).
## 2. h intermedio (h ≈ 0.18-0.21): Cruza el eje 0 en h ≈ 0.17 y alcanza máx ~ +0.005.
## Indica leve EXCESO/AGRUPACIÓN regular (anillo de vecinos próximos por acomodo).
## 3. h grande (h ≳ 0.22): Vuelve a valores levemente negativos (~ -0.016) en h = 0.25.
## Fluctuación típica de CSR esperable por el tamaño muestral reducido (n=42).

## ------------------------------------------------------------
## 2(c) Interpretación: ¿es razonable un modelo hard-core?
## ------------------------------------------------------------

## Del gráfico de b) podemos ver:
## 1. Comportamiento: Negativo a distancias cortas (mín ~ -0.09 en h ≈ 0.10), 
## rebota sobre cero (máx ~ +0.005 en h ≈ 0.20) y oscila cerca de CSR.
## 2. Sentido biológico: Las células tienen volumen y no pueden solaparse; 
## la geometría celular impone una distancia mínima de exclusión.
## 3. Conclusión: Es razonable un modelo HARD-CORE (núcleo duro) con distancia 
## de inhibición δ ≈ 0.10, consistente con la física del proceso y el déficit 
## de pares a escala corta.


## ============================================================
##  PROBLEMA 3: MODELAMIENTO PARAMÉTRICO Y VALIDACIÓN MONTE CARLO
## ============================================================

## ------------------------------------------------------------
## 3(b) Ajuste del modelo log-lineal con ppm() e interpretación de θ̂2
## ------------------------------------------------------------

modelo_loglineal <- ppm(ppp_no_homogeneo ~ x)

cat("--- Modelo ajustado: log λ(x, y) = θ1 + θ2 · x ---\n")
print(coef(summary(modelo_loglineal)))

theta_hat <- coef(modelo_loglineal)
theta1_hat <- theta_hat[["(Intercept)"]]
theta2_hat <- theta_hat[["x"]]

cat(sprintf("θ̂1 (intercepto) = %.4f  →  λ̂(0, y) = exp(θ̂1) = %.2f\n",
            theta1_hat, exp(theta1_hat)))
cat(sprintf("θ̂2 (pendiente)  = %.4f\n\n", theta2_hat))

## Podemos ver que de θ2:

## El signo es consistente con el modelo real. Confirma que la intensidad decrece
## al aumentar x. Respecto a la magnitud, podemos decir que recupera bien la estructura
## poblacional, pues el IC contiene al valor real, lo cual es coherente con lo teórico.
## Finalmente, el Test de Wald rechaza H0 (theta2 = 0), luego la dependencia espacial
## de la intensidad respecto a x es altamente significativa.

## ------------------------------------------------------------
## 3(c) Validación mediante bandas de Monte Carlo para L(h) - h
## ------------------------------------------------------------

## Bandas de confianza al 95% bajo H0 (PPP no homogéneo con intensidad estimada).
## Simulación de Monte Carlo (99 realizaciones): envelope puntual con nivel
## de significancia global aproximado de alfa = 2 / (99 + 1) = 0.02.
## Evalúa si la estructura de segundo orden se explica solo por la intensidad espacial x.
envelope_L <- envelope(modelo_loglineal, fun = Lest, nsim = 99,
                       transform = expression(. - r),
                       correction = "Ripley", verbose = FALSE)

plot(envelope_L, main = "Bandas de Monte Carlo (95%) para L(h) - h\nbajo el modelo log-lineal ajustado",
     xlab = "h (distancia)", ylab = "L(h) - h")

## Del gráfico podemos ver que hay una tendencia creciente, tanto en la curva observada, como la
## media bajo H0 (roja), las cuales crecen con h. Es un efecto espurio de la no homogeneidad de primer
## orden, no se evidencia agregación.

## Respecto al comportamiento, podemos ver que la curva observada sigue de cerca a la media simulada
## y se mantiene siempre dentro de las bandas de confianza.

## En conclusión: No hay evidencia para rechazar H0. La variación espacial se explica por heterogeneidad
## de primer orden. No se requiere interacción de segundo orden. El ejercicio valida que no se confunde
## homogeneidad con interacción real.