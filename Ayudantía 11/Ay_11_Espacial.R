## ============================================================
##  AYUDANTÍA 11 — EYP3417 Estadística Espacial
##  Tema: Procesos Puntuales — PPP, Intensidad Kernel, CSR y Modelo Log-lineal
##  Ejercicio 3: Trampas Cámara en Bosque Templado (datos simulados)
## ============================================================

library(spatstat)
library(viridisLite)

set.seed(2026)


## ============================================================
##  SIMULACIÓN DE DATOS
## ============================================================

## Simulamos avistamientos de fauna en una zona boscosa de 2 × 2 km.
## La función de intensidad es log λ(x,y) = θ₁ + θ₂·x = -9.0 - 0.001·x,
## lo que produce mayor actividad faunística hacia el borde oeste (x = 0).

ventana      <- owin(xrange = c(0, 2000), yrange = c(0, 2000))
lambda_fauna <- function(x, y) { exp(-9.0 - 0.0010 * x) }
camaras_ppp  <- rpoispp(lambda = lambda_fauna, win = ventana)

cat(sprintf("Avistamientos simulados:          n = %d puntos\n", camaras_ppp$n))
cat(sprintf("Intensidad esperada total: μ(S) ≈ %.1f eventos\n\n",
            integral(as.im(lambda_fauna, W = ventana))))

## Se simularon n = 216 avistamientos, muy cercanos a los μ(S) ≈ 213.4 esperados bajo el
## modelo, lo que confirma que la variabilidad de conteo es consistente con Poisson(213.4).


## ============================================================
##  EJERCICIO 3: ANÁLISIS ESPACIAL DE TRAMPAS CÁMARA
## ============================================================

## ------------------------------------------------------------
##  Visualización de la intensidad teórica y el patrón simulado
## ------------------------------------------------------------

lambda_im <- as.im(lambda_fauna, W = ventana)

par(mfrow = c(1, 2), mar = c(3, 3, 3, 3))
plot(lambda_im,
     col  = viridis(128, option = "plasma"),
     main = "Intensidad teórica λ(x, y)",
     ribargs = list(las = 1))
plot(camaras_ppp,
     main = "Patrón simulado — Trampas Cámara",
     pch = 21, cex = 0.8, cols = "#1a3d1c", bg = "#5a9e5d")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## La intensidad teórica muestra un gradiente pronunciado de oeste (amarillo, λ ≈ 1.2×10⁻⁴)
## a este (azul oscuro, λ ≈ 1.7×10⁻⁵). El patrón simulado reproduce fielmente este gradiente:
## mayor densidad de puntos en el borde izquierdo, dilución progresiva hacia la derecha.


## ------------------------------------------------------------
## 3(a) Construcción del objeto ppp y visualización
## ------------------------------------------------------------

cat(sprintf("Patrón de puntos: n = %d avistamientos en ventana [0, 2000]² m\n\n",
            camaras_ppp$n))

par(mfrow = c(1, 1), mar = c(2, 2, 3, 2))
plot(camaras_ppp,
     main = "Avistamientos de fauna — Trampas Cámara (ventana 2 × 2 km)",
     pch = 21, cex = 0.9, cols = "#1a3d1c", bg = "#5a9e5d")

## El patrón presenta concentración visible de puntos hacia el borde oeste, con zonas
## progresivamente más vacías hacia el este. Esto anticipa un PPP no homogéneo con
## intensidad decreciente en x y motiva el modelo log-lineal del apartado (e).


## ------------------------------------------------------------
## 3(b) Estimación kernel de la intensidad con tres bandwidths
## ------------------------------------------------------------

bw_optimo <- bw.diggle(camaras_ppp)

intensidad_chica  <- density.ppp(camaras_ppp, sigma = bw_optimo / 5, kernel = "gaussian")
intensidad_kernel <- density.ppp(camaras_ppp, sigma = bw_optimo,      kernel = "gaussian")
intensidad_grande <- density.ppp(camaras_ppp, sigma = bw_optimo * 5, kernel = "gaussian")

cat(sprintf("Bandwidth óptimo (regla de Diggle): sigma = %.2f m\n\n", bw_optimo))

paleta_intensidad <- viridis(128, option = "plasma")

par(mfrow = c(1, 3), mar = c(1, 1, 3, 3))
plot(intensidad_chica,
     main = sprintf("Bandwidth pequeño\n(σ = %.1f m)", bw_optimo / 5),
     col = paleta_intensidad, ribargs = list(las = 1))
points(camaras_ppp, pch = 21, cex = 0.4, col = "white", bg = "white")
plot(intensidad_kernel,
     main = sprintf("Bandwidth óptimo\n(σ = %.1f m)", bw_optimo),
     col = paleta_intensidad, ribargs = list(las = 1))
points(camaras_ppp, pch = 21, cex = 0.4, col = "white", bg = "white")
plot(intensidad_grande,
     main = sprintf("Bandwidth grande\n(σ = %.1f m)", bw_optimo * 5),
     col = paleta_intensidad, ribargs = list(las = 1))
points(camaras_ppp, pch = 21, cex = 0.4, col = "white", bg = "white")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## σ = 24.4 m (pequeño): sobreajuste evidente — picos aislados centrados en cada punto,
## fondo azul oscuro de intensidad casi nula. El mapa retrata el ruido muestral,
## no la estructura subyacente (alta varianza, bajo sesgo).

## σ = 121.8 m (óptimo, bw.diggle): el mapa revela el gradiente oeste-este de la
## intensidad real, con manchas suaves de mayor actividad en el borde izquierdo.
## Es el único que recupera visualmente la forma de λ(x,y) = exp(-9 - 0.001x).

## σ = 609.1 m (grande): sobresuavizado — el gradiente se diluye en una superficie
## casi constante de izquierda a derecha. Se pierde la heterogeneidad real (alto sesgo).


## ------------------------------------------------------------
## 3(c) Tests de CSR: funciones G y F con bandas de Monte Carlo
## ------------------------------------------------------------

n_sim      <- 99
alpha      <- 0.05
rank_value <- alpha * (n_sim + 1)

set.seed(2026)
G_env <- envelope(camaras_ppp, fun = Gest, correction = "km",
                  nsim = n_sim, rank = rank_value, verbose = FALSE)

set.seed(2026)
F_env <- envelope(camaras_ppp, fun = Fest, correction = "km",
                  nsim = n_sim, rank = rank_value, verbose = FALSE)

cat(sprintf("Nivel efectivo del test Monte Carlo: 2 / (%d + 1) = %.4f\n\n",
            n_sim, 2 / (n_sim + 1)))

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
plot(G_env, main = "Función G — Trampas Cámara", xlab = "r (metros)", ylab = "G(r)")
plot(F_env, main = "Función F — Trampas Cámara", xlab = "r (metros)", ylab = "F(r)")
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))

## Nivel efectivo: 2/100 = 0.02, más conservador que α = 0.05.

## Función G: la curva observada se mantiene dentro de las bandas de confianza en
## todo el rango recomendado. No se rechaza CSR desde la perspectiva de distancias
## entre vecinos.

## Función F: la curva observada también permanece dentro de las bandas. No hay
## evidencia de exceso de espacio vacío ni de ocupación más uniforme que CSR.

## Ambas funciones son consistentes: no rechazan CSR a escala local (r pequeño).
## Esto es esperable para un PPP — la no homogeneidad es un efecto de primer orden
## que G y F con bandas bajo CSR homogéneo no alcanzan a capturar de forma robusta
## a las escalas de vecindad más próximas.


## ------------------------------------------------------------
## 3(d) Análisis de segundo orden: función L(h) - h
## ------------------------------------------------------------

set.seed(2026)
L_env <- envelope(camaras_ppp, fun = Lest, correction = "Ripley",
                  nsim = n_sim, rank = rank_value,
                  transform = expression(. - r), verbose = FALSE)

plot(L_env,
     main = "Función L(h) - h con bandas de Monte Carlo (CSR)\nTrampas Cámara",
     xlab = "h (metros)", ylab = "L(h) - h")
abline(h = 0, lty = 2, col = "gray50")

## Para h ≲ 100 m la curva oscila dentro de las bandas, compatible con CSR a escala local.
## A partir de h ≈ 150 m la curva sale y se aleja progresivamente de las bandas, alcanzando
## L(h) - h ≈ 65 en h = 500 m. Esto indica un aparente exceso de pares a distancias medias
## y grandes que parece sugerir clustering.
## Sin embargo, esta desviación NO es evidencia de interacción entre puntos: es un artefacto
## conocido de la heterogeneidad de primer orden. El envelope bajo CSR homogéneo no contempla
## la variación de λ(x,y), por lo que sobreestima los pares esperados a gran escala.
## El apartado (e) permite discriminar entre ambas hipótesis ajustando por la no homogeneidad.


## ------------------------------------------------------------
## 3(e) Modelo log-lineal de intensidad y validación Monte Carlo
## ------------------------------------------------------------

modelo_loglineal <- ppm(camaras_ppp ~ x)

cat("--- Modelo ajustado: log λ(x, y; θ) = θ̂₁ + θ̂₂·x ---\n")
print(coef(summary(modelo_loglineal)))

theta_hat  <- coef(modelo_loglineal)
theta1_hat <- theta_hat[["(Intercept)"]]
theta2_hat <- theta_hat[["x"]]

cat(sprintf("\nθ̂₁ = %7.4f  [real: -9.0000]  →  λ̂(0,y) = %.5f eventos/m²\n",
            theta1_hat, exp(theta1_hat)))
cat(sprintf("θ̂₂ = %7.4f  [real: -0.0010]  →  cada 100 m adicionales en x multiplican λ por %.4f\n\n",
            theta2_hat, exp(100 * theta2_hat)))

## θ̂₁ = -8.9285 (real: -9.0) y θ̂₂ = -0.0011 (real: -0.0010): ambos estimadores
## recuperan los valores poblacionales con buena precisión. Los ICs al 95% contienen
## los valores reales en ambos casos.
## θ̂₂ < 0 y altamente significativo (Z = -8.26, p < 0.001): la intensidad de
## avistamientos decrece significativamente con la coordenada este. Por cada 100 m
## adicionales en x, λ se reduce en un factor de exp(100·θ̂₂) ≈ 0.897, es decir,
## aproximadamente un 10% menos de avistamientos esperados por cada 100 m hacia el este.

set.seed(2026)
envelope_L_mod <- envelope(modelo_loglineal, fun = Lest, nsim = n_sim,
                            transform = expression(. - r),
                            correction = "Ripley", verbose = FALSE)

plot(envelope_L_mod,
     main = "Bandas de Monte Carlo (95%) para L(h) - h\nbajo el modelo log-lineal ajustado",
     xlab = "h (metros)", ylab = "L(h) - h")
abline(h = 0, lty = 2, col = "gray50")

## La curva observada (negra) sigue de cerca la media simulada (roja) y permanece
## dentro de las bandas grises en todo el rango de h. La tendencia creciente de
## L(h) - h es compartida por el envelope, confirmando que es un artefacto esperado
## bajo el propio modelo no homogéneo ajustado, no evidencia de clustering real.
## Conclusión: no se rechaza el modelo log-lineal. La estructura espacial del patrón
## queda completamente explicada por la heterogeneidad de primer orden en x, sin
## necesidad de invocar interacción de segundo orden entre puntos — exactamente lo
## esperado dado que los datos fueron generados como un PPP.
