## Ay_04_Espacial.R -------------------------------------------------------
## EYP3417 — Estadística Espacial
## Ayudantía 04: Kriging Simple y Ordinario — Acuífero Wolfcamp
## Profesor: Alfredo Alegría  |  Ayudante: Juan Pino

# 0. Setup ----------------------------------------------------------------
library(geoR)
library(gstat)
library(tidyverse)
library(patchwork)

## Retomamos el flujo completo de las sesiones anteriores: tendencia lineal
## → residuos → objeto geodata → estimación por MV. El código es idéntico
## al de Ay_03 para mantener consistencia en los parámetros estimados.
data(wolfcamp)

df_wolf <- data.frame(
  x        = wolfcamp$coords[, 1],
  y        = wolfcamp$coords[, 2],
  pressure = wolfcamp$data
)

model_trend <- lm(pressure ~ x + y, data = df_wolf)

data_res <- df_wolf |>
  mutate(residuals = residuals(model_trend))

geo_res <- as.geodata(data_res, coords.col = c(1, 2), data.col = 4)

## Re-estimamos θ por MV — mismo modelo que Ay_03 (Parte c)
fit_ml <- suppressMessages(
  likfit(
    geodata      = geo_res,
    cov.model    = "exponential",
    ini.cov.pars = c(2000, 100),
    nugget       = 500,
    lik.method   = "ML"
  )
)

# Parte a: Grilla de predicción y vector θ = (σ², φ, τ²) ----------------
## El primer paso del kriging es definir el dominio sobre el que vamos a
## predecir. Construimos una grilla regular que cubra el bounding box de
## los datos de wolfcamp con 50 × 50 = 2.500 puntos de predicción.
##
## Luego extraemos los parámetros estimados en Ay_03:
##   σ² (sigmasq): varianza del proceso (sill parcial)
##   φ  (phi):     escala del Exponencial; rango práctico = 3φ
##   τ² (tausq):   nugget (variabilidad a escala sub-muestral o error de medición)
##   μ  (beta):    media estimada por GLS — usada como "media conocida" en KS

x_seq    <- seq(floor(min(df_wolf$x)),   ceiling(max(df_wolf$x)),   length.out = 50)
y_seq    <- seq(floor(min(df_wolf$y)),   ceiling(max(df_wolf$y)),   length.out = 50)
grid_df  <- expand.grid(x = x_seq, y = y_seq)
grid_mat <- as.matrix(grid_df)

## Parámetros θ
sigmasq <- fit_ml$sigmasq
phi     <- fit_ml$phi
tausq   <- fit_ml$tausq
mu_ks   <- fit_ml$beta[1]   ## media GLS de los residuos (~14 ft) — "conocida" en KS

cat("=== Vector θ estimado por MV ===\n")
cat("σ²      =", round(sigmasq, 2), "ft²\n")
cat("φ       =", round(phi,     2), "km   (rango práctico = 3φ =", round(3*phi, 1), "km)\n")
cat("τ²      =", round(tausq,   2), "ft²\n")
cat("μ (GLS) =", round(mu_ks,   4), "ft\n")
cat("Grilla:  ", nrow(grid_mat), "puntos de predicción\n\n")

## Los parámetros son los mismos de Ay_03, confirmados aquí:
##
##   σ² = 3.782 ft²: varianza del proceso espacial de los residuos de
##     tendencia. Representa la variabilidad total del campo ε(s) una vez
##     removida la tendencia lineal SW → NE. Es la "altura" del sill en el
##     variograma: γ(h) → σ² + τ² = 3782 + 673 ≈ 4.455 ft² cuando h → ∞,
##     consistente con lo observado en el variograma experimental de Ay_03.
##
##   φ = 37.28 km: parámetro de escala del modelo Exponencial. El rango
##     práctico 3φ ≈ 111.8 km es la distancia a la que la correlación
##     entre dos sitios cae al ~5%. Dos pozos separados más de ~112 km
##     se consideran prácticamente independientes bajo este modelo.
##
##   τ² = 672.75 ft²: nugget. Representa la variabilidad a escala menor
##     que la distancia mínima entre pozos (~15 km) más el error de
##     medición de las sondas piezométricas. Es el ~15% del sill total
##     (672 / 4455 ≈ 0.15): moderado, indica que parte de la variabilidad
##     es ruido no estructurado, pero la mayor parte (~85%) es dependencia
##     espacial modelable.
##
##   μ (GLS) = 14.02 ft: media de los residuos estimada por GLS dentro de
##     likfit. Recordar que los residuos del lm() tienen media 0 por
##     construcción OLS, pero GLS pondera por Σθ⁻¹ y puede diferir
##     ligeramente. Este valor se usa como la media "conocida" en KS.
##
##   Grilla: 50 × 50 = 2.500 puntos uniformemente distribuidos sobre el
##     dominio x ∈ [−218, 178] km, y ∈ [−148, 139] km. Del gráfico se
##     aprecia que los pozos no cubren uniformemente el dominio: hay
##     clusters en el centro-este y zonas sin datos en los bordes oeste
##     y sur — esas regiones tendrán la mayor varianza de kriging (Parte d).

## Visualizamos el dominio: grilla + ubicaciones de los pozos
ggplot() +
  geom_point(data = grid_df, aes(x = x, y = y),
             size = 0.3, color = "grey80") +
  geom_point(data = df_wolf, aes(x = x, y = y),
             size = 2.5, color = "steelblue") +
  coord_equal() +
  theme_minimal() +
  labs(
    title    = "Grilla de predicción y ubicaciones observadas",
    subtitle = "Gris: 2.500 puntos de la grilla · Azul: 85 pozos de wolfcamp",
    x        = "X (km)",
    y        = "Y (km)"
  )

## La grilla cubre el bounding box del acuífero. Zonas alejadas de los pozos
## (bordes del dominio) tendrán mayor varianza de predicción — lo veremos
## en la Parte d.

# Parte b: Kriging Simple (KS) -------------------------------------------
## En Kriging Simple asumimos que la media del campo aleatorio es constante
## y CONOCIDA: E[Z(x)] = μ para todo x ∈ D.
##
## El predictor es Z*(x₀) = λ₀ + Σλᵢ Z(xᵢ), donde los pesos λᵢ se
## obtienen minimizando el ECM bajo la condición de insesgamiento:
##
##   E[Z*(x₀) − Z(x₀)] = 0  ⟹  λ₀ = μ(1 − Σλᵢ)
##
## Sustituyendo λ₀ en el predictor:
##   Z*(x₀) = μ + Σλᵢ [Z(xᵢ) − μ]
##
## y minimizando el ECM se obtiene:
##   λ = Σ⁻¹ k,  donde k = [C(x₀−x₁), …, C(x₀−xₙ)]ᵀ
##
## ECM mínimo: C(0) − kᵀ Σ⁻¹ k
##
## Usamos krige.conv de geoR con type.krige = "SK" y beta = μ conocida.

krige_sk <- krige.conv(
  geodata   = geo_res,
  locations = grid_mat,
  krige     = krige.control(
    type.krige = "SK",
    obj.model  = fit_ml,
    beta       = mu_ks      ## media conocida
  )
)

grid_df <- grid_df |>
  mutate(
    pred_sk = krige_sk$predict,
    var_sk  = krige_sk$krige.var
  )

## Mapa de predicción KS
ggplot(grid_df, aes(x = x, y = y, fill = pred_sk)) +
  geom_raster() +
  geom_point(data = data_res, aes(x = x, y = y),
             size = 2, color = "white", shape = 21, stroke = 0.5,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    title    = "Kriging Simple — Mapa de predicción",
    subtitle = paste0("Media conocida μ = ", round(mu_ks, 2), " ft · Modelo Exponencial (MV)"),
    x        = "X (km)", y = "Y (km)",
    fill     = "Residuo\npredicho (ft)"
  )

## El mapa muestra los residuos predichos por KS sobre el dominio del acuífero.
## El rango de colores va de ~−70 ft (púrpura, zonas de residuos negativos)
## a ~+150 ft (amarillo, zona de residuo positivo alto en el extremo oeste,
## x ≈ −220 km, y ≈ −10 km). Ese punto amarillo corresponde a un pozo con
## un residuo alto — la presión observada está muy por encima de lo que
## predice la tendencia lineal en esa zona, y el kriging lo captura localmente.
##
## En zonas alejadas de cualquier dato (bordes norte, sur y este de la grilla),
## el predictor KS converge hacia μ = 14.02 ft — la media conocida. Es el
## comportamiento esperado: cuando no hay información local, el mejor
## predictor es la media del proceso. Esto se ve como la transición gradual
## hacia tonos homogéneos (magenta/rosado) en los bordes.
##
## La estructura espacial es suave: el kriging produce superficies continuas
## que interpolan entre los datos (puntos blancos) respetando el rango de
## correlación φ = 37.28 km (rango práctico ~112 km). Pozos cercanos entre
## sí comparten colores similares — la correlación espacial está siendo
## incorporada explícitamente en la predicción.

## ── Verificación analítica de λ₀ en el centroide geográfico ────────────
## Construimos manualmente la matriz de covarianza Σ y el vector k para
## verificar la condición de insesgamiento λ₀ = μ(1 − Σλᵢ) en KS.

## Punto de predicción: centroide geográfico de la grilla
## Nota: round(nrow/2) NO da el centro espacial porque expand.grid varía
## x primero — la fila del medio tiene x = x_max, no x_centro.
## Usamos el punto de la grilla más cercano al centroide (x̄, ȳ).
x_centro   <- mean(range(x_seq))
y_centro   <- mean(range(y_seq))
idx_centro <- which.min(abs(grid_df$x - x_centro) + abs(grid_df$y - y_centro))
x0         <- grid_mat[idx_centro, , drop = FALSE]

## Función de covarianza Exponencial
cov_exp <- function(h, sigmasq, phi, tausq, nugget_en_cero = TRUE) {
  ## C(h) = sigmasq · exp(−h/phi) + τ² · I(h = 0)
  ## El nugget solo aparece en la diagonal de Σ (cuando h = 0 exacto)
  sigmasq * exp(-h / phi) + tausq * nugget_en_cero * (h == 0)
}

## Matriz Σ (85 × 85) y vector k (85 × 1)
coords <- geo_res$coords
n      <- nrow(coords)
D_obs  <- as.matrix(dist(coords))
Sigma  <- cov_exp(D_obs, sigmasq, phi, tausq, nugget_en_cero = TRUE)

d_x0  <- sqrt((coords[, 1] - x0[1, 1])^2 + (coords[, 2] - x0[1, 2])^2)
k_x0  <- cov_exp(d_x0, sigmasq, phi, tausq, nugget_en_cero = FALSE)

## λ = Σ⁻¹ k
lambda_ks <- solve(Sigma, k_x0)

## λ₀ = μ(1 − Σλᵢ)
lambda0_ks <- mu_ks * (1 - sum(lambda_ks))

cat("=== Verificación analítica λ₀ (KS, punto central de la grilla) ===\n")
cat("x₀            = (", round(x0[1], 1), ",", round(x0[2], 1), ") km\n")
cat("Σλᵢ           =", round(sum(lambda_ks), 6), "\n")
cat("1 − Σλᵢ       =", round(1 - sum(lambda_ks), 6), "\n")
cat("μ             =", round(mu_ks, 4), "ft\n")
cat("λ₀ = μ(1−Σλᵢ) =", round(lambda0_ks, 4), "ft\n")
cat("Predictor KS  =", round(lambda0_ks + sum(lambda_ks * geo_res$data), 4), "ft\n")

## La verificación se realizó en el centroide geográfico x₀ = (−30.2, −7.4) km,
## punto interior al dominio bien rodeado de pozos.
##
## Resultados:
##   Σλᵢ = 0.9399: muy cerca de 1. El centroide está bien cubierto por los
##     85 pozos — los pesos absorben casi toda la variabilidad local y los
##     datos hacen el 94% del trabajo en la predicción. Contrastar con el
##     borde este (x₀ ≈ 182 km, antes del fix): allí Σλᵢ = 0.671,
##     reflejando la mayor escasez de información en la periferia.
##
##   1 − Σλᵢ = 0.060: fracción "no explicada" por los datos. Es pequeña
##     porque x₀ tiene pozos relativamente cercanos en todas las direcciones.
##
##   λ₀ = μ · (1 − Σλᵢ) = 14.02 × 0.060 = 0.843 ft: contribución de la
##     media global al predictor. Es mínima — el predictor casi no necesita
##     "recordar" μ porque los datos locales ya son suficientes. En zonas
##     de borde (sin datos), λ₀ crece y el predictor va convergiendo a μ.
##
##   Predictor KS = λ₀ + Σλᵢ·Z(xᵢ) = 0.843 + 50.317 ≈ 51.16 ft:
##     residuo predicho en el centroide. Al sumarlo a la tendencia lineal
##     evaluada en (−30.2, −7.4) se obtendría la presión piezométrica final.

# Parte c: Kriging Ordinario (KO) ----------------------------------------
## En Kriging Ordinario la media es constante pero DESCONOCIDA.
## La condición de insesgamiento λ₀ + μ(Σλᵢ − 1) = 0 debe cumplirse
## para cualquier valor de μ, lo que exige:
##
##   λ₀ = 0  y  Σλᵢ = 1
##
## El predictor queda Z*(x₀) = Σλᵢ Z(xᵢ), y el problema de minimización
## del ECM sujeto a Σλᵢ = 1 se resuelve con multiplicadores de Lagrange.
##
## La función de Lagrange es:
##   L(λ, η) = λᵀΣλ − 2λᵀk + 2η(Σλᵢ − 1)
##
## Derivando e igualando a cero:
##   Σλ + η1 = k   (condición de optimalidad)
##   1ᵀλ = 1       (restricción de insesgamiento)
##
## Esto conduce al sistema lineal aumentado de tamaño (n+1) × (n+1):
##   [Σ  1][λ]   [k]
##   [1ᵀ 0][η] = [1]
##
## ECM mínimo: C(0) − λᵀk + η
##
## krige.conv con type.krige = "OK" resuelve este sistema internamente.

krige_ok <- krige.conv(
  geodata   = geo_res,
  locations = grid_mat,
  krige     = krige.control(
    type.krige = "OK",
    obj.model  = fit_ml
  )
)

grid_df <- grid_df |>
  mutate(
    pred_ok = krige_ok$predict,
    var_ok  = krige_ok$krige.var
  )

## ── Verificación manual del sistema de Lagrange ─────────────────────────
## Resolvemos el sistema aumentado para el mismo punto central x₀

ones <- rep(1, n)
A_ok <- rbind(cbind(Sigma, ones), c(ones, 0))   ## (n+1) × (n+1)
b_ok <- c(k_x0, 1)
sol  <- solve(A_ok, b_ok)

lambda_ok <- sol[1:n]       ## pesos KO
eta_ok    <- sol[n + 1]     ## multiplicador de Lagrange

cat("\n=== Verificación Sistema Kriging Ordinario (Lagrange) ===\n")
cat("Σλᵢ            =", round(sum(lambda_ok), 8), "(debe ser exactamente 1)\n")
cat("η (Lagrange)   =", round(eta_ok, 6), "\n")
cat("Predictor KO   =", round(sum(lambda_ok * geo_res$data), 4), "ft\n")
cat("ECM = C(0)−λᵀk+η =",
    round(cov_exp(0, sigmasq, phi, tausq) - sum(lambda_ok * k_x0) + eta_ok, 4), "ft²\n")

## Verificación del sistema de Lagrange en x₀ = (−30.2, −7.4) km:
##
##   Σλᵢ = 1 exacto: la restricción de insesgamiento se satisface a
##     precisión de máquina — el sistema lineal aumentado está bien
##     planteado y solve() lo resuelve sin problemas.
##
##   η = −15.77: el multiplicador de Lagrange puede ser negativo — no hay
##     ningún problema. Su signo depende de la geometría local: η < 0
##     indica que la restricción Σλᵢ = 1 "ayuda" a reducir el ECM en
##     este punto (los pesos sin restricción sumarían más de 1, y la
##     restricción los "comprime").
##
##   Predictor KO = 51.1595 ft = Predictor KS: son idénticos. Esto no
##     es una coincidencia ni un error — es un resultado teórico. Cuando
##     la media "conocida" usada en KS es el estimador GLS de μ (que es
##     exactamente lo que likfit devuelve en fit_ml$beta), KS y KO son
##     predictores equivalentes. El estimador GLS es el mejor estimador
##     lineal insesgado de μ bajo correlación espacial, y es justamente
##     lo que KO estima implícitamente vía la restricción Σλᵢ = 1.
##
##   ECM = 2.512 ft² × 10³ = 2511.97 ft²: varianza de predicción de KO
##     en el centroide. C(0) = σ² + τ² = 3782 + 673 = 4455 ft²,
##     y el kriging reduce esa incertidumbre un ~44% gracias a los datos
##     cercanos. El error estándar de predicción es √2511.97 ≈ 50.1 ft,
##     que es sustancial pero esperable para un punto sin dato propio.
##     La varianza KO ≥ varianza KS (propiedad general) — KO paga el
##     costo de no conocer μ con mayor incertidumbre de predicción.

## ── Mapas comparativos KS vs KO ─────────────────────────────────────────
lim_pred <- range(c(grid_df$pred_sk, grid_df$pred_ok))

p_ks <- ggplot(grid_df, aes(x = x, y = y, fill = pred_sk)) +
  geom_raster() +
  geom_point(data = data_res, aes(x = x, y = y),
             size = 1.5, color = "white", shape = 21, stroke = 0.4,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "plasma", limits = lim_pred) +
  coord_equal() + theme_minimal() +
  labs(title = "Kriging Simple",
       subtitle = "Media conocida μ",
       x = "X (km)", y = "Y (km)", fill = "Residuo (ft)")

p_ko <- ggplot(grid_df, aes(x = x, y = y, fill = pred_ok)) +
  geom_raster() +
  geom_point(data = data_res, aes(x = x, y = y),
             size = 1.5, color = "white", shape = 21, stroke = 0.4,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "plasma", limits = lim_pred) +
  coord_equal() + theme_minimal() +
  labs(title = "Kriging Ordinario",
       subtitle = "Media desconocida, Σλᵢ = 1",
       x = "X (km)", y = "Y (km)", fill = "Residuo (ft)")

p_ks + p_ko +
  plot_annotation(title = "Comparación KS vs KO — Residuos de tendencia lineal")

## ── Mapa de diferencias KO − KS ─────────────────────────────────────────
grid_df <- grid_df |> mutate(dif = pred_ok - pred_sk)

ggplot(grid_df, aes(x = x, y = y, fill = dif)) +
  geom_raster() +
  geom_point(data = data_res, aes(x = x, y = y),
             size = 1.5, color = "grey20", shape = 21, stroke = 0.4,
             inherit.aes = FALSE) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato",
                       midpoint = 0) +
  coord_equal() + theme_minimal() +
  labs(
    title    = "Diferencia KO − KS",
    subtitle = "Las mayores discrepancias aparecen en zonas sin datos (bordes del dominio)",
    x        = "X (km)", y = "Y (km)",
    fill     = "KO − KS (ft)"
  )

## Los dos mapas son visualmente indistinguibles, lo que es consistente
## con los predictores numéricos idénticos (51.16 ft en el centroide).
##
## El mapa de diferencias KO − KS lo confirma: el fondo es prácticamente
## blanco en todo el dominio — las diferencias entre ambos predictores
## son cercanas a cero en todos los 2.500 puntos de la grilla. Esto
## ocurre porque μ_KS = β̂_GLS (el estimador GLS de likfit), que es
## exactamente lo que KO estima internamente. Cuando la media "conocida"
## de KS es el estimador óptimo de μ bajo correlación espacial, ambos
## métodos convergen a la misma predicción.
##
## Si hubiéramos usado una media incorrecta en KS (por ejemplo μ = 0 o
## μ = media aritmética de los residuos), el mapa de diferencias mostraría
## gradientes de color intensos en los bordes, donde KS extrapolaría hacia
## esa μ errónea mientras que KO se mantendría anclado a los datos.

# Parte d: Varianza de predicción e interpolador exacto ------------------
## ── Mapas de varianza ───────────────────────────────────────────────────
## La varianza de kriging σ²ₖ(x₀) = C(0) − kᵀΣ⁻¹k (KS) o C(0)−λᵀk+η (KO)
## cuantifica la incertidumbre del predictor en cada punto de la grilla.
## Es una función solo del diseño de muestreo y del modelo de covarianza
## — no depende de los valores observados. Por eso el mapa de varianza
## tiene la misma forma independientemente de la realización observada.
##
## La varianza es mínima cerca de los pozos (información directa) y
## máxima en los bordes del dominio (lejos de cualquier dato).

lim_var <- range(c(grid_df$var_sk, grid_df$var_ok))

p_var_ks <- ggplot(grid_df, aes(x = x, y = y, fill = var_sk)) +
  geom_raster() +
  geom_point(data = data_res, aes(x = x, y = y),
             size = 1.5, color = "white", shape = 21, stroke = 0.4,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "magma", limits = lim_var) +
  coord_equal() + theme_minimal() +
  labs(title    = "Varianza de predicción — KS",
       x = "X (km)", y = "Y (km)", fill = "Var (ft²)")

p_var_ko <- ggplot(grid_df, aes(x = x, y = y, fill = var_ok)) +
  geom_raster() +
  geom_point(data = data_res, aes(x = x, y = y),
             size = 1.5, color = "white", shape = 21, stroke = 0.4,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "magma", limits = lim_var) +
  coord_equal() + theme_minimal() +
  labs(title    = "Varianza de predicción — KO",
       x = "X (km)", y = "Y (km)", fill = "Var (ft²)")

p_var_ks + p_var_ko +
  plot_annotation(title = "Mapas de varianza de kriging — KS y KO")

## Los dos mapas muestran el mismo patrón espacial: varianza mínima
## (negro/púrpura oscuro) en las ubicaciones de los pozos y varianza
## máxima (amarillo) en los bordes del dominio sin datos, especialmente
## en el extremo oeste (x < −150 km) donde hay muy poca cobertura.
##
## La escala va de ~1.700 a ~4.400 ft². El mínimo ~1.700 ft² aparece
## en los clusters de pozos del centro-este — allí la información local
## reduce la incertidumbre en un ~62% respecto de C(0) = 4.455 ft².
## El máximo ~4.400 ft² ≈ C(0) en los bordes: sin datos cercanos,
## el kriging no puede reducir la incertidumbre inicial.
##
## La varianza de kriging depende únicamente del diseño de muestreo
## (posición de los pozos) y del modelo de covarianza — NO de los
## valores observados. Por eso el mapa de varianza tiene la misma
## forma independientemente de la realización observada de la presión.
##
## Comparando KS y KO: ambos mapas son casi idénticos (la escala es
## compartida y los patrones coinciden). La diferencia teórica
## var_KO ≥ var_KS existe pero es pequeña en este caso porque μ_KS
## es el estimador GLS óptimo — KO paga un costo mínimo por no conocer μ.

## ── Interpolador exacto (propiedad teórica, nugget = 0) ─────────────────
## El kriging es teóricamente un INTERPOLADOR EXACTO cuando τ² = 0:
## si predecimos en un sitio xᵢ donde tenemos un dato Z(xᵢ), el predictor
## devuelve Z*(xᵢ) = Z(xᵢ) y la varianza de kriging es σ²ₖ(xᵢ) = 0.
##
## ¿Por qué? Cuando x₀ = xᵢ y τ² = 0, el vector k tiene en la posición i
## el valor C(0) = σ², y el peso λᵢ = 1 (los demás → 0). El predictor
## colapsa al dato observado y el ECM C(0) − kᵀΣ⁻¹k = 0.
##
## Con τ² > 0 (nugget), el kriging es un smoother: predice un valor
## intermedio entre la media y el dato, y la varianza en datos es τ² > 0.
## Esto modela el caso en que parte de la variabilidad es error de medición.
##
## Demostramos empíricamente con τ² = 0:

krige_exacto <- krige.conv(
  geodata   = geo_res,
  locations = geo_res$coords,   ## predicción en las coordenadas observadas
  krige     = krige.control(
    type.krige = "OK",
    cov.model  = "exponential",
    cov.pars   = c(sigmasq, phi),
    nugget     = 0              ## τ² = 0 para demostrar la propiedad exacta
  )
)

exact_df <- data.frame(
  observado = geo_res$data,
  predicho  = krige_exacto$predict,
  varianza  = krige_exacto$krige.var
)

## Verificación puntual en el sitio i = 1
cat("\n=== Interpolador Exacto (τ² = 0, Kriging Ordinario) ===\n")
cat("Sitio i = 1:\n")
cat("  Z(x₁) observado  =", round(exact_df$observado[1], 4), "ft\n")
cat("  Z*(x₁) predicho  =", round(exact_df$predicho[1], 4), "ft\n")
cat("  Error |Z*−Z|     =", round(abs(exact_df$predicho[1] - exact_df$observado[1]), 6), "ft\n")
cat("  Varianza σ²ₖ(x₁) =", round(exact_df$varianza[1], 6), "ft²\n")
cat("\nResumen sobre los 85 sitios:\n")
cat("  Error máximo =", round(max(abs(exact_df$predicho - exact_df$observado)), 6), "ft\n")
cat("  Varianza máx =", round(max(exact_df$varianza), 6), "ft²\n")

## Propiedad demostrada sobre los 85 sitios observados (τ² = 0):
##
##   Z(x₁) = −22.9078 ft → Z*(x₁) = −22.9078 ft: el predictor
##     reproduce exactamente el valor observado en el primer pozo.
##     Error = 0 ft y varianza = 0 ft² a precisión de máquina.
##
##   Error máximo = 0 ft y Varianza máxima = 0 ft² sobre los 85 sitios:
##     la propiedad se cumple en TODOS los pozos, no solo en el primero.
##     El scatter plot observado vs predicho muestra los 85 puntos
##     perfectamente alineados sobre la diagonal y = x, sin ninguna
##     dispersión.
##
## ¿Por qué ocurre esto algebraicamente? Cuando x₀ = xᵢ y τ² = 0:
##   - El vector k tiene kᵢ = C(0) = σ² y k_j = C(xᵢ − x_j) para j ≠ i
##   - k es exactamente la i-ésima columna de Σ (ya que Σᵢⱼ = C(xᵢ − x_j))
##   - Por lo tanto Σ⁻¹k = eᵢ (vector canónico): λᵢ = 1, λⱼ = 0 ∀j ≠ i
##   - Predictor = Σλⱼ Z(xⱼ) = Z(xᵢ) ✓
##   - Varianza = C(0) − kᵀ Σ⁻¹k = C(0) − kᵀeᵢ = C(0) − C(0) = 0 ✓
##
## Con τ² > 0 (nuestro modelo real, τ² = 672.75 ft²), el kriging deja
## de ser interpolador exacto: la varianza en datos es τ² > 0 y el
## predictor "suaviza" hacia la media en lugar de reproducir el dato.
## Esto modela explícitamente el error de medición de las sondas.

## Gráfico de dispersión: observado vs predicho (nugget = 0)
ggplot(exact_df, aes(x = observado, y = predicho)) +
  geom_point(size = 2.5, color = "steelblue", alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "tomato", linewidth = 1) +
  theme_minimal() +
  labs(
    title    = "Interpolador exacto — KO con τ² = 0",
    subtitle = "Los puntos deben caer exactamente sobre la diagonal",
    x        = "Z(xᵢ) observado (ft)",
    y        = "Z*(xᵢ) predicho (ft)"
  )

# Parte e: Costo computacional O(n³) -------------------------------------
## El paso central del kriging es resolver el sistema Σλ = k (KS) o el
## sistema aumentado de tamaño (n+1)×(n+1) (KO). Ambos requieren la
## factorización de Σ, que en la práctica se realiza por descomposición
## de Cholesky: Σ = LLᵀ, con L triangular inferior.
##
## El costo de la factorización de Cholesky es O(n³/3) ≈ O(n³):
##   - Para cada columna j = 1, …, n se realizan O(n²) operaciones
##   - En total: n × O(n²) = O(n³)
##
## Una vez que L está disponible, resolver Lx = b (o Lᵀx = b) cuesta
## O(n²) por sustitución hacia adelante/atrás. Pero la factorización
## domina para n grande.
##
## Consecuencias prácticas:
##   n = 85   (wolfcamp):  Σ es 85×85  → ~614.000 operaciones → instantáneo
##   n = 1.000:            Σ es 1000×1000 → ~3×10⁸ operaciones → ~0.1 s
##   n = 10.000:           Σ es 10⁴×10⁴ → ~3×10¹¹ operaciones → minutos
##   n = 100.000:          Σ es 10⁵×10⁵ → ~3×10¹⁴ operaciones → días
##   n = 1.000.000:        Σ es 10⁶×10⁶ → inviable en memoria (8 TB solo para Σ)
##
## Para n grande se usan aproximaciones:
##   - Vecinos más cercanos (NNGP): solo se usa un subconjunto de n' ≪ n vecinos
##   - Métodos low-rank: Σ ≈ UUᵀ + τ²I con U de rango r ≪ n
##   - Lattice approximations (SPDE/INLA): aprovechan estructura de grilla
##   - Tapering: se fuerza C(h) = 0 para h > umbral → Σ es sparse

## Demostración empírica del crecimiento O(n³):
## Medimos el tiempo de solve() (= Cholesky + sustitución) para distintos n.
ns      <- c(50, 85, 200, 500, 1000)
tiempos <- sapply(ns, function(n_sim) {
  set.seed(42)
  ## Matriz simétrica definida positiva de tamaño n_sim × n_sim
  A <- matrix(rnorm(n_sim^2), n_sim, n_sim)
  A <- A %*% t(A) + diag(n_sim)
  b <- rnorm(n_sim)
  system.time(solve(A, b))["elapsed"]
})

timing_df <- data.frame(n = ns, tiempo_s = pmax(tiempos, 1e-4))

ggplot(timing_df, aes(x = n, y = tiempo_s)) +
  geom_point(size = 3, color = "steelblue") +
  geom_line(color = "steelblue", linewidth = 0.9) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10(labels = scales::scientific) +
  theme_minimal() +
  labs(
    title    = "Costo computacional de la inversión matricial",
    subtitle = "Escala log-log: pendiente ≈ 3 confirma O(n³)",
    x        = "Tamaño n (escala log₁₀)",
    y        = "Tiempo de cómputo en segundos (escala log₁₀)"
  )

## El gráfico tiene dos tramos bien diferenciados:
##
## Tramo plano (n = 50, 85, 200 → tiempo = 1e-4 s): es un artefacto de
##   medición. system.time() en Windows tiene resolución de ~10 ms —
##   operaciones que terminan en < 1 ms devuelven "elapsed = 0", que
##   pmax(..., 1e-4) eleva al piso de 1e-4. El kriging con n ≤ 200 es
##   prácticamente instantáneo en cualquier hardware moderno, incluyendo
##   nuestros 85 pozos de Wolfcamp.
##
## Tramo creciente (n = 500 → n = 1000):
##   - n = 500:  ~0.03 s
##   - n = 1000: ~0.20 s
##   - Razón: 0.20 / 0.03 ≈ 6.7 (esperado teórico: 2³ = 8)
##   El valor observado (6.7) está por debajo de 8 debido a ruido de
##   medición, efectos de caché y las optimizaciones BLAS de R (LAPACK
##   usa DGEMM vectorizado). Aun así, la pendiente en escala log-log es
##   consistente con O(n³).
##
## Extrapolando desde n = 1000 (t ≈ 0.2 s) con factor 8 por cada
## duplicación de n:
##   n = 2.000:  ~1.6 s
##   n = 10.000: ~0.2 × (10)³ ≈ 200 s  (~3 minutos)
##   n = 50.000: ~0.2 × (50)³ ≈ 250.000 s (~3 días)
##
## Esto confirma que el kriging exacto es inviable para datasets
## modernos de teledetección o sensores continuos (n > 10⁴), y motiva
## el uso de las aproximaciones mencionadas arriba (NNGP, low-rank,
## SPDE/INLA, tapering).
