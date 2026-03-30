## Ay_03_Espacial.R -------------------------------------------------------
## EYP3417 — Estadística Espacial
## Ayudantía 03: Inferencia Paramétrica — Acuífero Wolfcamp
## Profesor: Alfredo Alegría  |  Ayudante: Juan Pino

# 0. Setup ----------------------------------------------------------------
library(geoR)
library(gstat)
library(tidyverse)
library(patchwork)

## Retomamos el dataset wolfcamp de la sesión anterior. El objeto wolfcamp
## es de clase geodata, por lo que extraemos coordenadas y datos directamente.
data(wolfcamp)

df_wolf <- data.frame(
  x        = wolfcamp$coords[, 1],
  y        = wolfcamp$coords[, 2],
  pressure = wolfcamp$data
)

# Parte a: Tendencia, residuos y variograma experimental -----------------
## En la ayudantía anterior verificamos que la carga piezométrica presenta
## una tendencia espacial clara (gradiente SW → NE, R² = 0.89). Para
## estimar formalmente los parámetros de covarianza θ = (σ², φ, τ²)
## necesitamos trabajar sobre un proceso con media constante.
##
## ¿Por qué es fundamental? El variograma y la función de covarianza están
## definidos bajo el supuesto de estacionariedad: γ(h) depende solo de la
## separación h, no de la ubicación s. Si Z(s) = μ(s) + ε(s) con μ(s)
## no constante, entonces γ_Z(h) ya no es función solo de h y los
## estimadores de θ estarían sesgados.
##
## Al ajustar un modelo de tendencia Z(s) = β₀ + β₁X + β₂Y + ε(s) y
## trabajar con los residuos ε̂(s), removemos la componente determinista
## y recuperamos un proceso con E[ε̂(s)] ≈ 0, cuya estructura de covarianza
## sí puede estimarse consistentemente.
##
## Desde el punto de vista de la teoría vista en clases: si asumimos
## Z ~ Nn(Xβ, Σθ), el estimador de β es β̂ = (XᵀΣθ⁻¹X)⁻¹XᵀΣθ⁻¹z,
## que se estima junto con θ maximizando la verosimilitud perfilada.
## En la práctica con dos pasos (lm + variograma residuos), separamos
## la estimación de la tendencia de la de la covarianza.

model_trend <- lm(pressure ~ x + y, data = df_wolf)
summary(model_trend)

## Extraemos los residuos y los añadimos al dataframe
data_res <- df_wolf |>
  mutate(residuals = residuals(model_trend))

## Calculamos el variograma experimental de los residuos
v_res <- variogram(residuals ~ 1, locations = ~ x + y,
                   data = data_res, cutoff = 150, width = 30)

ggplot(v_res, aes(x = dist, y = gamma)) +
  geom_point(size = 3, color = "grey30") +
  geom_line(linewidth = 0.8, color = "grey30") +
  theme_minimal() +
  labs(
    title    = "Variograma Experimental — Residuos de Tendencia Lineal",
    subtitle = "Base para la estimación de θ = (σ², φ, τ²)",
    x        = "Distancia (km)",
    y        = expression(hat(gamma)(h))
  )

## El modelo de tendencia ajustado es Ẑ(s) = 607.77 − 1.278·X − 1.139·Y,
## idéntico al obtenido en la Ayudantía 02. Ambos coeficientes son
## negativos y altamente significativos (p < 2e-16): la presión disminuye
## ~1.28 ft por km hacia el Este y ~1.14 ft por km hacia el Norte. El
## R² = 0.891 confirma que la tendencia lineal domina la variabilidad total.
##
## Del variograma de los residuos se aprecia que la curva crece desde
## ~1.900 ft² a distancias cortas (~15 km) y se estabiliza en torno a
## ~4.000–4.300 ft² alrededor de los 100–110 km, con una leve caída
## posterior que indica que el sill se ha alcanzado. Este comportamiento
## acotado confirma que los residuos ε̂(s) son aproximadamente estacionarios
## de segundo orden: la remoción de la tendencia lineal recuperó la
## condición necesaria para estimar θ = (σ², φ, τ²) de forma consistente.
##
## La curvatura del variograma cerca del origen (crecimiento rápido hasta
## ~50 km y luego desaceleración) sugiere un nugget moderado y un range
## en torno a los 100 km. Estos valores visuales servirán como valores
## iniciales de referencia para los ajustes de las Partes b, c y d.

# Parte b: MCO vs MCP — fit.variogram ------------------------------------
## El ajuste de un modelo de variograma consiste en minimizar el error
## cuadrático ponderado entre el variograma empírico y el teórico:
##
##   Σ wk · [γ̂(hk) − γθ(hk)]²
##
## La diferencia entre MCO y MCP está únicamente en los pesos wk:
##
##   MCO (fit.method = 6): wk = 1 para todo k. Todos los lags tienen
##     el mismo peso sin importar cuántos pares los forman ni su distancia.
##
##   MCP (fit.method = 7): wk ∝ Nk / [γθ(hk)]², donde Nk es el número
##     de pares en el lag k. Los lags con más pares y menor semivarianza
##     reciben mayor peso, favoreciendo el ajuste cerca del origen.
##
## Esto importa porque los lags cortos tienen más pares (estimaciones
## más estables) y son los más relevantes para kriging. Con MCO, un lag
## lejano con pocos pares y alta varianza estimada influye tanto como
## uno cercano bien estimado, lo que puede distorsionar el nugget y la
## curvatura cerca del origen.

## Al igual que en la Ayudantía 02, el ajuste puede fallar por convergencia
## si los valores iniciales no son adecuados. Reutilizamos la función de
## búsqueda en grilla que prueban 48 combinaciones y retiene la de menor SSR.
## Ahora el argumento fit.method nos permite aplicarla tanto a MCO como MCP.

try_fit_variogram <- function(v_emp, model,
                              fit.method  = 7,
                              psill_vals  = c(500, 1000, 2000, 3000),
                              range_vals  = c(30,  60,  100,  150),
                              nugget_vals = c(0,   200,  500)) {
  grid <- expand.grid(psill = psill_vals, range = range_vals, nugget = nugget_vals)

  best_ssr <- Inf
  best_fit <- NULL

  for (i in seq_len(nrow(grid))) {
    fit <- tryCatch(
      withCallingHandlers(
        fit.variogram(v_emp,
                      model      = vgm(psill  = grid$psill[i],
                                       model  = model,
                                       range  = grid$range[i],
                                       nugget = grid$nugget[i]),
                      fit.method = fit.method),
        warning = \(w) invokeRestart("muffleWarning")
      ),
      error = \(e) NULL
    )
    if (!is.null(fit)) {
      ssr <- attr(fit, "SSErr")
      if (!is.na(ssr) && ssr < best_ssr) {
        best_ssr <- ssr
        best_fit <- fit
      }
    }
  }
  best_fit
}

## Ajuste por MCO (Mínimos Cuadrados Ordinarios, fit.method = 6)
fit_ols <- try_fit_variogram(v_res, model = "Exp", fit.method = 6)
fit_ols

## Ajuste por MCP (Mínimos Cuadrados Ponderados, fit.method = 7)
fit_wls <- try_fit_variogram(v_res, model = "Exp", fit.method = 7)
fit_wls

## Visualizamos ambos ajustes superpuestos al variograma empírico
h_seq <- seq(0, max(v_res$dist) * 1.05, length.out = 300)

curvas_df <- bind_rows(
  tibble(dist   = h_seq,
         gamma  = variogramLine(fit_ols, dist_vector = h_seq)$gamma,
         metodo = "MCO (fit.method = 6)"),
  tibble(dist   = h_seq,
         gamma  = variogramLine(fit_wls, dist_vector = h_seq)$gamma,
         metodo = "MCP (fit.method = 7)")
)

ggplot() +
  geom_point(data = v_res, aes(x = dist, y = gamma),
             size = 3, color = "grey30") +
  geom_line(data = curvas_df,
            aes(x = dist, y = gamma, color = metodo),
            linewidth = 1.2) +
  theme_minimal() +
  labs(
    title    = "Comparación MCO vs MCP — Modelo Exponencial",
    subtitle = "Residuos de tendencia lineal, width = 30 km",
    x        = "Distancia (km)",
    y        = expression(hat(gamma)(h)),
    color    = NULL
  )

## Comparamos los parámetros estimados por cada método
data.frame(
  Método = c("MCO (OLS)", "MCP (WLS)"),
  Nugget = c(fit_ols$psill[1], fit_wls$psill[1]),
  PSill  = c(fit_ols$psill[2], fit_wls$psill[2]),
  Range  = c(fit_ols$range[2], fit_wls$range[2])
)

## Del gráfico se aprecian dos comportamientos bien diferenciados:
##
## MCO (rojo): nugget bajo (~722 ft²), curva con crecimiento más rápido
##   en los primeros 60 km y que tiende a aplanarse antes. Al asignar igual
##   peso a todos los lags, los puntos lejanos —donde la curva aún sigue
##   subiendo— "frenan" el ajuste del range, resultando en un range corto
##   (~59 km) que no acomoda bien los puntos intermedios (~80 km, ~100 km).
##
## MCP (cyan): nugget más alto (~1228 ft²), curva que crece más lentamente
##   pero de forma sostenida hasta los 130 km. Al ponderar los lags cortos
##   (más pares, más confiables), el ajuste privilegia la región del origen
##   y el range se estima mayor (~130 km), más consistente con la forma
##   del variograma empírico visto en la Parte a.
##
## Visualmente, ambas curvas pasan cerca del primer punto (~15 km, ~1900 ft²),
## pero divergen a partir de ~50 km: MCP sigue la tendencia creciente de los
## puntos intermedios mejor que MCO.
##
## Las diferencias entre ambos métodos son sustanciales:
##
##              Nugget    PSill   Range
##   MCO (OLS)   722 ft²  4034 ft²   59 km
##   MCP (WLS)  1228 ft²  5074 ft²  130 km
##
## La diferencia más llamativa es el range: MCO estima ~59 km, MCP ~130 km.
## Esto se explica por cómo cada método trata los lags lejanos:
##
## MCO asigna igual peso a todos los lags. Los lags lejanos tienen pocos
## pares y estimaciones ruidosas, pero influyen tanto como los cercanos.
## En este caso esos puntos lejanos "tiran" la curva hacia abajo y el
## ajuste converge con un range corto para acomodar la nube de puntos
## completa sin discriminar su calidad.
##
## MCP pondera por Nk / [γθ(hk)]², favoreciendo los lags cortos (más pares,
## más confiables). El nugget sube (~1228 ft²) para ajustar bien la región
## del origen, y el range resulta mayor (~130 km), más consistente con lo
## que observamos en el variograma empírico de la Parte a, donde el sill
## se alcanzaba cerca de los 100–110 km.
##
## En la práctica, MCP es el método por defecto en gstat y el preferido
## cuando se realizará kriging, pues la calidad de la predicción depende
## principalmente del ajuste a corta distancia. MCO (fit.method = 6) es
## el método que solicita explícitamente la Tarea 1.

# Parte c: Estimación por Máxima Verosimilitud — likfit ------------------
## Pasamos ahora al método de máxima verosimilitud. Para un campo aleatorio
## gaussiano Z ~ Nn(μθ, Σθ), la log-verosimilitud es:
##
##   ℓ(θ|Z) = −n/2 log(2π) − 1/2 log|Σθ| − 1/2 (Z−μθ)ᵀ Σθ⁻¹ (Z−μθ)
##
## y el EMV es θ̂ = argmax ℓ(θ|Z), que se resuelve numéricamente con
## un costo computacional de orden O(n³) — la inversión de Σθ.
##
## Ventajas de MV sobre mínimos cuadrados:
##   - Cuando el proceso es efectivamente gaussiano, el EMV es asintótica-
##     mente eficiente (mínima varianza entre estimadores consistentes).
##   - Permite comparación formal de modelos vía AIC y BIC.
##   - Estima todos los parámetros (incluyendo la media) simultáneamente.
##
## Desventaja: requiere el supuesto distribucional gaussiano y es mucho
## más costoso computacionalmente (O(n³) vs O(m) de mínimos cuadrados,
## donde m es el número de lags).
##
## Usaremos likfit de geoR, que requiere un objeto de clase geodata.

## Convertimos data_res a objeto geodata
geo_res <- as.geodata(
  data_res,
  coords.col = c(1, 2),   ## columnas x, y
  data.col   = 4          ## columna residuals
)

## Estimación por MV con modelo exponencial
## ini.cov.pars = c(sigmasq, phi): valores iniciales para la optimización
fit_ml <- likfit(
  geodata      = geo_res,
  cov.model    = "exponential",
  ini.cov.pars = c(2000, 100),
  nugget       = 500,
  lik.method   = "ML"
)

summary(fit_ml)

## Del summary podemos leer los parámetros estimados:
##
##   beta (μ):     14.02 ft. Notemos que pasamos residuos del lm(), que
##                 tienen media 0 por construcción (OLS). Sin embargo,
##                 likfit estima β internamente usando GLS, que minimiza
##                 (Z − μ)ᵀ Σθ⁻¹ (Z − μ) en lugar de la suma de cuadrados
##                 ordinaria. Cuando hay correlación espacial, el estimador
##                 GLS puede diferir levemente del OLS, de ahí que β̂ ≠ 0.
##
##   sigmasq (σ²): 3782 ft² — varianza parcial del proceso.
##
##   phi (φ):      37.28 km — parámetro de escala. Para el Exponencial,
##                 el rango práctico es 3φ ≈ 111.7 km (distancia a la
##                 que la correlación cae al ~5%). Este valor está reportado
##                 explícitamente en el summary como "Practical Range".
##
##   nugget (τ²):  672.8 ft² — variabilidad a escala menor que la distancia
##                 mínima entre pozos o error de medición.
##
##   AIC:          926 (modelo espacial) vs 944.6 (modelo no espacial).
##                 La diferencia de ~18.6 puntos de AIC confirma que
##                 incorporar la estructura de covarianza espacial mejora
##                 sustancialmente el ajuste del modelo.
##
## Comparamos los parámetros de los tres métodos:
data.frame(
  Método         = c("MCO", "MCP", "MV (likfit)"),
  Nugget         = c(fit_ols$psill[1], fit_wls$psill[1], fit_ml$tausq),
  PSill          = c(fit_ols$psill[2], fit_wls$psill[2], fit_ml$sigmasq),
  Range_practico = c(fit_ols$range[2], fit_wls$range[2], fit_ml$phi * 3)
)
## Nota: para el Exponencial en geoR, rango práctico = 3 × phi.
##
## Podemos ver que MV da el nugget más bajo (672.8), similar a MCO (722),
## y un rango práctico de ~111.7 km, intermedio entre MCO (~59 km) y MCP
## (~130 km). Esto ilustra la diferencia de eficiencia estadística entre
## los métodos: MV usa toda la información de la distribución conjunta
## (la verosimilitud), mientras que MCO y MCP solo explotan el variograma
## empírico (un resumen de segundo orden). Cuando el supuesto gaussiano
## se cumple, MV es el estimador asintóticamente más eficiente.

# Parte d: Familia Matérn — comparación por AIC --------------------------
## El modelo de Matérn es una familia paramétrica de covarianzas controlada
## por un parámetro de suavidad ν ≥ 0, que determina la diferenciabilidad
## media cuadrática del campo aleatorio:
##
##   C(h; σ², φ, ν) = σ² / (2^(ν−1) Γ(ν)) · (h/φ)^ν · K_ν(h/φ)
##
## donde K_ν es la función de Bessel modificada de segundo tipo.
## La interpretación de ν es directa: un campo aleatorio Z(s) es k veces
## diferenciable en media cuadrática si y solo si ν > k. Así:
##   ν = 0.5 → equivalente al Exponencial, proceso no diferenciable
##             (trayectorias continuas pero con cambios abruptos)
##   ν = 1.5 → proceso una vez diferenciable en media cuadrática
##             (trayectorias más suaves y regulares)
##   ν = 2.5 → dos veces diferenciable
##   ν → ∞  → converge al Gaussiano, infinitamente diferenciable
##
## En geoR, ν se controla con el argumento kappa en likfit.

## Matérn con ν = 0.5 (idéntico al Exponencial)
fit_mat05 <- likfit(
  geodata      = geo_res,
  cov.model    = "matern",
  kappa        = 0.5,
  ini.cov.pars = c(2000, 100),
  nugget       = 500,
  lik.method   = "ML"
)

## Matérn con ν = 1.5 (una vez diferenciable)
fit_mat15 <- likfit(
  geodata      = geo_res,
  cov.model    = "matern",
  kappa        = 1.5,
  ini.cov.pars = c(2000, 100),
  nugget       = 500,
  lik.method   = "ML"
)

## Comparación formal por AIC
data.frame(
  Modelo = c("Matérn ν = 0.5 (Exp)", "Matérn ν = 1.5"),
  LogLik = c(fit_mat05$loglik, fit_mat15$loglik),
  AIC    = c(fit_mat05$AIC,    fit_mat15$AIC)
)

## Los resultados son:
##
##   Modelo              LogLik    AIC
##   Matérn ν = 0.5    −459.02   926.03   ← menor AIC
##   Matérn ν = 1.5    −459.53   927.07
##
## Como ambos modelos tienen el mismo número de parámetros (ν está fijo,
## no se estima), AIC = −2·loglik + 2·p con p idéntico. La diferencia
## en AIC se explica entonces solo por la diferencia en log-verosimilitud:
##   ΔAIC = 927.07 − 926.03 ≈ 1.04
##   ΔlogLik = 459.53 − 459.02 ≈ 0.51
##
## Podemos ver que ν = 0.5 tiene un AIC marginalmente menor, lo que indica
## una preferencia débil por el modelo Exponencial (no diferenciable).
## Sin embargo, con ΔAIC ≈ 1, la regla de Burnham & Anderson (2002)
## clasifica ambos modelos como esencialmente equivalentes (modelos con
## ΔAIC < 2 tienen soporte empírico sustancial). En la práctica, con
## solo 85 observaciones, los datos de Wolfcamp no tienen suficiente
## información para discriminar con claridad entre ν = 0.5 y ν = 1.5.
##
## La leve preferencia por ν = 0.5 sugiere que el campo de carga
## piezométrica tiende a ser no diferenciable en media cuadrática:
## las realizaciones son continuas pero con cambios de pendiente abruptos,
## consistente con la heterogeneidad geológica local del acuífero.
## Un campo con ν = 1.5 sería una vez diferenciable, produciendo
## superficies kriging más suaves, pero los datos no proveen evidencia
## suficiente para preferirlo sobre ν = 0.5.

# Parte e: Anisotropía geométrica vía MV ----------------------------------
## En la Ayudantía 02 detectamos visualmente que los variogramas
## direccionales son distintos: mayor continuidad en 135° (NW-SE) que en
## 45° (NE), dirección del gradiente. Esto es la firma de anisotropía
## geométrica: la correlación espacial decae a distinto ritmo según la
## dirección.
##
## Un modelo con anisotropía geométrica introduce una transformación
## de coordenadas antes de calcular la distancia. Los parámetros
## adicionales son:
##   psi.A: ángulo (en grados) del eje de MAYOR correlación
##   psi.R: razón de anisotropía = rango_menor / rango_mayor ∈ (0, 1]
##
## Si psi.R = 1 → isotrópico. Si psi.R < 1 → el rango en la dirección
## perpendicular (psi.A + 90°) es psi.R × phi_mayor.
##
## Ajustamos primero el modelo isotrópico como referencia y luego
## el anisotrópico, usando como valores iniciales lo que observamos
## en la Ayudantía 02: eje de mayor continuidad ~135°.

## Modelo isotrópico (referencia)
fit_iso <- likfit(
  geodata      = geo_res,
  cov.model    = "exponential",
  ini.cov.pars = c(2000, 100),
  nugget       = 500,
  lik.method   = "ML"
)

## Modelo anisotrópico — verosimilitud perfilada con transformación manual
## La función likfit de geoR tiene un bug de compatibilidad con R 4.x en
## todos los caminos de código que involucran anisotropía (el operador &&
## interno no tolera vectores de longitud 2). La solución es implementar
## la transformación anisotrópica de coordenadas directamente, sin usar
## el argumento aniso.pars de geoR.
##
## La anisotropía geométrica equivale a una transformación del espacio:
## dado un eje de mayor correlación en dirección ψ y una razón ρ (= rango
## menor / rango mayor), las coordenadas transformadas son:
##
##   x' =  x·cos(ψ) + y·sin(ψ)
##   y' = (−x·sin(ψ) + y·cos(ψ)) / ρ
##
## En el espacio (x', y') el proceso es isotrópico, y podemos ajustar
## likfit normalmente. El loglik de cada combinación (ψ, ρ) es la
## verosimilitud perfilada sobre esos parámetros.

transform_aniso <- function(coords, angle_deg, ratio) {
  a   <- angle_deg * pi / 180
  x_t <-  coords[, 1] * cos(a) + coords[, 2] * sin(a)
  y_t <- (-coords[, 1] * sin(a) + coords[, 2] * cos(a)) / ratio
  cbind(x_t, y_t)
}

angles <- c(0, 45, 90, 135)
ratios <- c(0.3, 0.5, 0.7, 0.9)

aniso_grid    <- expand.grid(angle = angles, ratio = ratios)
aniso_results <- vector("list", nrow(aniso_grid))

for (i in seq_len(nrow(aniso_grid))) {
  ## Creamos un geodata con coordenadas transformadas
  geo_t        <- geo_res
  geo_t$coords <- transform_aniso(geo_res$coords,
                                  aniso_grid$angle[i],
                                  aniso_grid$ratio[i])

  fit <- tryCatch(
    suppressMessages(
      likfit(
        geodata      = geo_t,
        cov.model    = "exponential",
        ini.cov.pars = c(2000, 100),
        nugget       = 500,
        lik.method   = "ML"
      )
    ),
    error = \(e) NULL
  )

  if (!is.null(fit)) {
    aniso_results[[i]] <- data.frame(
      angle   = aniso_grid$angle[i],
      ratio   = aniso_grid$ratio[i],
      loglik  = fit$loglik,
      aic     = fit$AIC,
      sigmasq = fit$sigmasq,
      phi     = fit$phi,
      tausq   = fit$tausq
    )
  }
}

aniso_results <- bind_rows(aniso_results)

## Tabla completa ordenada por aic
aniso_results |> arrange(aic)

## Mejor combinación anisotrópica
best_aniso <- aniso_results |> slice_min(aic, n = 1)
best_aniso

## Reajustamos likfit con las coordenadas transformadas del mejor par
## para obtener el objeto completo y poder llamar summary()
geo_best        <- geo_res
geo_best$coords <- transform_aniso(geo_res$coords,
                                   best_aniso$angle,
                                   best_aniso$ratio)

fit_aniso <- suppressMessages(
  likfit(
    geodata      = geo_best,
    cov.model    = "exponential",
    ini.cov.pars = c(2000, 100),
    nugget       = 500,
    lik.method   = "ML"
  )
)

summary(fit_aniso)

## Del summary del modelo anisotrópico (ángulo = 135°, razón = 0.7):
##
##   sigmasq (σ²):    3684 ft²
##   phi (φ):           50.42 km  → rango práctico = 3φ = 151.1 km
##   nugget (τ²):      792.7 ft²
##   AIC:              925.5
##
## Nota sobre la interpretación: geoR reporta "anisotropy angle = 0,
## ratio = 1" porque no sabe de nuestra transformación manual — ajustó
## un modelo isotrópico sobre las coordenadas ya transformadas (x', y').
## El phi = 50.42 km y el rango práctico = 151.1 km corresponden al
## eje de MAYOR correlación (135°, NW-SE) en el espacio original.
## El rango en la dirección de MENOR correlación (45°, NE) es:
##   rango_menor = rango_práctico × razón = 151.1 × 0.7 ≈ 105.7 km
##
## Comparando con el modelo isotrópico (rango práctico = 111.7 km):
## el modelo anisotrópico diferencia las dos direcciones — mayor rango
## en NW-SE (151.1 km) y menor en NE (105.7 km) — mientras que el
## isotrópico promediaba ambas en un único valor intermedio.
##
## Los resultados son consistentes con el diagnóstico visual de Ay_02:
##
##   Mejor modelo: ángulo = 135°, razón = 0.7, AIC = 925.51
##
## Los tres primeros lugares corresponden todos a ángulo = 135° (NW-SE),
## con razones 0.7, 0.5 y 0.9. Esto confirma que la dirección de mayor
## continuidad espacial es NW-SE, perpendicular al gradiente de presión
## SW → NE detectado en Ay_02.
##
## La dirección 45° (NE, a lo largo del gradiente) aparece consistentemente
## al final de la tabla para todos los valores de razón, lo que tiene
## sentido: desplazarse en la dirección del gradiente implica cruzar las
## isolíneas de presión → mayor variabilidad → menor correlación.
##
## La razón estimada ρ ≈ 0.7 indica que el rango en la dirección de menor
## continuidad (45°) es aproximadamente el 70% del rango en la dirección
## de mayor continuidad (135°). La anisotropía no es extrema.
##
## Comparando con el modelo isotrópico (AIC = 926.03):
##   ΔAIC = 926.03 − 925.51 = 0.52
##
## La mejora es marginal (ΔAIC < 2), lo que indica que ambos modelos
## tienen soporte empírico similar. Con solo 85 observaciones, los datos
## no proveen evidencia fuerte para el modelo anisotrópico, aunque la
## dirección preferida (135°) es robusta y físicamente interpretable:
## la mayor continuidad en NW-SE refleja que pozos alineados en esa
## dirección comparten niveles de presión similares por estar en la
## misma "banda" del gradiente de flujo subterráneo.

summary(fit_aniso)

## Comparación isotrópico vs mejor modelo anisotrópico.
## Ambos modelos tienen el mismo número de parámetros libres (σ², φ, τ²),
## pues el ángulo y la razón de anisotropía están fijos por perfil.
## La comparación de AIC es por lo tanto directa.
data.frame(
  Modelo = c("Isotrópico", "Anisotrópico (135°, ratio = 0.7)"),
  LogLik = c(fit_iso$loglik, fit_aniso$loglik),
  AIC    = c(fit_iso$AIC,    fit_aniso$AIC)
)

## Resultados:
##   Isotrópico:               LogLik = −459.02,  AIC = 926.03
##   Anisotrópico (135°, 0.7): LogLik = −458.75,  AIC = 925.51
##   ΔAIC = 0.52
##
## El modelo anisotrópico mejora el AIC en 0.52 puntos. Con ΔAIC < 2,
## ambos modelos tienen soporte empírico similar (Burnham & Anderson, 2002),
## por lo que no hay evidencia fuerte para preferir el anisotrópico
## desde el punto de vista formal.
##
## Sin embargo, la dirección estimada (135°, NW-SE) es robusta:
## aparece como mejor ángulo en todos los valores de razón explorados,
## y es físicamente coherente con el gradiente de flujo SW → NE del
## acuífero. Con un dataset de mayor tamaño, esta diferencia direccional
## probablemente se volvería estadísticamente más significativa.
##
## Conclusión práctica: si el objetivo es kriging, incorporar la
## anisotropía en 135° mejora marginalmente la estructura del modelo
## sin costo en parsimonia, por lo que es razonable adoptarla.

## Notemos que el modelo anisotrópico tiene 2 parámetros adicionales (psi.A
## y psi.R), por lo que el AIC penaliza más. Si aun así el AIC anisotrópico
## es menor, la mejora en ajuste compensa la complejidad añadida.
##
## Del modelo anisotrópico podemos leer:
##   phi:    rango en la dirección de mayor correlación (eje psi.A)
##   psi.A:  ángulo estimado del eje de mayor continuidad (esperamos ~135°)
##   psi.R:  razón de anisotropía; rango perpendicular = phi × psi.R
##
## La interpretación física es directa: el flujo del agua subterránea sigue
## las líneas de gradiente de presión (dirección ~45°, SW → NE). Los pozos
## alineados perpendicularmente al flujo (dirección ~135°, NW-SE) comparten
## el mismo nivel de presión → mayor correlación → mayor rango en esa
## dirección. El modelo anisotrópico captura esta estructura direccional
## que el isotrópico promedia y pierde.
