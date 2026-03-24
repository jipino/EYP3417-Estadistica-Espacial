## Ay_02_Espacial.R -------------------------------------------------------
## EYP3417 — Estadística Espacial
## Ayudantía 02: Análisis de la Carga Piezométrica — Acuífero Wolfcamp
## Profesor: Alfredo Alegría  |  Ayudante: Juan Pino

# 0. Setup ----------------------------------------------------------------
library(geoR)
library(gstat)
library(tidyverse)
library(patchwork)

## Cargamos el dataset wolfcamp del paquete geoR. Este corresponde a mediciones
## de carga piezométrica (presión del agua subterránea, en pies) en 85 pozos
## del acuífero Wolfcamp, Texas. El objeto wolfcamp es de clase geodata (lista),
## por lo que extraemos manualmente las coordenadas y los datos.
data(wolfcamp)

df_wolf <- data.frame(
  x        = wolfcamp$coords[, 1],
  y        = wolfcamp$coords[, 2],
  pressure = wolfcamp$data
)

# Parte a: Visualización espacial -----------------------------------------
## Antes de calcular cualquier variograma, es importante explorar los datos
## visualmente. El objetivo es detectar si la carga piezométrica cambia
## sistemáticamente en el espacio, lo que se conoce como tendencia de primer orden.
##
## Recordemos que un proceso Z(s) es estacionario de primer orden si su media
## es constante en todo el dominio: E[Z(s)] = μ para todo s. Si el mapa muestra
## que los valores altos se concentran en una zona y los bajos en otra, esa
## condición no se cumple, y cualquier análisis variográfico sobre los datos
## brutos estará contaminado por esa tendencia.

plot_spatial_data <- function(data) {
  ggplot(data, aes(x = x, y = y, size = pressure, color = pressure)) +
    geom_point(alpha = 0.7) +
    scale_color_viridis_c(option = "magma") +
    scale_size_continuous(range = c(1, 6)) +
    theme_minimal() +
    labs(
      title    = "Carga Piezométrica — Acuífero Wolfcamp",
      subtitle = "¿Existe una dirección preferencial de cambio en la presión?",
      x        = "Coordenada X (km)",
      y        = "Coordenada Y (km)",
      color    = "Presión (ft)",
      size     = "Presión (ft)"
    )
}

p_spatial <- plot_spatial_data(df_wolf)
p_spatial

## Del gráfico podemos ver que la carga piezométrica es alta en el SW (esquina inf. izq.) (~900 ft,
## puntos amarillos) y disminuye progresivamente hacia el NE (esquina sup. der.) (~400–500 ft, puntos
## morados). Hay un gradiente espacial claro en dirección SW → NE.
##
## Esto indica que E[Z(s)] = μ(s) ≠ cte, es decir, el supuesto de media constante
## no se cumple: existe una tendencia determinista de primer orden. Si calculamos
## γ̂(h) directamente sobre los datos, estaríamos incluyendo en la semivarianza
## tanto la variabilidad del proceso aleatorio ε(s) como el efecto del gradiente
## μ(s), lo que resultaría en un variograma artificialmente inflado que nunca
## alcanza un sill. Las Partes b y c mostrarán este efecto, y la Parte e mostrará
## cómo remover la tendencia recupera la estacionariedad en los residuos.

# Parte b: Sensibilidad del variograma al ancho de lag --------------------
## El variograma γ(h) mide la semivarianza promedio entre pares de puntos
## separados por una distancia h. Para estimarlo, agrupamos los pares en clases
## de distancia de ancho `width`, de acuerdo al estimador de Matheron:
##
##   γ̂(h_k) = 1/(2·|N(h_k)|) · Σ [Z(sᵢ) − Z(sⱼ)]²
##
## donde N(h_k) es el conjunto de pares con distancia en la clase k. A mayor
## |N(h_k)| (más pares por clase), más estable es la estimación.
##
## Existe un trade-off en la elección del width: un width pequeño genera clases
## estrechas con pocos pares, lo que produce estimaciones ruidosas. Un width
## grande genera clases anchas con más pares, la curva es más suave, pero puede
## ocultar estructura a corta distancia (sesgo). La regla práctica de Cressie
## (1993) indica que cada clase debería tener al menos ~30 pares para que la
## estimación sea confiable.
##
## El argumento cutoff indica la distancia máxima considerada. Se recomienda
## no superar la mitad del diámetro del dominio, para que los lags lejanos
## aún tengan pares suficientes.

compute_variograms <- function(data, widths, cutoff = 150) {
  ## Nombramos la lista con el valor de width para identificar cada variograma
  set_names(widths, widths) |>
    map(\(w) variogram(pressure ~ 1, locations = ~ x + y,
                       data = data, width = w, cutoff = cutoff))
}

plot_variogram_sensitivity <- function(v_list) {
  ## Apilamos los variogramas en un data frame largo para comparar en un gráfico
  bind_rows(v_list, .id = "lag_width") |>
    mutate(lag_width = paste0("width = ", lag_width, " km")) |>
    ggplot(aes(x = dist, y = gamma, color = lag_width)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    theme_minimal() +
    labs(
      title    = "Sensibilidad del Variograma al Ancho de Lag",
      subtitle = "Width pequeño → más ruido; Width grande → suaviza estructura",
      x        = "Distancia (km)",
      y        = expression(hat(gamma)(h)),
      color    = NULL
    )
}

v_list        <- compute_variograms(df_wolf, widths = c(10, 30, 60))
p_sensitivity <- plot_variogram_sensitivity(v_list)
p_sensitivity

## Podemos ver que ninguno de los tres variogramas alcanza un sill: γ̂(h) sigue
## creciendo hasta el cutoff (150 km) sin estabilizarse. Esto confirma la no
## estacionariedad detectada en la Parte a, pues un proceso estacionario debería
## producir un variograma acotado.
##
## Respecto al efecto del width, con width = 10 km se obtiene mayor resolución
## a distancias cortas, pero a distancias grandes hay pocas estaciones tan
## separadas entre sí, lo que produce estimaciones muy variables (ver las
## oscilaciones entre 90–120 km). Con width = 30 km se logra un balance
## adecuado: la curva es suave y tiene suficiente resolución, por lo que
## usaremos este en las Partes d y e. Con width = 60 km la curva es muy suave,
## casi lineal, y puede ocultar curvatura real a cortas distancias.
##
## Los tres variogramas coinciden bien hasta ~50 km. Las diferencias se amplifican
## a distancias grandes, donde el ruido domina cuando el width es pequeño.

# Parte c: Variogramas direccionales --------------------------------------
## Hasta ahora calculamos variogramas omnidireccionales: promediamos todos los
## pares sin importar su orientación. Sin embargo, si la estructura espacial
## depende de la dirección en que nos desplazamos, ese promedio oculta información
## importante. Un variograma direccional restringe los pares a aquellos que
## apuntan aproximadamente hacia un ángulo α (con una tolerancia angular de
## 22.5° por defecto en gstat).
##
## Si los variogramas direccionales son similares entre sí, el proceso es
## isótropo (la dependencia espacial es igual en todas las direcciones). Si
## difieren en su ritmo de crecimiento o en el valor del sill, el proceso
## es anisótropo. La convención de ángulos en gstat es: 0° → Norte (N-S),
## 90° → Este (E-W), 45° → Noreste y 135° → Sureste (equivalente a NW-SE).
##
## Nota: la formulación matemática formal de la anisotropía geométrica,
## incluyendo la razón de anisotropía y la transformación de coordenadas que
## la corrige, será vista en la clase del miércoles. Hoy nos concentramos
## en la detección visual.

compute_directional_variograms <- function(data,
                                           angles = c(0, 45, 90, 135),
                                           cutoff = 150) {
  variogram(pressure ~ 1, locations = ~ x + y, data = data,
            alpha = angles, cutoff = cutoff)
}

plot_directional <- function(v_dir) {
  ggplot(v_dir, aes(x = dist, y = gamma, color = factor(dir.hor))) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    facet_wrap(
      ~ dir.hor,
      labeller = labeller(dir.hor = \(x) paste0(x, "°"))
    ) +
    theme_minimal() +
    labs(
      title    = "Variogramas Direccionales",
      subtitle = "Si el range difiere entre direcciones → anisotropía geométrica",
      x        = "Distancia (km)",
      y        = expression(hat(gamma)(h)),
      color    = "Dirección (°)"
    )
}

v_dir <- compute_directional_variograms(df_wolf)
p_dir <- plot_directional(v_dir)
p_dir

## Se aprecia que los cuatro variogramas direccionales se comportan de forma
## muy distinta, lo que indica que la dependencia espacial no es igual en todas
## las direcciones: el proceso es anisótropo.
##
## En la dirección 45° (NE), el variograma crece más pronunciadamente. Esto
## tiene sentido, pues es la dirección del gradiente principal identificado en
## la Parte a: al desplazarnos hacia el NE cruzamos el gradiente completo, los
## pares tienen diferencias grandes y γ̂(h) no muestra sill. En la dirección
## 135° (NW-SE, perpendicular al gradiente), en cambio, se tiene el variograma
## más plano y es el único que sugiere un sill (~6.000–7.000 ft²). Esto es
## intuitivo: al movernos en sentido NW-SE vamos "paralelos" a las isolíneas
## de presión, los valores cambian poco y la continuidad espacial es mayor.
## Las direcciones 0° y 90° presentan un comportamiento intermedio.
##
## Una forma de entender esto es pensar en las isolíneas de presión como curvas
## de nivel en un mapa topográfico: cruzar las curvas (dirección 45°) implica
## mucho cambio de valor, mientras que caminar a lo largo de una curva (135°)
## implica poco cambio.

# Parte d: Ajuste de modelos teóricos -------------------------------------
## El variograma empírico γ̂(h) es ruidoso; necesitamos un modelo teórico suave
## y matemáticamente válido (condicionalmente definido negativo) para usar en
## kriging u otras predicciones. Los tres parámetros principales de cualquier
## modelo de variograma son: el nugget (c₀), que representa la varianza a
## distancia cero y captura variabilidad a micro-escala o error de medición;
## el sill (c₀ + c), que es el valor máximo que alcanza el variograma y
## corresponde a la varianza total del proceso estacionario; y el range (a),
## que es la distancia a la que se alcanza el sill (o ~95% de él), más allá
## de la cual los puntos son prácticamente independientes.
##
## Los modelos difieren principalmente en el comportamiento cerca del origen:
##   Exponencial: γ(h) = c₀ + c·[1 − exp(−h/a)], crece linealmente cerca del
##     origen, lo que implica un proceso irregular (no diferenciable).
##   Esférico: γ(h) = c₀ + c·[3h/(2a) − h³/(2a³)] para h ≤ a (y c₀+c si h > a),
##     comportamiento intermedio, alcanza el sill exactamente en h = a.
##   Gaussiano: γ(h) = c₀ + c·[1 − exp(−h²/a²)], crece parabólicamente, lo
##     que implica un proceso muy suave (infinitamente diferenciable).
##
## Como criterio de selección usaremos el SSR (suma de cuadrados ponderada de
## los residuos entre el variograma empírico y el teórico): menor SSR = mejor ajuste.
##
## Importante: intentar ajustar estos modelos acotados sobre el variograma de
## los datos originales (pressure ~ 1) falla por convergencia, pues la curva
## no tiene sill y el algoritmo no puede encontrar parámetros que la acomoden.
## Esto confirma la no estacionariedad de las Partes a–c. La solución es
## calcular el variograma usando `pressure ~ x + y` en lugar de `pressure ~ 1`,
## lo que le indica a gstat que remueva la tendencia lineal antes de estimar
## el variograma. Los residuos sí son aproximadamente estacionarios y su
## variograma debería mostrar un sill definido.

try_fit_variogram <- function(v_emp, model,
                              psill_vals  = c(500, 1000, 2000, 3000),
                              range_vals  = c(30,  60,  100,  150),
                              nugget_vals = c(0,   200,  500)) {
  ## Probamos una grilla de 4×4×3 = 48 combinaciones de valores iniciales
  ## y nos quedamos con la que entrega el menor SSR, para evitar problemas
  ## de convergencia con valores iniciales mal especificados.
  grid <- expand.grid(psill = psill_vals, range = range_vals, nugget = nugget_vals)

  best_ssr <- Inf
  best_fit <- NULL

  for (i in seq_len(nrow(grid))) {
    fit <- tryCatch(
      withCallingHandlers(
        fit.variogram(v_emp,
                      model = vgm(psill  = grid$psill[i],
                                  model  = model,
                                  range  = grid$range[i],
                                  nugget = grid$nugget[i])),
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

fit_theoretical_models <- function(v_emp, models = c("Exp", "Sph", "Gau")) {
  set_names(models) |>
    map(\(m) try_fit_variogram(v_emp, m))
}

## Comparación visual: puntos empíricos + curvas ajustadas superpuestas
plot_fitted_models <- function(v_emp, fitted_list) {
  h_seq <- seq(0, max(v_emp$dist) * 1.05, length.out = 300)

  curves_df <- map(fitted_list, \(fit) {
    tibble(dist = h_seq,
           gamma = variogramLine(fit, dist_vector = h_seq)$gamma)
  }) |>
    bind_rows(.id = "model")

  ggplot() +
    geom_point(data = v_emp, aes(x = dist, y = gamma),
               size = 2.5, color = "grey30") +
    geom_line(data = curves_df, aes(x = dist, y = gamma, color = model),
              linewidth = 1.2) +
    theme_minimal() +
    labs(
      title    = "Ajuste de Modelos Teóricos al Variograma Empírico",
      subtitle = "Gau: parabólico (suave) | Exp: lineal | Sph: intermedio",
      x        = "Distancia (km)",
      y        = expression(hat(gamma)(h)),
      color    = "Modelo"
    )
}

## Tabla con SSR de cada modelo (menor SSR = mejor ajuste)
model_ssr <- function(fitted_list) {
  map_dbl(fitted_list, \(fit) attr(fit, "SSErr")) |>
    enframe(name = "modelo", value = "SSR") |>
    arrange(SSR)
}

v_base        <- variogram(pressure ~ x + y, locations = ~ x + y, data = df_wolf, width = 30)
models_fitted <- fit_theoretical_models(v_base)

plot_fitted_models(v_base, models_fitted)
model_ssr(models_fitted)

## Podemos ver que el Esférico obtiene el menor SSR (~19.097), seguido del
## Exponencial (~29.924), mientras que el Gaussiano es claramente inadecuado
## (SSR ~385.643). El Gaussiano alcanza el sill en apenas ~30 km, dejando
## todos los puntos de distancias grandes muy por encima de la curva, y su
## forma parabólica cerca del origen supone un proceso infinitamente diferenciable,
## que no se condice con la variabilidad observada en los pozos.
##
## El Esférico captura bien la curvatura intermedia cerca del origen y alcanza
## el sill de forma exacta, por lo que es el modelo recomendado para este dataset.
## El Exponencial es una segunda opción razonable, aunque su convergencia asintótica
## al sill lo hace ligeramente peor.
##
## En cuanto a los parámetros estimados, el nugget es de aproximadamente
## 1.200–1.300 ft², lo que refleja variabilidad a escalas menores que la
## distancia mínima entre pozos o error de medición. El sill se estima en
## ~4.000–4.500 ft², que corresponde a la varianza del proceso estacionario ε(s),
## y el range en ~100–150 km, distancia más allá de la cual dos pozos tienen
## carga piezométrica prácticamente independiente (dado el trend).

# Parte e: Análisis de tendencia y residuos -------------------------------
## El modelo que subyace a este análisis es Z(s) = μ(s) + ε(s), donde μ(s) es
## la media determinista (tendencia) y ε(s) es el proceso estacionario de media
## cero. Si μ(s) = β₀ + β₁X + β₂Y, ajustamos un plano por mínimos cuadrados y
## los residuos ε̂(s) = Z(s) − μ̂(s) deberían ser estacionarios.
##
## Para entender por qué el variograma original crece sin sill, notemos que:
##
##   Z(s+h) − Z(s) = [μ(s+h) − μ(s)] + [ε(s+h) − ε(s)]
##                 = [β₁Δx + β₂Δy]   + [ε(s+h) − ε(s)]
##
## Al tomar la semivarianza, el término de tendencia [β₁Δx + β₂Δy]² crece
## cuadráticamente con ||h||, llevando γ_Z(h) → ∞. Al calcular los residuos,
## ese término desaparece y γ_ε(h) puede alcanzar un sill finito.

analyze_residuals <- function(data) {
  model_trend <- lm(pressure ~ x + y, data = data)

  data_res <- data |>
    mutate(residuals = residuals(model_trend))

  v_res <- variogram(residuals ~ 1, locations = ~ x + y,
                     data = data_res, cutoff = 150, width = 30)

  list(data = data_res, variogram = v_res, model = model_trend)
}

compare_original_residual <- function(v_original, v_residual) {
  bind_rows(
    mutate(v_original, tipo = "Original (con tendencia)"),
    mutate(v_residual, tipo = "Residuos (tendencia removida)")
  ) |>
    ggplot(aes(x = dist, y = gamma, color = tipo)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.8) +
    theme_minimal() +
    labs(
      title    = "Variograma: Original vs Residuos",
      subtitle = "Remover la tendencia lineal debería recuperar estacionariedad",
      x        = "Distancia (km)",
      y        = expression(hat(gamma)(h)),
      color    = NULL
    )
}

res_analysis <- analyze_residuals(df_wolf)
summary(res_analysis$model)

## El modelo ajustado es Ẑ(s) = 607.77 − 1.278·X − 1.139·Y. Ambos coeficientes
## son negativos y altamente significativos (p < 2e-16): la presión cae ~1.28 ft
## por km hacia el Este y ~1.14 ft por km hacia el Norte, lo que es consistente
## con el gradiente SW → NE observado en la Parte a. El R² = 0.891 indica que
## el plano lineal explica el 89% de la varianza total, lo que confirma que la
## tendencia espacial es la fuente dominante de variabilidad en estos datos.

## Comparamos el variograma original (pressure ~ 1) con el de los residuos
v_original <- variogram(pressure ~ 1, locations = ~ x + y,
                        data = df_wolf, width = 30)
compare_original_residual(v_original, res_analysis$variogram)

## Del gráfico se aprecia claramente el efecto de remover la tendencia. El
## variograma original crece sin sill hasta ~25.000 ft² a 160 km, mientras
## que el de los residuos se estabiliza en ~4.000–4.500 ft² a partir de ~100 km,
## mostrando un sill definido. La escala cae en un factor ~6, lo que cuantifica
## cuánta varianza era atribuible a la tendencia y no al proceso estocástico.
## La presencia del sill en los residuos confirma que ε̂(s) es aproximadamente
## estacionario de segundo orden, validando el supuesto del modelo.


