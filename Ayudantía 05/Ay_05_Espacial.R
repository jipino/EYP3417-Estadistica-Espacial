## Ay_05_Espacial.R -------------------------------------------------------
## EYP3417 — Estadística Espacial
## Ayudantía 05: Pipeline geoestadístico completo — Isla Rongelap
## Profesor: Alfredo Alegría  |  Ayudante: Juan Pino

# 0. Setup ----------------------------------------------------------------
library(geoR)
library(gstat)
library(tidyverse)
library(patchwork)

## Cargamos el dataset rongelap (paquete geoR). Contiene 157 mediciones de
## radiación gamma en la isla Rongelap (Pacífico), contaminada por pruebas
## nucleares estadounidenses en los años 1950.
##
## Variables relevantes:
##   coords:  matriz 157 × 2 con coordenadas espaciales (metros)
##   data:    conteos de partículas detectadas en cada sitio
##   units.m: tiempo de exposición en segundos en cada sitio
##
## Construimos la tasa de radiación:
##   r(xᵢ) = data_i / units.m_i   (conteos por segundo)
## Esta normalización elimina el efecto del tiempo de exposición variable
## entre sitios y es la variable de interés para todos los ítems siguientes.
data(rongelap)

df_ron <- data.frame(
  x    = rongelap$coords[, 1],
  y    = rongelap$coords[, 2],
  rate = rongelap$data / rongelap$units.m
)

## Dataset: 157 sitios en un dominio de ~6.000 × 3.430 m. La tasa varía
## entre 0.25 y 15.10 ct/s (media = 7.60 ct/s), evidenciando heterogeneidad
## espacial marcada que motivará el análisis variográfico.

# Parte a: Análisis exploratorio y tendencia espacial --------------------
## Comenzamos con un gráfico de dispersión espacial coloreado por r(xᵢ).
## Si existe un gradiente espacial, la transición de colores seguirá una
## dirección preferencial — lo que motivaría ajustar una tendencia lineal
## antes de modelar la estructura de dependencia espacial (variograma).
## Sin tendencia removida, el variograma mezcla la variabilidad de gran
## escala (drift) con la dependencia local, sobreestimando el sill.

ggplot(df_ron, aes(x = x, y = y, color = rate)) +
  geom_point(size = 2.5) +
  scale_color_viridis_c(option = "plasma") +
  coord_equal() +
  theme_minimal() +
  labs(
    title    = "Dispersión espacial — Tasa de radiación gamma",
    subtitle = "Isla Rongelap · 157 sitios de monitoreo",
    x        = "X (m)", y = "Y (m)",
    color    = "Tasa\n(ct/s)"
  )

## El gráfico revela tres rasgos importantes:
##
## 1. Geometría de la isla: Rongelap tiene una forma muy irregular — una
##    franja diagonal orientada SW→NE (de x≈−6.000, y≈−3.300 hasta
##    x≈−200, y≈0) más un ramal denso al noreste (x≈−800 a −50,
##    y≈−500 a −1.800). Esta geometría NO es rectangular: el bounding
##    box de la grilla de predicción (Parte d) cubrirá muchos puntos fuera
##    del dominio real de la isla — hay que tenerlo en cuenta al interpretar
##    los mapas de kriging.
##
## 2. Patrón espacial marcado: sí existe un gradiente visual. Los valores
##    más bajos (púrpura, ~2–4 ct/s) se concentran en el extremo suroeste
##    (x < −5.500) y en algunos puntos aislados del noreste. Los valores
##    más altos (amarillo/naranja, 10–15 ct/s) aparecen en dos clusters
##    bien definidos: uno al noreste (x≈−600, y≈−1.600) y otro dentro del
##    extremo SW (x≈−5.700). La diagonal central tiene tasas intermedias
##    (salmón, ~5–8 ct/s) que aumentan suavemente de oeste a este.
##
## 3. Gradiente E-W no lineal: la tendencia lineal en x captura la tendencia
##    global (tasas crecen hacia el este en la franja central), pero no
##    puede capturar los dos clusters de alta contaminación ni la variabilidad
##    local. Esto explica el R² bajo (0.15): la mayor parte de la variabilidad
##    es de naturaleza espacial local, no un drift global lineal.

## Ajustamos el modelo de tendencia lineal:
##   r(x) = β₀ + β₁x + β₂y + ε(x)
## por Mínimos Cuadrados Ordinarios. Revisamos el summary:
##   R²: fracción de variabilidad total de r explicada por (x, y)
##   p-valores de β₁, β₂: si son < 0.05, la tendencia es estadísticamente
##     significativa → trabajar con residuos ε̂(xᵢ) en los ítems siguientes

model_trend <- lm(rate ~ x + y, data = df_ron)
summary(model_trend)

## Resultados:
##   β₀ = 6.40  (p < 0.001): intercepto
##   β₁ = −4.5×10⁻⁴ (p = 0.071): posible gradiente E-W, marginalmente
##     significativo al 10% pero no al 5%. Sugiere una leve disminución de
##     la tasa hacia el este (coordenadas más negativas) — consistente con
##     el patrón de contaminación conocido de la isla.
##   β₂ = 2.8×10⁻⁵ (p = 0.967): sin gradiente N-S significativo.
##   R² = 0.149: la tendencia lineal explica solo el 15% de la variabilidad
##     total. La mayor parte de la variabilidad es espacialmente estructurada
##     y no capturada por un plano.
##
## Decisión: aunque β₂ no es significativo y β₁ es marginalmente significativo,
## procedemos con los RESIDUOS para los ítems siguientes. Justificación:
##   (1) Garantiza estacionariedad en media del campo residual ε̂(xᵢ).
##   (2) Los residuos OLS tienen media exactamente 0, lo que nos permitirá
##       contrastar la media "conocida" μ_KS = 0 vs la estimada por GLS
##       en la Parte d (punto pedagógico central de la ayudantía).
##   (3) El variograma de los residuos y de la tasa son prácticamente
##       idénticos (R² bajo → remoción de tendencia mínima), por lo que
##       la decisión no afecta la estimación variográfica.

df_ron <- df_ron |>
  mutate(residuals = residuals(model_trend))

geo_res <- as.geodata(df_ron, coords.col = c(1, 2), data.col = 4)

## Confirmación: media residuos OLS = 0 (exacto por construcción MCO con intercepto).
## Varianza residuos = 6.38 ct²/s², levemente menor que var(rate) = 7.50 ct²/s²
## (diferencia explicada por la tendencia removida).

# Parte b: Variograma experimental y ajuste por mínimos cuadrados --------
## El variograma experimental estima la semivarianza:
##   γ̂(hₖ) = (1 / 2|Nₖ|) Σ_{(i,j)∈Nₖ} [ε̂(xᵢ) − ε̂(xⱼ)]²
## donde Nₖ = pares cuya distancia cae en el bin k (hₖ ± width/2).
##
## Parámetros clave:
##   width:  ancho de cada clase de distancia (bin)
##   cutoff: distancia máxima a incluir (regla práctica: ≤ mitad del
##            diámetro del dominio, garantizando suficientes pares por clase)

max_dist <- max(dist(geo_res$coords))   ## diámetro del dominio: 6.702 m

## Variograma omnidireccional base (cutoff = 50% del diámetro ≈ 3.350 m)
vario_emp <- variog(
  geo_res,
  max.dist = max_dist / 2,
  uvec     = seq(0, max_dist / 2, length.out = 20)
)

plot(vario_emp,
     main = "Variograma experimental — residuos de tasa",
     xlab = "Distancia (m)", ylab = expression(hat(gamma)(h)))

## El variograma exhibe una forma claramente NO MONÓTONA, con tres fases:
##
##   Fase 1 — Crecimiento (h = 0 → 882 m): γ̂ parte en ~4.5 ct²/s² y sube
##     hasta el máximo de ~9.2 ct²/s² en h ≈ 882 m. Este valor inicial
##     elevado (~4.5) ya es alto respecto al máximo (~9.2), lo que anticipa
##     un nugget importante — variabilidad a escala menor que la distancia
##     mínima entre sitios, coherente con los clusters de mediciones muy
##     próximas vistos en el scatter.
##
##   Fase 2 — Caída (h = 882 → 1.750 m): el variograma DESCIENDE hasta
##     γ̂ ≈ 4.0 ct²/s², casi el mismo nivel que el valor inicial. Esto es
##     el "efecto hoyo": pares de sitios separados ~1.750 m son más
##     similares entre sí que pares separados ~882 m. Conecta con la
##     geometría de la isla: los dos extremos de la franja diagonal (SW y NE)
##     tienen tasas parecidas (moderadas), mientras que los sitios a ~900 m
##     son los más heterogéneos (mezcla de zonas de alta y baja contaminación).
##
##   Fase 3 — Oscilación (h > 1.750 m): el variograma oscila entre 3.5 y
##     6.5 ct²/s² sin estabilizarse en un sill definido. Esto refleja la
##     geometría irregular de la isla: para distancias grandes, los pocos
##     pares disponibles corresponden a puntos en distintos ramales de la
##     isla con comportamientos muy distintos → estimación inestable.
##
## Consecuencia para el ajuste: el modelo exponencial (monotónamente
## creciente) no puede capturar esta forma. Por eso usamos un cutoff
## reducido a 1.500 m para el ajuste MCO/MCP (solo la Fase 1, donde el
## variograma se comporta de forma más estándar). La varianza total de
## los residuos es 6.38 ct²/s², inferior al pico del variograma (9.2) —
## esto es posible porque el variograma omnidireccional promedia
## direcciones con estructuras muy distintas.

## --- Sensibilidad al width (ancho de banda) ------------------------------
## width pequeño (muchas clases): curva detallada pero ruidosa porque
##   cada bin contiene pocos pares → estimación de γ̂ inestable
## width grande (pocas clases): curva suavizada pero pierde resolución
##   a distancias cortas → puede enmascarar el nugget y el rango inicial

vario_w_pocos <- variog(geo_res, max.dist = max_dist / 2,
                         uvec = seq(0, max_dist / 2, length.out = 10))
vario_w_muchos <- variog(geo_res, max.dist = max_dist / 2,
                          uvec = seq(0, max_dist / 2, length.out = 35))

tibble(
  dist  = c(vario_w_pocos$u, vario_emp$u, vario_w_muchos$u),
  gamma = c(vario_w_pocos$v, vario_emp$v, vario_w_muchos$v),
  tipo  = rep(
    c("Pocas clases (n=10)", "Base (n=20)", "Muchas clases (n=35)"),
    c(length(vario_w_pocos$u), length(vario_emp$u), length(vario_w_muchos$u))
  )
) |>
  ggplot(aes(x = dist, y = gamma, color = tipo, group = tipo)) +
  geom_point() + geom_line() +
  theme_minimal() +
  labs(
    title = "Sensibilidad al ancho de banda (width)",
    x = "Distancia (m)", y = expression(hat(gamma)(h)), color = NULL
  )

## Con pocas clases (n=10): la curva es más suave pero los bins son tan
## anchos que se promedian distancias muy dispares — el pico de ~882 m
## queda difuminado. Con muchas clases (n=35): se distingue mejor la
## estructura a cortas distancias y se aprecia el nugget (~4.5 ct²/s²),
## pero las últimas clases son ruidosas (n ≈ 62 pares en la última clase
## con cutoff = 3.350 m). La configuración base (n=20) ofrece un equilibrio.
##
## El gráfico muestra tres comportamientos bien diferenciados:
##
##   Pocas clases (n=10, azul): la curva más suave de las tres. El pico
##     queda atenuado (~8.3 ct²/s² vs ~9.2 en la base) porque el bin ancho
##     promedia el rango 700–1.100 m junto con distancias adyacentes más
##     "tranquilas". Para h > 1.500 m es la más estable — pero a costa de
##     perder resolución: la caída post-pico es gradual y no permite
##     identificar bien el mínimo ni la posible estructura oscilatoria.
##
##   Base (n=20, roja): captura el pico con claridad (~9.2, h ≈ 882 m),
##     la caída posterior y la oscilación moderada para h > 1.500 m. Es
##     el mejor equilibrio entre resolución y estabilidad, y confirma las
##     tres fases del variograma identificadas en el gráfico anterior.
##
##   Muchas clases (n=35, verde): la más volátil. El pico sube a ~10.2
##     porque los bins estrechos capturan de forma más "pura" las distancias
##     750–900 m donde los pares son especialmente heterogéneos. Pero para
##     h > 1.500 m la curva oscila entre 2.5 y 7.5 ct²/s² sin patrón
##     interpretable — son bins con muy pocos pares y la estimación de γ̂
##     es inestable. Esta configuración no es recomendable para ajuste de
##     modelos paramétricos.
##
## Conclusión: las tres curvas COINCIDEN en la región h < 700 m (la más
## informada, con hasta 1.200 pares por clase) y divergen solo a grandes
## distancias. Esto confirma que la estructura a corta distancia (nugget
## ~4.5, crecimiento inicial) es robusta al width. El width afecta
## principalmente la nitidez del pico y la estabilidad a grandes lags.

## --- Sensibilidad al cutoff ----------------------------------------------
## Cutoff pequeño: solo distancias cortas → no se aprecia el sill
## Cutoff grande: incluye pares muy lejanos con pocos pares → inestabilidad

vario_c_chico <- variog(geo_res, max.dist = max_dist / 4,
                         uvec = seq(0, max_dist / 4, length.out = 20))
vario_c_grande <- variog(geo_res, max.dist = max_dist * 0.75,
                          uvec = seq(0, max_dist * 0.75, length.out = 20))

tibble(
  dist  = c(vario_c_chico$u, vario_emp$u, vario_c_grande$u),
  gamma = c(vario_c_chico$v, vario_emp$v, vario_c_grande$v),
  tipo  = rep(
    c("Cutoff 25%", "Cutoff 50% (base)", "Cutoff 75%"),
    c(length(vario_c_chico$u), length(vario_emp$u), length(vario_c_grande$u))
  )
) |>
  ggplot(aes(x = dist, y = gamma, color = tipo, group = tipo)) +
  geom_point() + geom_line() +
  theme_minimal() +
  labs(
    title = "Sensibilidad al cutoff",
    x = "Distancia (m)", y = expression(hat(gamma)(h)), color = NULL
  )

## El gráfico muestra el efecto del cutoff sobre la forma observada del variograma:
##
##   Cutoff 25% (rojo, ≈ 1.675 m): la curva tiene el mayor número de puntos
##     por clase (bins más cortos, mismo n=20 clases en un rango menor), lo
##     que la hace la más ruidosa de las tres en su propio rango. Captura el
##     pico claramente (~10.9 ct²/s² en h ≈ 882 m) y la caída inicial, pero
##     termina alrededor de h = 1.675 m — justo cuando el variograma empieza
##     su segundo ascenso. No se puede determinar si el variograma estabiliza
##     en un sill o continúa oscilando. Insuficiente para decisiones de modelo.
##
##   Cutoff 50% (verde, base, ≈ 3.350 m): muestra el ciclo completo: subida,
##     pico, caída y primera recuperación. Coincide bien con el cutoff 25%
##     en la zona h < 1.500 m (los puntos se superponen en los valores de γ̂),
##     confirmando que esa región del variograma es robusta al cutoff. Para
##     h > 2.000 m comienza la oscilación descrita en el gráfico anterior.
##
##   Cutoff 75% (azul, ≈ 5.025 m): revela que para h > 3.500 m el variograma
##     SUBE de nuevo, alcanzando ~7.3 ct²/s² cerca de h = 4.000 m. Esto
##     refuerza la interpretación de "efecto hoyo" con estructura oscilatoria:
##     el variograma no colapsa a cero sino que tiene otro ciclo a gran escala.
##     Sin embargo, con las últimas clases teniendo muy pocos pares (la última
##     clase a ~5.000 m tiene ≈ 62 pares), esta recuperación podría ser un
##     artefacto estadístico.
##
## Conclusión: las tres curvas coinciden perfectamente en h < 600 m (la zona
## mejor estimada), divergen en el pico (~882 m) y muestran comportamientos
## distintos para h > 1.500 m según cuántos pares incluyen. El cutoff 50%
## (base) es la elección más equilibrada: captura la estructura principal
## sin arrastrar inestabilidades de grandes distancias. Para el ajuste
## paramétrico se usa cutoff = 1.500 m (zona creciente monótona), evitando
## que el efecto hoyo distorsione los parámetros estimados.

## --- Variogramas direccionales -------------------------------------------
## Si la correlación espacial varía según la dirección → anisotropía
## geométrica. Calculamos 0° (E-W), 45° (NE-SW), 90° (N-S), 135° (NW-SE).

vario_dir <- variog4(
  geo_res,
  max.dist = max_dist / 2,
  uvec     = seq(0, max_dist / 2, length.out = 15)
)

plot(vario_dir,
     main = "Variogramas direccionales (0°, 45°, 90°, 135°)",
     xlab = "Distancia (m)", ylab = expression(hat(gamma)(h)))

## El gráfico es inusual y requiere interpretación cuidadosa. Hay tres
## problemas evidentes originados en la geometría irregular de la isla:
##
##   135° (NW-SE): AUSENTE del gráfico. La dirección perpendicular al eje
##     principal de la isla (SW-NE) casi no tiene pares dentro del ángulo
##     de tolerancia (±22.5°) — la isla es demasiado estrecha en esa
##     dirección para generar bins con suficientes observaciones.
##
##   0° (E-W, negro): comportamiento DEGENERADO. Arranca con un valor
##     negativo (~−0.5, matemáticamente imposible para un variograma válido)
##     y salta inmediatamente a ~14 ct²/s², donde se mantiene plana hasta
##     el cutoff. La dirección E-W pura también tiene muy pocos pares
##     (la isla corre en diagonal SW-NE, no horizontalmente), por lo que
##     los bins tienen estimaciones erráticas con altísima varianza. El
##     valor ~14 corresponde probablemente a un par de clases con pocos
##     pares muy heterogéneos — no es interpretable como un sill real.
##
##   45° (NE-SW, rojo) y 90° (N-S, verde): curvas planas en ~6.4 y ~2.7
##     ct²/s² respectivamente, sin estructura clara de nugget-a-sill.
##     La dirección 45° (eje principal de la isla) tiene más pares y su
##     sill (~6.4) es coherente con la varianza total (6.38 ct²/s²). La
##     dirección N-S tiene un sill menor (~2.7) porque sitios alineados
##     N-S dentro de la estrecha franja de la isla tienden a ser más
##     similares entre sí (menos variabilidad N-S que NE-SW).
##
## Conclusión: los sills muy distintos entre direcciones (2.7 vs 6.4 vs 14)
## constituyen evidencia de ANISOTROPÍA GEOMÉTRICA FUERTE — la variabilidad
## espacial no es la misma en todas las direcciones. Sin embargo, la
## calidad de las estimaciones direccionales es muy pobre para este dominio
## irregular, y modelar la anisotropía requeriría técnicas más avanzadas.
## En este análisis procedemos con el modelo isótropo solicitado en el
## enunciado, siendo conscientes de esta limitación.

## --- Ajuste MCO y MCP ----------------------------------------------------
## Nota: el ajuste se realiza sobre un variograma con cutoff = 1.500 m.
## Razón: más allá de 1.500 m el variograma cae (efecto hoyo), y un modelo
## exponencial monótono no puede capturarlo. Usar la región creciente inicial
## (hasta el pico en h ≈ 882 m) da ajustes más estables y parámetros
## interpretables (nugget visible, rango apreciable).
##
## MCO (weights = "equal"): minimiza Σₖ [γ̂(hₖ) − γ_θ(hₖ)]² sin ponderación.
## MCP (weights = "npairs"): pondera por el número de pares |Nₖ|. Las clases
## con más pares (distancias cortas e intermedias) pesan más, dando un ajuste
## más robusto ante clases con pocas observaciones.

vario_fit <- variog(
  geo_res,
  max.dist = 1500,
  uvec     = seq(0, 1500, length.out = 20)
)

fit_mco <- variofit(
  vario_fit,
  cov.model    = "exponential",
  ini.cov.pars = c(4, 300),
  nugget       = 4.5,
  weights      = "equal"
)

fit_mcp <- variofit(
  vario_fit,
  cov.model    = "exponential",
  ini.cov.pars = c(4, 300),
  nugget       = 4.5,
  weights      = "npairs"
)

cat("=== Parámetros MCO vs MCP (modelo Exponencial, cutoff = 1.500 m) ===\n")
cat("MCO: σ² =", round(fit_mco$cov.pars[1], 4),
    " φ =", round(fit_mco$cov.pars[2], 1), "m",
    " τ² =", round(fit_mco$nugget, 4), "\n")
cat("MCP: σ² =", round(fit_mcp$cov.pars[1], 4),
    " φ =", round(fit_mcp$cov.pars[2], 1), "m",
    " τ² =", round(fit_mcp$nugget, 4), "\n")

## Gráfico comparativo MCO vs MCP
plot(vario_fit,
     main = "Variograma experimental + MCO + MCP (cutoff = 1.500 m)",
     xlab = "Distancia (m)", ylab = expression(hat(gamma)(h)))
lines(fit_mco, col = "steelblue", lwd = 2)
lines(fit_mcp, col = "firebrick", lwd = 2)
legend("bottomright",
       legend = c("MCO", "MCP"),
       col    = c("steelblue", "firebrick"),
       lwd    = 2, bty = "n")

## El gráfico muestra las curvas ajustadas sobre el variograma empírico
## con cutoff = 1.500 m (20 clases). Observaciones clave:
##
## 1. Diferencia en el nugget: la diferencia más visible entre MCO y MCP
##    es el punto de partida de la curva (nugget τ²), mostrado como el
##    segmento horizontal a la izquierda de h = 0. MCO (azul) tiene un
##    nugget considerablemente mayor que MCP (rojo). Esto refleja que MCO
##    con pesos iguales es arrastrado por los bins de h > 700 m (con pocos
##    pares pero valores extremos γ̂ ≈ 10–11), elevando la estimación global
##    y concentrando más variabilidad en el nugget. MCP pondera más los
##    bins cortos (h < 300 m, con hasta 668 pares), que tienen γ̂ ≈ 4–6,
##    y obtiene un nugget más bajo que refleja mejor esa región.
##
## 2. Sill similar: ambas curvas convergen aproximadamente al mismo sill
##    (~7.5 ct²/s²), coherente con la varianza total de los residuos
##    (6.38 ct²/s²). Que el sill sea similar pero los nuggets difieran
##    implica que MCO asigna más variabilidad al componente no estructurado
##    (noise) y MCP asigna más al componente espacial estructurado (σ²).
##
## 3. Ajuste a los datos: ambas curvas ajustan razonablemente bien en el
##    rango h ≈ 300–500 m (los puntos empíricos y las curvas coinciden
##    alrededor de γ̂ ≈ 6.5–7). Pero ninguna puede capturar el pico en
##    h ≈ 789–868 m (γ̂ ≈ 10.5): el modelo exponencial es monótono creciente
##    y simplemente satura en el sill antes de llegar a esos valores.
##    Para h > 1.000 m las curvas sobreestiman ligeramente el variograma
##    empírico (oscila por debajo del sill ajustado), consecuencia del
##    efecto hoyo que el modelo paramétrico no puede reproducir.

# Parte c: MV y selección de modelo --------------------------------------
## likfit() maximiza la log-verosimilitud gaussiana del campo completo:
##   ℓ(θ) = −(n/2)log(2π) − (1/2)log|Σθ| − (1/2)(z−μ1)ᵀΣθ⁻¹(z−μ1)
## A diferencia de MCO/MCP (que usan solo el variograma empírico), MV
## incorpora toda la información del vector de datos simultáneamente.
##
## Comparamos Exponencial vs Esférico con AIC = −2ℓ̂ + 2p.
## Usamos los parámetros MCP como valores iniciales para la optimización.

fit_exp <- suppressMessages(
  likfit(
    geodata      = geo_res,
    cov.model    = "exponential",
    ini.cov.pars = c(fit_mcp$cov.pars[1], fit_mcp$cov.pars[2]),
    nugget       = fit_mcp$nugget,
    lik.method   = "ML"
  )
)

fit_sph <- suppressMessages(
  likfit(
    geodata      = geo_res,
    cov.model    = "spherical",
    ini.cov.pars = c(fit_mcp$cov.pars[1], fit_mcp$cov.pars[2]),
    nugget       = fit_mcp$nugget,
    lik.method   = "ML"
  )
)

cat("=== Selección de modelo por AIC ===\n")
cat("Exponencial: AIC =", round(fit_exp$AIC, 2), "\n")
cat("Esférico:    AIC =", round(fit_sph$AIC, 2), "\n")
cat("Modelo seleccionado:",
    ifelse(fit_exp$AIC <= fit_sph$AIC, "Exponencial", "Esférico"), "\n\n")

fit_ml <- if (fit_exp$AIC <= fit_sph$AIC) fit_exp else fit_sph

cat("=== Parámetros θ̂ del modelo seleccionado (MV) ===\n")
cat("σ²      =", round(fit_ml$sigmasq, 4), "\n")
cat("φ       =", round(fit_ml$phi,     2), "m  (rango =",
    round(fit_ml$phi, 1), "m para Esférico)\n")
cat("τ²      =", round(fit_ml$tausq,   4), "\n")
cat("μ (GLS) =", round(fit_ml$beta[1], 4), "\n\n")

cat("--- Comparación MCO / MCP / MV ---\n")
cat(sprintf("MCO: σ² = %.4f  φ = %.1f m  τ² = %.4f\n",
            fit_mco$cov.pars[1], fit_mco$cov.pars[2], fit_mco$nugget))
cat(sprintf("MCP: σ² = %.4f  φ = %.1f m  τ² = %.4f\n",
            fit_mcp$cov.pars[1], fit_mcp$cov.pars[2], fit_mcp$nugget))
cat(sprintf(" MV: σ² = %.4f  φ = %.1f m  τ² = %.4f\n",
            fit_ml$sigmasq, fit_ml$phi, fit_ml$tausq))

## Resultados:
##   Exponencial: AIC = 738.86
##   Esférico:    AIC = 738.42  ← modelo seleccionado (ΔAIC = 0.44)
##
## La diferencia de AIC es mínima (< 2): ambos modelos son prácticamente
## equivalentes en términos de ajuste. En la práctica, una diferencia < 2
## no es evidencia contundente para preferir uno sobre otro — podríamos
## elegir el Esférico por parsimonia interpretativa (el rango φ es la
## distancia exacta donde la correlación cae a 0, más fácil de comunicar).
##
## Parámetros MV (Esférico seleccionado):
##   σ² = 1.83: varianza del proceso espacial estructurado
##   φ  = 249.5 m: rango del modelo Esférico (la correlación se anula
##     exactamente para h > 249.5 m). Sitios a más de 250 m son
##     prácticamente independientes bajo este modelo.
##   τ² = 4.72: nugget dominante. τ²/(σ²+τ²) = 4.72/6.55 ≈ 72% del sill
##     total es nugget. Esto refleja la alta variabilidad local entre sitios
##     cercanos — posiblemente por heterogeneidad sub-escala de la
##     contaminación o variabilidad en la medición.
##   μ (GLS) = −0.248: media GLS de los residuos. Difiere de 0 (media OLS)
##     porque GLS pondera por Σθ⁻¹. Esta diferencia es clave para la Parte d.
##
## Comparación MCO/MCP/MV: los tres métodos capturan un nugget alto (~3-5)
## y un rango moderado (178–250 m). MV da el nugget más alto porque optimiza
## la verosimilitud globalmente, ponderando correctamente la incertidumbre
## en todas las distancias sin depender de los bins del variograma empírico.

# Parte d: Mapas de predicción KS y KO -----------------------------------
## Construimos una grilla de alta resolución (100 × 100 = 10.000 puntos)
## sobre el bounding box del dominio de la isla.
##
## Para KS: media conocida = promedio aritmético de los residuos OLS = 0
## (exactamente 0 por construcción de MCO con intercepto).
## Para KO: media desconocida, restricción Σλᵢ = 1.
##
## La diferencia clave: μ_GLS (likfit) = −0.248 ≠ 0 = μ_arit.
## Esto generará diferencias visibles en la Parte e, especialmente en
## zonas alejadas de todos los datos donde la media domina la predicción.

x_seq    <- seq(floor(min(df_ron$x)),   ceiling(max(df_ron$x)),   length.out = 100)
y_seq    <- seq(floor(min(df_ron$y)),   ceiling(max(df_ron$y)),   length.out = 100)
grid_df  <- expand.grid(x = x_seq, y = y_seq)
grid_mat <- as.matrix(grid_df)

mu_arit <- mean(geo_res$data)   ## = 0 exacto (residuos OLS)

cat("Media aritmética de los residuos (μ_KS):", round(mu_arit, 8), "\n")
cat("Media GLS de likfit              (μ_GLS):", round(fit_ml$beta[1], 4), "\n\n")

## KS con media aritmética μ = 0
krige_sk <- krige.conv(
  geodata   = geo_res,
  locations = grid_mat,
  krige     = krige.control(
    type.krige = "SK",
    obj.model  = fit_ml,
    beta       = mu_arit
  )
)

## KO: media estimada implícitamente vía Σλᵢ = 1
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
    pred_sk = krige_sk$predict,
    var_sk  = krige_sk$krige.var,
    pred_ok = krige_ok$predict,
    var_ok  = krige_ok$krige.var
  )

## Mapas de predicción KS vs KO
p_sk <- ggplot(grid_df, aes(x = x, y = y, fill = pred_sk)) +
  geom_raster() +
  geom_point(data = df_ron, aes(x = x, y = y),
             size = 1, color = "white", shape = 21, stroke = 0.3,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() + theme_minimal() +
  labs(
    title    = "Kriging Simple — Predicción",
    subtitle = paste0("μ_KS = 0  (promedio aritmético de residuos OLS)"),
    x = "X (m)", y = "Y (m)", fill = "ε̂ (ct/s)"
  )

p_ok <- ggplot(grid_df, aes(x = x, y = y, fill = pred_ok)) +
  geom_raster() +
  geom_point(data = df_ron, aes(x = x, y = y),
             size = 1, color = "white", shape = 21, stroke = 0.3,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(option = "plasma") +
  coord_equal() + theme_minimal() +
  labs(
    title = "Kriging Ordinario — Predicción",
    x = "X (m)", y = "Y (m)", fill = "ε̂ (ct/s)"
  )

p_sk + p_ok

## Los dos mapas son casi idénticos visualmente, lo que anticipa que la
## diferencia KO − KS será pequeña (confirmado en la Parte e). Se observan
## tres rasgos importantes:
##
## 1. Dominio "vacío" uniforme: la mayor parte del bounding box (el
##    rectángulo) corresponde al mar, no a la isla. En esas zonas no hay
##    datos y la predicción colapsa RÁPIDAMENTE hacia la media asumida.
##    Para KS → μ_KS = 0 (salmón neutro). Para KO → μ_GLS = −0.248
##    (levemente más púrpura). El rango del modelo Esférico es φ = 249 m:
##    a más de ~250 m de cualquier sitio observado, la covarianza es 0 y
##    la predicción es puramente la media. Con un dominio de 6.000 m,
##    la gran mayoría del mapa es "fondo de media".
##
## 2. Estructura a lo largo de la isla: siguiendo la franja diagonal SW→NE
##    se aprecian los dos clusters ya observados en el scatter:
##    - Zonas amarillas (residuos positivos, ~+1.5 ct/s): el NE (x≈−200,
##      y≈−1.500) y parcialmente el SW — tasas por encima de la tendencia.
##    - Zonas púrpuras intensas (residuos negativos, hasta ~−3 ct/s): el
##      extremo SW (x≈−6.000, y≈−3.300) y un punto aislado al NE —
##      tasas por debajo de la tendencia lineal ajustada.
##    - La diagonal central tiene residuos moderados (salmón/naranja suave).
##
## 3. Localidad extrema de la predicción: la influencia de cada sitio decae
##    en un radio de ~250 m. Por eso los "manchones" de color en el mapa
##    son pequeños y bien localizados — el kriging no interpola suavemente
##    sobre distancias largas sino que vuelve rápido a la media. Esto es
##    consecuencia directa del nugget dominante (τ²/sill ≈ 72%): la mayor
##    parte de la variabilidad es no estructurada espacialmente.

## Mapas de varianza de predicción
p_vsk <- ggplot(grid_df, aes(x = x, y = y, fill = var_sk)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  coord_equal() + theme_minimal() +
  labs(
    title = "KS — Varianza de predicción",
    x = "X (m)", y = "Y (m)", fill = "σ²_KS"
  )

p_vok <- ggplot(grid_df, aes(x = x, y = y, fill = var_ok)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  coord_equal() + theme_minimal() +
  labs(
    title = "KO — Varianza de predicción",
    x = "X (m)", y = "Y (m)", fill = "σ²_KO"
  )

p_vsk + p_vok

## Los mapas de varianza tienen cuatro rasgos destacables:
##
## 1. Fondo uniformemente alto (amarillo): la mayor parte del dominio
##    rectangular (el mar) tiene varianza máxima — KS: ~6.3, KO: ~6.6
##    ct²/s². A más de 249 m de cualquier sitio observado, la covarianza
##    C(x₀, xᵢ) = 0 (rango del Esférico), y la varianza de kriging alcanza
##    su máximo C(0) = σ² + τ² = sill total. El rectángulo "de mar" está
##    completamente fuera del alcance de cualquier observación.
##
## 2. Reducción de varianza solo a lo largo de la isla: la franja diagonal
##    y los dos clusters NE y SW muestran colores más oscuros (menor
##    varianza). La varianza MÍNIMA (~5.4 ct²/s²) se alcanza en los dos
##    clusters densos de sitios (SW y NE): muchos vecinos cercanos reducen
##    k(x₀)ᵀ Σ⁻¹ k(x₀), aunque el nugget dominate impide bajar más.
##    El patrón pixelado en los clusters refleja la resolución de la grilla
##    (100×100) respecto al espaciado interno de los sitios observados.
##
## 3. Varianza alta en todas partes: el rango [5.4, 6.6] ct²/s² es
##    estrecho y elevado. El piso es impuesto por el nugget τ² ≈ 4.72:
##    incluso pegado a un sitio observado, la varianza de predicción no
##    puede bajar de τ² porque el nugget representa variabilidad irreducible
##    a cualquier distancia. Con τ²/(σ²+τ²) ≈ 72%, la mayor parte de la
##    varianza del proceso es nugget — el kriging nunca puede ser muy preciso.
##
## 4. σ²_KO > σ²_KS en zonas alejadas de datos: el mapa KO tiene escala
##    hasta 6.6 vs 6.3 en KS. La diferencia se concentra en el fondo
##    (mar), donde la restricción Σλᵢ = 1 de KO es más "costosa" —
##    obliga a repartir peso entre sitios lejanos con covarianza ≈ 0,
##    aumentando la varianza respecto a KS que simplemente usa μ conocida.
##    Cerca de los datos, ambas varianzas son prácticamente iguales (η ≈ 0).
##
## ¿Por qué la varianza depende del diseño de muestreo y no de los valores
## observados? Porque σ²_KS(x₀) = C(0) − kᵀΣ⁻¹k solo involucra la
## función de covarianza evaluada en distancias entre sitios — no Z(xᵢ).
## Dos diseños con las mismas ubicaciones pero valores distintos tendrán
## idénticos mapas de varianza.

# Parte e: Comparación de predictores ------------------------------------
## El mapa de diferencias Z*_KO − Z*_KS revela el costo de asumir una
## media incorrecta en KS. En este caso:
##   μ_KS = 0 (promedio aritmético, exactamente)
##   μ_GLS = −0.248 (estimado por MV via GLS)
##
## En zonas con datos densos: ambos predictores convergen porque los pesos
## λᵢ se concentran en los vecinos cercanos y la media tiene poco efecto.
## En zonas alejadas de todos los datos: KS predice hacia 0 (μ_KS),
## mientras KO predice hacia −0.248 (media estimada implícitamente).
## → Diferencia teórica en zonas sin datos ≈ μ_GLS − μ_KS = −0.248 ct/s.

grid_df <- grid_df |>
  mutate(diff_ko_ks = pred_ok - pred_sk)

ggplot(grid_df, aes(x = x, y = y, fill = diff_ko_ks)) +
  geom_raster() +
  geom_point(data = df_ron, aes(x = x, y = y),
             size = 1, color = "black", shape = 21, stroke = 0.3,
             inherit.aes = FALSE) +
  scale_fill_gradient2(
    low      = "steelblue",
    mid      = "white",
    high     = "firebrick",
    midpoint = 0
  ) +
  coord_equal() + theme_minimal() +
  labs(
    title    = "Diferencia de predictores: Z*_KO − Z*_KS",
    subtitle = "Rojo: KO > KS  ·  Azul: KO < KS  ·  Blanco: coinciden",
    x = "X (m)", y = "Y (m)",
    fill = "KO − KS\n(ct/s)"
  )

## El mapa es una ilustración empírica perfecta del resultado teórico:
##
## 1. Todo el mapa es AZUL — no hay rojo en ningún punto: KO < KS en todo
##    el dominio sin excepción. Esto confirma que cuando μ_KS = 0 > μ_GLS =
##    −0.248, KO siempre predice por debajo de KS. La diferencia KO − KS
##    es negativa en todas partes, con signo determinado por μ_GLS − μ_KS.
##
## 2. Fondo azul oscuro uniforme (el mar): diferencia ≈ −0.248 ct/s en toda
##    zona alejada de los datos. Es exactamente |μ_GLS − μ_arit| = 0.248
##    ct/s. En esas zonas, KS predice hacia μ_KS = 0 y KO predice hacia
##    μ_GLS = −0.248 — la diferencia es la constante μ_GLS − μ_KS. El azul
##    es uniforme porque todos esos puntos están fuera del rango del modelo
##    (>249 m de cualquier dato) y reciben idéntica influencia de la media.
##
## 3. Manchones blancos exactamente sobre la isla: a lo largo de la franja
##    diagonal y en los dos clusters densos, la diferencia se acerca a 0
##    (blanco). Cuanto más cerca de un sitio observado, mayor es Σλᵢ k(x₀)
##    y menor el peso de la media en ambos predictores → ambos convergen al
##    mismo valor observado. Los clusters NE y SW muestran los parches
##    blancos más grandes: más sitios = radio de influencia mayor = zona
##    donde la media importa poco, más extensa.
##
## 4. Gradiente suave data→mar: entre los puntos de datos y el fondo marino
##    se ve una transición gradual de blanco a azul oscuro en ~249 m (el
##    rango del modelo). Es la "zona de influencia" de cada sitio donde
##    la diferencia crece desde 0 hasta el máximo de −0.248.
##
## ¿Qué pasaría si se usara el estimador GLS (μ = −0.248) como media en KS?
## Si μ_KS = μ_GLS, KS y KO son teóricamente equivalentes (resultado de
## Ay_04): el mapa sería completamente blanco. La única diferencia entre
## los dos predictores es la elección de μ — y este mapa hace eso visible.

# Parte f: Validación cruzada leave-one-out (KO) -------------------------
## xvalid() implementa LOO-CV en un solo llamado: para cada sitio xᵢ,
## retira la observación y calcula Z*₋ᵢ(xᵢ) y σ²_{KO,−i}(xᵢ) usando
## los n−1 datos restantes con los mismos parámetros θ̂ (no re-estima).
##
## RMSE: error de predicción promedio, en unidades de la variable.
##   Cuantifica la precisión predictiva absoluta del modelo.
##
## MSSE: error cuadrático medio estandarizado por la varianza de kriging.
##   MSSE ≈ 1 → el modelo está bien calibrado: la varianza de kriging
##     captura correctamente la incertidumbre predictiva.
##   MSSE ≫ 1 → el modelo subestima la varianza real (demasiado "confiado"):
##     la varianza de kriging es pequeña pero los errores son grandes →
##     el modelo de covarianza sobreestima la dependencia espacial.
##   MSSE ≪ 1 → el modelo sobreestima la varianza (demasiado conservador):
##     la varianza de kriging es grande pero los errores son pequeños →
##     el modelo subestima la dependencia espacial.

xv <- xvalid(
  geodata = geo_res,
  model   = fit_ml,
  krige   = krige.control(type.krige = "OK", obj.model = fit_ml)
)

errors  <- xv$data - xv$predicted
var_loo <- xv$krige.var

RMSE <- sqrt(mean(errors^2))
MSSE <- mean(errors^2 / var_loo)

cat("=== Validación cruzada LOO — Kriging Ordinario ===\n")
cat("RMSE =", round(RMSE, 4), "ct/s\n")
cat("MSSE =", round(MSSE, 4), "(ideal ≈ 1)\n")

## Resultados:
##   RMSE = 2.49 ct/s: el error cuadrático medio de predicción es 2.49
##     conteos por segundo. Para contexto, la tasa media es 7.60 ct/s →
##     el error es ~33% de la media. Refleja la alta variabilidad local
##     capturada por el nugget dominante (τ² = 4.72).
##
##   MSSE = 1.019 ≈ 1: el modelo está muy bien calibrado. La varianza de
##     kriging predice casi perfectamente la magnitud real del error de
##     predicción. Un MSSE tan cercano a 1 confirma que el modelo Esférico
##     seleccionado por AIC es apropiado para describir la estructura
##     espacial de la contaminación gamma en Rongelap.

## Gráfico de diagnóstico: errores vs valores predichos
## Un buen modelo muestra errores sin patrón sistemático alrededor de 0.
ggplot(
  data.frame(pred = xv$predicted, err = errors),
  aes(x = pred, y = err)
) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick") +
  theme_minimal() +
  labs(
    title    = "LOO-CV KO — Errores de predicción",
    subtitle = paste0("RMSE = ", round(RMSE, 4),
                      " ct/s  ·  MSSE = ", round(MSSE, 4)),
    x = "Predicción Z*₋ᵢ (ct/s)",
    y = "Error: r(xᵢ) − Z*₋ᵢ(xᵢ)"
  )

## El gráfico muestra cuatro aspectos relevantes del comportamiento del modelo:
##
## 1. Sin sesgo sistemático: la nube de puntos se distribuye simétricamente
##    alrededor de la línea y = 0 (punteada roja). No hay tendencia creciente
##    ni decreciente de los errores con el valor predicho — el modelo no
##    sobreestima sistemáticamente en ninguna zona del rango de predicción.
##    El error medio es −0.013 ≈ 0, confirma la ausencia de sesgo global.
##
## 2. Dos outliers extremos negativos (error ≈ −7.8 y −8.2 ct/s, predicho
##    ≈ −1 ct/s): corresponden a sitios con residuos muy negativos que el
##    modelo no puede predecir bien desde sus vecinos. Cuando se retira ese
##    sitio, la predicción LOO sube hacia la media local (~−1 ct/s), pero
##    el valor real es mucho más bajo → error grande negativo. Son los sitios
##    de baja radiación aislados visibles en el scatter (zonas oscuras del
##    mapa de predicción). La localidad extrema del modelo (φ=249 m) implica
##    que si un sitio atípico no tiene vecinos dentro de ~250 m con valores
##    similares, el kriging no puede predecirlo bien.
##
## 3. Sin heteroscedasticidad clara: la dispersión de errores es
##    aproximadamente homogénea a lo largo de todo el rango de predichos
##    (−1.5 a +1.5 ct/s). No se observa el patrón en "embudo" que indicaría
##    que la varianza del error depende del nivel predicho. Esto es coherente
##    con la propiedad de la varianza de kriging de depender solo del diseño,
##    no de los valores observados.
##
## 4. Dispersión amplia (RMSE = 2.49 ct/s) pero MSSE = 1.019 ≈ 1: la
##    mayoría de los errores están entre ±3 ct/s. El RMSE parece grande
##    respecto a la media de la tasa (7.60 ct/s), pero el MSSE ≈ 1 confirma
##    que la varianza de kriging (σ²_KO ≈ 5.4–6.6, es decir σ_KO ≈ 2.3–2.6
##    ct/s) predice correctamente esa dispersión. El modelo sabe que es
##    impreciso (nugget alto) y lo cuantifica bien — está bien calibrado.
