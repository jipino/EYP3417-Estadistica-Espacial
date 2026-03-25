# Resumen Conceptual — Ayudantía 02
### EYP3417 Estadística Espacial · Ayudante: Juan Pino · Profesor: Alfredo Alegría

---

## 1. Visualización de datos espaciales

El primer paso en cualquier análisis geoestadístico es visualizar los datos en el espacio. Se grafica cada observación en su ubicación geográfica usando color y tamaño para representar el valor de la variable de interés.

**¿Qué buscar en este gráfico?**

- **Gradiente espacial**: si los valores altos se concentran en una zona y los bajos en otra, hay una tendencia de primer orden. La media no es constante en el dominio.
- **Agrupamientos**: si valores similares aparecen cerca entre sí, hay dependencia espacial (autocorrelación positiva).
- **Valores atípicos**: puntos con valores muy distintos a sus vecinos pueden ser errores de medición o fenómenos locales.

En el caso del acuífero Wolfcamp, se observó un gradiente claro: la carga piezométrica era alta en la zona occidental (X negativo) y disminuía hacia el norte y el este, sugiriendo una tendencia en dirección W/SW → NE. Esto fue confirmado posteriormente por el modelo de regresión lineal (R² = 0.89).

> **Regla práctica:** si en el mapa ves que los colores siguen un patrón geográfico (un lado más claro, otro más oscuro), ya tienes evidencia visual de que la media no es constante, y eso tendrá consecuencias en el variograma.

---

## 2. El variograma experimental

El variograma empírico estima cómo varía la semivarianza en función de la distancia entre pares de puntos:

$$\hat{\gamma}(h) = \frac{1}{2|N(h)|} \sum_{(i,j) \in N(h)} [Z(s_i) - Z(s_j)]^2$$

donde $N(h)$ es el conjunto de pares de puntos separados aproximadamente por una distancia $h$.

### Parámetros clave al calcular el variograma

**`width`** (ancho del lag): define cuánta tolerancia en distancia se usa para agrupar pares en una clase. Por ejemplo, `width = 30` agrupa todos los pares cuya distancia está entre 0 y 30 km en el primer lag, entre 30 y 60 km en el segundo, etc.

- `width` pequeño → más clases, mayor resolución a cortas distancias, pero menos pares por clase → estimación más ruidosa.
- `width` grande → menos clases, curva más suave, pero puede ocultar estructura a corta distancia.
- Regla práctica: cada clase de distancia debería contener al menos ~30 pares de puntos.

**`cutoff`** (distancia máxima): distancia hasta la cual se calculan los lags. Una regla habitual es usar como cutoff entre la mitad y dos tercios de la distancia máxima entre puntos del dominio.

### Estructura del variograma: qué observar

Un variograma de un proceso **estacionario de segundo orden** tiene esta forma característica:

```
γ̂(h)
  |         sill ___________
  |             /
  |            /
  |           /
  | nugget ●/
  |___________________________ h
             ↑
           range
```

| Parámetro | Qué representa | Cómo se lee en el gráfico |
|-----------|---------------|--------------------------|
| **Nugget** | Variabilidad a escala muy pequeña (menor que la distancia mínima entre observaciones) o error de medición | Valor de γ̂ cuando h → 0 (salto en el origen) |
| **Sill** | Varianza total del proceso estacionario. A partir de aquí, la distancia ya no agrega más diferencia entre puntos | Valor donde la curva se aplana horizontalmente |
| **Range** | Distancia a partir de la cual los puntos ya no están correlacionados | Distancia donde la curva alcanza el sill |

> **Punto clave:** si el variograma **no alcanza un sill** y sigue creciendo indefinidamente, eso es evidencia de **no estacionariedad**. El proceso tiene una tendencia que infla artificialmente la semivarianza a distancias grandes.

---

## 3. Por qué el variograma crece indefinidamente cuando hay tendencia

Supongamos que el proceso tiene la forma:

$$Z(s) = \mu(s) + \varepsilon(s)$$

donde $\mu(s) = \beta_0 + \beta_1 X + \beta_2 Y$ es una tendencia determinista y $\varepsilon(s)$ es un proceso estacionario de segundo orden.

Al calcular la diferencia entre dos puntos:

$$Z(s+h) - Z(s) = \underbrace{[\mu(s+h) - \mu(s)]}_{\beta_1 \Delta x + \beta_2 \Delta y} + [\varepsilon(s+h) - \varepsilon(s)]$$

El primer término **crece con la distancia**: a mayor separación entre puntos, mayor es la diferencia en la media. Al elevar al cuadrado y tomar esperanza:

$$2\hat{\gamma}(h) \approx (\beta_1 \Delta x + \beta_2 \Delta y)^2 + 2\gamma_\varepsilon(h)$$

El término cuadrático hace que $\hat{\gamma}(h) \to \infty$ sin cota superior. El variograma empírico no está midiendo solo la dependencia espacial de $\varepsilon(s)$, sino también la tendencia.

---

## 4. Sensibilidad al ancho del lag (`width`)

Al calcular el variograma con distintos valores de `width`, se obtiene un análisis de sensibilidad que permite:

1. **Detectar si la estructura varía mucho con el width**: si las curvas son muy distintas, la estimación es inestable y habría que revisar el número de pares por clase.
2. **Elegir un width apropiado**: se busca el menor width que produzca una curva razonablemente suave.
3. **Identificar artefactos**: un lag muy estrecho puede producir oscilaciones espurias que no reflejan estructura real.

Las tres curvas deberían ser *cualitativamente similares* (misma tendencia general). Si una difiere mucho, puede deberse a pocas observaciones en alguna clase de distancia.

---

## 5. Variogramas direccionales y anisotropía

El variograma omnidireccional promedia sobre todas las direcciones. Si la dependencia espacial no es igual en todas las direcciones, ese promedio puede ocultar información importante.

Los variogramas direccionales calculan $\hat{\gamma}(h)$ solo para pares de puntos separados en una dirección específica (con cierta tolerancia angular). Las direcciones convencionales son 0°, 45°, 90° y 135°, donde:

- **0°** corresponde al eje Norte-Sur (eje Y)
- **90°** corresponde al eje Este-Oeste (eje X)
- **45°** y **135°** corresponden a las diagonales NE y NW respectivamente

**¿Qué indica la anisotropía?**

Si los variogramas direccionales muestran **curvas con distinto ritmo de crecimiento** o que se aplanan a distintas distancias, el proceso es **anisotrópico**: la dependencia espacial depende de la dirección.

- La dirección donde la curva **crece más lentamente** (o se aplana antes) es la de **mayor continuidad espacial**: los valores cambian poco en esa dirección.
- La dirección donde la curva **crece más rápido** es la de **menor continuidad**: los valores cambian mucho al moverse en esa dirección.

En el caso Wolfcamp, la dirección 135° (NW-SE) mostraba una curva casi plana, mientras que 45° (NE) crecía muy pronunciadamente. Esto es consistente con el gradiente de presión en dirección W/SW → NE: moverse en la dirección del gradiente produce grandes cambios, mientras que moverse perpendicular a él (NW-SE) produce cambios mucho menores.

> **Nota:** la formulación matemática rigurosa de la anisotropía geométrica (que involucra una transformación lineal del espacio para convertir el proceso anisotrópico en uno isótropo) será vista en clases. Lo importante por ahora es saber *identificarla* en los variogramas direccionales.

---

## 6. Modelos teóricos de variograma

El variograma empírico es solo una estimación ruidosa. Para usar el variograma en interpolación (kriging) se necesita un modelo teórico válido. Los tres modelos más usados son:

### Modelo Exponencial

$$\gamma(h) = \tau^2 + \sigma^2 \left(1 - e^{-h/\phi}\right)$$

- Crece de forma aproximadamente lineal cerca del origen.
- Nunca alcanza el sill exactamente (convergencia asintótica).
- El range efectivo se define como la distancia donde $\gamma(h) = 0.95 \cdot \text{sill}$, que corresponde a $h \approx 3\phi$.

### Modelo Esférico

$$\gamma(h) = \tau^2 + \sigma^2 \left[\frac{3h}{2\phi} - \frac{h^3}{2\phi^3}\right] \cdot \mathbf{1}(h \leq \phi) + \sigma^2 \cdot \mathbf{1}(h > \phi)$$

- Crecimiento intermedio cerca del origen (entre lineal y parabólico).
- Alcanza el sill **exactamente** en $h = \phi$ (el range).
- Es el modelo más usado en geociencias.

### Modelo Gaussiano

$$\gamma(h) = \tau^2 + \sigma^2 \left(1 - e^{-h^2/\phi^2}\right)$$

- Crecimiento **parabólico** (muy suave) cerca del origen.
- Implica un proceso infinitamente diferenciable — muy suave espacialmente.
- Puede ser numéricamente inestable en kriging para datos irregulares.
- El range efectivo es $h \approx \sqrt{3}\,\phi$.

### ¿Cómo elegir el mejor modelo?

Se ajustan los modelos por mínimos cuadrados (minimizando la suma de cuadrados de los residuos entre el variograma empírico y el teórico). El modelo con menor SSR es el que mejor reproduce la forma del variograma empírico.

Cualitativamente, el comportamiento **cerca del origen** es el criterio más importante:
- ¿La curva sube suavemente (parabólico) → Gaussiano
- ¿Sube de forma lineal → Exponencial
- ¿Sube de forma intermedia → Esférico

---

## 7. Remoción de tendencia y variograma de residuos

Cuando la variable de interés tiene una tendencia espacial, la estrategia estándar es:

1. **Ajustar un modelo de tendencia**: por ejemplo, una superficie lineal $\hat{\mu}(s) = \hat{\beta}_0 + \hat{\beta}_1 X + \hat{\beta}_2 Y$ mediante regresión lineal múltiple.
2. **Calcular los residuos**: $\hat{\varepsilon}(s) = Z(s) - \hat{\mu}(s)$.
3. **Estimar el variograma sobre los residuos**: $\hat{\gamma}_\varepsilon(h)$.

Los residuos tienen media aproximadamente cero en todo el dominio y su variograma debería mostrar un sill definido, lo que indica que la estacionariedad de segundo orden ha sido recuperada.

**¿Cómo se ve este efecto en el gráfico?**

Al comparar el variograma original (sobre $Z(s)$) con el variograma de residuos (sobre $\hat{\varepsilon}(s)$), se espera ver:
- El variograma original sigue creciendo sin límite.
- El variograma de residuos se aplana en un sill mucho más bajo.
- La diferencia entre ambos cuantifica cuánto inflaba la tendencia al variograma original.

En el caso Wolfcamp, el R² = 0.89 indicó que la tendencia lineal explicaba el 89% de la varianza total. El variograma original alcanzaba ~25.000 ft² a 160 km, mientras que el de residuos se estabilizaba en ~4.000–4.500 ft² después de 100 km — una reducción de escala de aproximadamente 6 veces.

---

## 8. Tabla resumen: señales de diagnóstico

| Lo que observas | Lo que indica |
|----------------|---------------|
| Variograma crece sin sill | No estacionariedad (tendencia en los datos) |
| Variograma con sill claro | Proceso estacionario de segundo orden |
| Nugget alto | Variabilidad a escala sub-muestral o error de medición |
| Variogramas direccionales distintos | Anisotropía (dependencia espacial depende de la dirección) |
| Variogramas direccionales similares | Isotropía (dependencia espacial igual en todas las direcciones) |
| Variograma de residuos con sill | La remoción de tendencia recuperó la estacionariedad |
| width pequeño → curva oscilatoria | Pocos pares por clase de distancia |
| width grande → curva suave | Promediado excesivo, pérdida de resolución |
