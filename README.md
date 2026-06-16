# EYP3417 - Estadística Espacial

Repositorio de ayudantías del curso **EYP3417 - Estadística Espacial** de la Pontificia Universidad Católica de Chile.

| | |
|---|---|
| **Profesor** | Alfredo Alegría |
| **Ayudante** | Juan Pino |
| **Año** | 2026 |

---

## Contenidos

| # | Tema | Materiales |
|---|------|------------|
| [01](./Ayudantía%2001/) | Semivariograma Esférico | [Enunciado](./Ayudantía%2001/Ay-01.pdf) · [Solución](./Ayudantía%2001/Ay_01_Sol.pdf) · [Código R](./Ayudantía%2001/Ay_01_Espacial_R.ipynb) |
| [02](./Ayudantía%2002/) | Análisis Variográfico — Acuífero Wolfcamp | [Enunciado](./Ayudantía%2002/Ay_02.pdf) · [Código R](./Ayudantía%2002/Ay_02_Espacial.R) |
| [03](./Ayudantía%2003/) | Inferencia Paramétrica — MCO, MCP, MV, Matérn, Anisotropía | [Enunciado](./Ayudantía%2003/Ay_03.pdf) · [Código R](./Ayudantía%2003/Ay_03_Espacial.R) |
| [04](./Ayudantía%2004/) | Kriging Simple y Ordinario — Mapas predictivos, varianza, interpolador exacto, O(n³) | [Enunciado](./Ayudantía%2004/Ay_04.pdf) · [Código R](./Ayudantía%2004/Ay_04_Espacial.R) |
| [05](./Ayudantía%2005/) | Pipeline geoestadístico completo — Isla Rongelap (EDA, variograma, MV, KS/KO, LOO-CV) | [Enunciado](./Ayudantía%2005/Ay_05.pdf) · [Código R](./Ayudantía%2005/Ay_05_Espacial.R) |
| [06](./Ayudantía%2006/) | Datos Areales — Modelos Ising y CAR | [Enunciado](./Ayudantía%2006/Ay_06.pdf) · [Código R](./Ayudantía%2006/Ay_06_Espacial.R) |
| [07](./Ayudantía%2007/) | Campos Aleatorios de Markov, Modelos CAR y Procesos de Poisson Espaciales | [Enunciado](./Ayudantía%2007/Ay_07.pdf) · [Solución](./Ayudantía%2007/Ay_07_Sol.pdf) |
| [08](./Ayudantía%2008/) | Evaluación Crítica de Dependencia Espacial en Delincuencia Comunal | [Enunciado](./Ayudantía%2008/Ay_08.pdf) · [Código R](./Ayudantía%2008/Ay_08_Diagnostico_Espacial.R) |
| [09](./Ayudantía%2009/) | Procesos Puntuales — Intensidad, PPP, Funciones K/L y Modelos Log-lineales | [Enunciado](./Ayudantía%2009/Ay_09.pdf) · [Código R](./Ayudantía%2009/Ay_09_Espacial.R) · [Solución teórica](./Ayudantía%2009/Ayudantía%2009%20-%20Estadística%20Espacial.md) |
| [10](./Ayudantía%2010/) | Tests de CSR — Funciones G, F, K y L con Datos Reales | [Enunciado](./Ayudantía%2010/Ay_10.pdf) · [Código R](./Ayudantía%2010/Ay_10_Espacial.R) |
| [11](./Ayudantía%2011/) | Procesos Puntuales — PPP Inhomogéneo, Modelo Matérn-I y Análisis de Trampas Cámara | [Enunciado](./Ayudantía%2011/Ay_11.pdf) · [Código R](./Ayudantía%2011/Ay_11_Espacial.R) · [Solución teórica](./Ayudantía%2011/Ayudantía%2011%20-%20Solución%20ejercicios%201%20y%202.md) |

---

## Herramientas utilizadas

- **R** (con paquetes `ggplot2`, `dplyr`, entre otros)
- **Jupyter Notebooks** con kernel de R

## Cómo usar los notebooks

Los notebooks `.ipynb` están pensados para ejecutarse en [Google Colab](https://colab.research.google.com/) con kernel de R, o localmente con `IRkernel`.

```bash
# Instalar IRkernel en R
install.packages('IRkernel')
IRkernel::installspec()
```

## Cómo clonar y mantenerse actualizado

Consulta el archivo [INSTRUCCIONES.md](./INSTRUCCIONES.md) para una guía paso a paso sobre cómo instalar Git, clonar el repositorio y actualizarlo cada semana con `git pull`.

---

*Repositorio mantenido por Juan Pino.*
