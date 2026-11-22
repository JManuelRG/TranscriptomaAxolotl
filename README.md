# Análisis bioinformático de genes asociados a la regeneración de extremidades en el ajolote (Ambystoma mexicanum): expresión diferencial y conservación evolutiva

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/TU_USUARIO/TU_REPOSITORIO/blob/main/NOMBRE_DEL_CUADERNO.ipynb)

> Proyecto de licenciatura en bioinformática, ejecutado íntegramente en Google Colab a partir de datos públicos, sin necesidad de equipo de cómputo especializado. Presentado como **póster** en [NOMBRE DEL CONGRESO].
>
> **Estado:** resultados preliminares, análisis *in silico*. Genera hipótesis; no las valida experimentalmente.

| | |
|---|---|
| **Autor(es)** | Álvarez González Nicole Joseline, González Cruz Daniela Andrea, Guzman Abalos Lucy Edith ,  Hernández Ornelas Nailea Fernanda, Rivera Garnica José Manuel|
| **Institución** | UPIIG-IPN |
| **Contacto** | jmriverag@ipn.mx|

---

## Contenido

1. [Resumen](#1-resumen)
2. [Introducción](#2-introducción)
3. [Marco teórico](#3-marco-teórico)
4. [Planteamiento, hipótesis y objetivos](#4-planteamiento-hipótesis-y-objetivos)
5. [Metodología](#5-metodología)
6. [Resultados](#6-resultados)
7. [Discusión y limitaciones](#7-discusión-y-limitaciones)
8. [Conclusiones](#8-conclusiones)
9. [Reproducibilidad y estructura del repositorio](#9-reproducibilidad-y-estructura-del-repositorio)
10. [Recursos de datos](#10-recursos-de-datos)
11. [Referencias](#11-referencias)

---

## 1. Resumen

El ajolote (*Ambystoma mexicanum*) puede regenerar extremidades completas tras una amputación mediante la formación de un **blastema**, un grupo de células progenitoras que reconstruye el miembro. Este proyecto reanaliza una matriz de conteos de RNA-seq publicada por Bryant et al. (2017) para comparar blastema (n = 4) y extremidad intacta (n = 10), e identifica genes y procesos biológicos diferencialmente expresados. Además, evalúa la conservación de las proteínas de tres genes (KAZALD1, CIRBP y PRRX1) frente a humano, ratón, pez cebra y *Xenopus*.

Con PyDESeq2 se evaluaron 37,009 genes; el blastema mostró enriquecimiento de procesos de remodelación de matriz extracelular, migración y proliferación celular, con participación de las vías BMP y Wnt, y menor expresión de genes de contracción muscular. Los controles biológicos (por ejemplo *KAZALD1*, *CIRBP*, *TBX5*, *FGF8* y *SHH*) se comportaron como se esperaba. En el análisis comparativo, el homeodominio de PRRX1 resultó idéntico en las cinco especies, el dominio RRM de CIRBP mostró alta conservación y KAZALD1 presentó divergencia moderada respecto a humano, aunque conserva su arquitectura de dominios.

**Palabras clave:** ajolote, regeneración, blastema, RNA-seq, expresión diferencial, genómica comparativa, Google Colab.

---

## 2. Introducción

Los mamíferos tienen una capacidad regenerativa muy limitada, mientras que los ajolotes pueden reemplazar extremidades completas durante toda su vida. Entender qué genes y procesos sostienen esa capacidad es relevante para la biología de la regeneración y, a largo plazo, para la medicina regenerativa.

Estudiar el ajolote a nivel genómico es difícil porque su genoma es muy grande (alrededor de 32 Gb) y estuvo incompletamente secuenciado durante años. Una alternativa es trabajar con **transcriptomas ensamblados *de novo***, que permiten medir expresión génica sin alinear contra un genoma de referencia. Bryant et al. (2017) generaron un transcriptoma de este tipo a partir de 42 muestras de 16 tejidos y pusieron a disposición pública las matrices de conteos y las anotaciones.

Este proyecto aprovecha ese recurso para practicar un flujo completo de bioinformática (datos públicos, estadística, interpretación biológica y genómica comparativa) con herramientas gratuitas y un costo computacional bajo, de modo que pueda reproducirse en un entorno educativo.

---

## 3. Marco teórico

**Regeneración y blastema.** Tras la amputación, las células cercanas a la herida se reorganizan y forman un blastema, un tejido de progenitores que reconstruye las estructuras perdidas. Bryant et al. (2017) muestrearon blastemas en estadio de yema media, derivados de amputaciones proximales y distales, tras retirar el epitelio de regeneración.

**Transcriptómica sin genoma de referencia.** El ensamblaje *de novo* (Trinity; Grabherr et al., 2011) reconstruye transcritos directamente de las lecturas. La cuantificación se hace con RSEM (Li y Dewey, 2011) y la anotación funcional con Trinotate, que asigna homologías a proteínas de UniProt. Una consecuencia importante es que los nombres de gen provienen de **homología por BLAST**, y por ello no demuestran ortología.

**Expresión diferencial.** Los modelos basados en la distribución binomial negativa, como DESeq2 (Love et al., 2014) y su implementación en Python, PyDESeq2 (Muzellec et al., 2023), estiman cambios de expresión entre condiciones y corrigen por pruebas múltiples.

**Enriquecimiento funcional.** Evalúa si una lista de genes contiene más miembros de una categoría (por ejemplo, un proceso de Gene Ontology) de lo esperado por azar, mediante la prueba hipergeométrica con corrección de tasa de falsos descubrimientos (Benjamini-Hochberg).

**Genes del panel comparativo.**

- ***KAZALD1*** (inhibidor de proteasas tipo Kazal): según Bryant et al. (2017) fue el transcrito más enriquecido en blastema, y su alteración afectó la regeneración de la extremidad.
- ***CIRBP*** (proteína de unión a RNA inducible por frío): enriquecida en blastema, con un papel citoprotector descrito en el mismo estudio. Es muy similar a su parálogo RBM3.
- ***PRRX1*** (factor de transcripción con homeodominio): gen asociado al desarrollo de la extremidad, incluido aquí como ejemplo de proteína con un dominio funcional muy conservado.

---

## 4. Planteamiento, hipótesis y objetivos

### Pregunta de investigación

¿Qué genes y procesos biológicos se expresan de forma diferencial en el blastema del ajolote respecto de la extremidad intacta, y qué tan conservadas están, frente a vertebrados no regenerativos, las proteínas de un panel de genes relacionados?

### Hipótesis de trabajo

**H1 (expresión).** El blastema presenta un perfil transcripcional distinto al de la extremidad intacta, caracterizado por una mayor expresión de genes asociados a remodelación de la matriz extracelular, migración y proliferación celular y a vías de señalización del desarrollo (por ejemplo BMP, Wnt, FGF y SHH), y por una menor expresión de genes de contracción del músculo maduro.

**H0 (expresión).** No existen diferencias sistemáticas en la expresión de esas categorías funcionales entre blastema y extremidad intacta.

**H2 (conservación).** Las proteínas de los genes del panel conservan sus dominios funcionales con una identidad de secuencia mayor dentro de los dominios que en la proteína completa, de modo que la capacidad regenerativa del ajolote se asociaría más con la regulación de la expresión que con cambios en la secuencia de estas proteínas.

> **Valoración tras el análisis:** H1 recibió apoyo en este conjunto de datos. H2 recibió apoyo parcial: se cumple para PRRX1 y CIRBP, pero no para KAZALD1, cuyos dominios también son divergentes.

### Objetivo general

Identificar y caracterizar una firma génica asociada al blastema del ajolote mediante el análisis de datos públicos de RNA-seq, y evaluar la conservación de un panel de proteínas.

### Objetivos específicos

1. Procesar una matriz de conteos pública y definir los grupos blastema y extremidad intacta.
2. Realizar un análisis de expresión diferencial y verificarlo con genes de control.
3. Efectuar enriquecimiento funcional (GO y KEGG) utilizando como universo los genes evaluados.
4. Comparar las proteínas de KAZALD1, CIRBP y PRRX1 con las de otras especies mediante alineamiento, análisis de dominios y un árbol exploratorio.
5. Comunicar los resultados en un póster y dejar un flujo reproducible.

---

## 5. Metodología

### 5.1 Datos

Matriz de conteos por gen con anotación del portal Axolotl-Omics (Broad Institute), recurso del estudio de Bryant et al. (2017): **1,388,798 "genes" de Trinity × 42 muestras de 16 tejidos**. Registro del estudio en GEO: GSE92429. Los enlaces se encuentran en la [sección 10](#10-recursos-de-datos).

### 5.2 Definición de grupos

Se utilizaron 14 de las 42 muestras:

| Grupo | Muestras | n |
|---|---|---|
| **Blastema** | `distal_*` (2) y `proximal_*` (2) | 4 |
| **Extremidad intacta** | `hand_*` (3), `forearm_*` (2), `elbow_*` (2), `upperarm_*` (3) | 10 |

Las otras 28 muestras (corazón, testículo, ovario, embriones, branquias, hueso, cartílago, músculo esquelético y vasos) no se usaron en la comparación principal, pero quedan pendientes para enriquecer el proyecto.

### 5.3 Flujo de análisis

| Fase | Procedimiento | Herramientas |
|---|---|---|
| 1. Filtrado | Conteos esperados de RSEM redondeados a enteros; se conservaron los genes con al menos 10 conteos en al menos 4 muestras (de 1,388,798 a 37,199) | pandas |
| 2. Exploración | PCA sobre log2 CPM | scikit-learn, seaborn |
| 3. Expresión diferencial | Contraste blastema vs. extremidad intacta | PyDESeq2 |
| 4. Anotación | Nombre de la mejor homología UniProt extraído del ID y convertido a símbolo de gen mediante la API de UniProt | requests |
| 5. Umbrales | Significativos: padj < 0.05 y \|log2FC\| > 1. Listas de enriquecimiento (más estrictas): padj < 0.001, \|log2FC\| > 2 y baseMean ≥ 50, con un registro por gen | pandas |
| 6. Enriquecimiento | Prueba hipergeométrica con corrección Benjamini-Hochberg, usando como universo los 13,023 genes con símbolo; bibliotecas GO Biological Process 2023 y KEGG 2021 Human (Enrichr) | scipy, statsmodels |
| 7. Genes de control | Revisión de genes esperados en blastema (*KAZALD1*, *CIRBP*, *TBX5*, *FGF8*, *SHH*, *PAX7*, *PRRX1*, *HAND2*) | pandas |
| 8. Genómica comparativa | Proteínas del ajolote (TransDecoder) frente a UniProt de humano, ratón, pez cebra y *Xenopus tropicalis*; alineamiento múltiple, identidad por dominio, conservación por posición y árbol de vecinos más cercanos (exploratorio) | Biopython, MAFFT |

### 5.4 Selección de la proteína del ajolote para cada gen

TransDecoder predice varios marcos de lectura por transcrito, y algunos no corresponden a la proteína real. Por eso la isoforma de cada gen **se eligió por su identidad con la proteína humana y no por su longitud**.

| Gen | Entrada | Proteína | Criterio |
|---|---|---|---|
| KAZALD1 | `c1076561_g4` | `m.178512` (307 aa) | Transcrito más inducido; aparece en la Tabla S4 de Bryant et al. (2017) anotado como KAZD1 |
| CIRBP | `c1090351_g5` | `m.295892` (166 aa) | Mayor identidad con CIRBP (81.8%) que con RBM3 (69.4%) |
| PRRX1 | `c1067048_g2` | `m.127462` (248 aa) | Mayor identidad con la proteína humana (86.1%) |

---

## 6. Resultados

### 6.1 Calidad y estructura de los datos

El Análisis de Componentes Principales (PCA) separa claramente el blastema (4 muestras) de la extremidad intacta en el primer componente. Las posiciones de la extremidad (mano, antebrazo, codo, brazo) no se separan entre sí, y tampoco lo hacen los blastemas proximales y distales; esto respalda agruparlos.

### 6.2 Expresión diferencial

| Resultado | Valor |
|---|---|
| Genes evaluados | 37,009 |
| Significativos (padj < 0.05, \|log2FC\| > 1) | 18,896 (10,642 arriba y 8,254 abajo en blastema) |
| Umbral estricto (padj < 0.001, \|log2FC\| > 2, baseMean ≥ 50) | 3,851 (1,064 arriba y 2,787 abajo) |

Que la mitad de los genes cambie es una señal para interpretar con prudencia (ver sección 7). PyDESeq2 determinó que el ajuste paramétrico de la curva de dispersión no convergió y usó una tendencia basada en la media.

### 6.3 Genes de control

| Gen | log2FC | padj | Comentario |
|---|---|---|---|
| *KAZALD1* (`c1076561_g4`) | 10.5 | 3×10⁻¹¹⁵ | Resultado más fuerte del análisis |
| *CIRBP* (`c1090351_g5`) | 1.79 | 7×10⁻¹¹ | Inducido |
| *CIRBP* (`c1091014_g1`) | 0.88 | 6×10⁻¹⁶ | Inducido, de magnitud moderada |
| *TBX5* | 4.5 | 8×10⁻⁴¹ | Coincide con lo descrito en el estudio original |
| *HAND2* | 5.9 | 1×10⁻³⁵ | Gen de desarrollo de la extremidad |
| *PRRX1* | 3.6 | 8×10⁻²² | |
| *FGF8* | 5.9 | 4×10⁻²⁸ | |
| *SHH* | 7.6 | 2×10⁻¹⁶ | |
| *PAX7* | 2.7 | 5×10⁻⁹ | |
| *FGFR1* | 1.8 | 3×10⁻⁸ | |
| *MMP9* | 0.04 | 0.95 | Sin cambio |
| *SOX2* | 0.5 | 0.26 | No significativo |

*MSX1*, *MSX2* y *LIN28* no aparecieron entre los resultados anotados; esto no implica que no se expresen (pueden haberse filtrado, estar anotados con otro nombre o fragmentados en el ensamblaje).

### 6.4 Enriquecimiento funcional

**Mayor expresión en blastema (570 genes):**

| Proceso | Genes | padj |
|---|---|---|
| Organización de la matriz extracelular | 36/144 | 4×10⁻¹⁵ |
| Interacción matriz extracelular-receptor (KEGG) | 18/77 | 5×10⁻⁷ |
| Desarrollo del sistema esquelético | 22/120 | 2×10⁻⁶ |
| Regulación de la migración celular | 40/349 | 3×10⁻⁶ |
| Organización de fibrillas de colágeno | 12/41 | 1×10⁻⁵ |
| Regulación positiva de la proliferación celular | 37/343 | 4×10⁻⁵ |
| Regulación de la vía BMP | 11/62 | 0.001 |
| Vía Wnt no canónica / Wnt | 7/35 y 9/59 | 0.006 y 0.008 |

**Menor expresión en blastema (1,360 genes):** contracción del músculo estriado (22/43, padj 2×10⁻⁸), ensamblaje de miofibrillas, organización del sarcómero, y también respuesta inflamatoria, señalización tipo NOD y necroptosis.

### 6.5 Análisis comparativo de proteínas

**Identidad con la proteína humana, por dominio:**

| Gen | Dominio (posición en humano) | Ajolote | Pez cebra | Ratón | *Xenopus* |
|---|---|---|---|---|---|
| PRRX1 | Homeodominio (94-153) | **100%** | 100% | 100% | 100% |
| CIRBP | RRM (6-84) | **86.1%** | 77.2% | 98.7% | 87.3% |
| KAZALD1 | IGFBP N-terminal (49-129) | 55.0% | no encontrado | 95.1% | 71.2% |
| KAZALD1 | Kazal (120-170) | 47.1% | no encontrado | 88.2% | 62.7% |
| KAZALD1 | Ig tipo C2 (172-269) | 49.0% | no encontrado | 88.8% | 76.5% |

- **CIRBP:** en un árbol exploratorio de vecinos más cercanos, los transcritos del ajolote agrupan del lado de CIRBP y no de RBM3. Con tan pocas secuencias el árbol no resuelve la filogenia entre especies y solo se usa para esta comparación.
- **KAZALD1:** la proteína del ajolote conserva el orden de dominios (péptido señal, IGFBP, Kazal, inmunoglobulina), pero su identidad con humano es de ~44%, frente a 85% del ratón y 65% de *Xenopus*. Existe otra entrada con mayor identidad (`c1049185_g1`, 64% con humano y 74% con *Xenopus*) que **no** cambia en blastema. Los cuatro candidatos se parecen más a KAZALD1 que a IGFBP7. No se sabe si hay un parálogo, una divergencia real del gen o una limitación de la predicción.

---

## 7. Discusión y limitaciones

**Interpretación.** El perfil del blastema (matriz extracelular, migración, proliferación y señales del desarrollo) es coherente con un tejido de progenitores en remodelación. La recuperación de *KAZALD1*, *CIRBP* y *TBX5*, descritos como enriquecidos en blastema por el estudio original, aporta una verificación de cordura del flujo de análisis, aunque no es evidencia independiente, porque varios de esos controles se tomaron de la misma publicación. Que el homeodominio de PRRX1 sea idéntico en las cinco especies y que el RRM de CIRBP esté muy conservado es compatible con que la diferencia regenerativa resida en la regulación de la expresión y no en la secuencia de estas proteínas, pero estos datos por sí solos no lo demuestran.

**Limitaciones.**

- **Posible efecto de lote.** Las muestras de blastema parecen provenir de otra corrida de secuenciación (los nombres incluyen identificadores de carril), y condición y lote coinciden por completo, por lo que no pueden separarse estadísticamente. La coherencia biológica de los resultados reduce, pero no elimina, esta preocupación.
- **Composición celular.** Se compara blastema (sin epitelio de regeneración) con extremidad completa, que contiene músculo, hueso y piel. Esto explica en parte los genes de menor expresión, especialmente los de contracción muscular. El estudio original plantea además que el blastema podría reprimir activamente un programa muscular; ambas explicaciones son posibles y no se distinguen con estos datos.
- **Pocas réplicas** (4 contra 10) y una fracción alta de genes significativos.
- **Anotación por homología.** Solo ~35% de los genes evaluados tenía símbolo asignado tras la conversión; los nombres provienen de mejores coincidencias BLAST, a veces de especies distintas, y se compararon contra bases humanas.
- **Genes con varias entradas** en el ensamblaje (fragmentos, isoformas o parálogos): se usó el registro más significativo por gen.
- **Umbrales** de corte elegidos por los autores del proyecto, que afectan a las listas de enriquecimiento.
- El árbol filogenético es **exploratorio** (vecinos más cercanos, pocas secuencias).
- Análisis *in silico*: sin validación experimental.

---

## 8. Conclusiones

1. Con datos públicos y herramientas gratuitas en Google Colab fue posible reproducir una firma transcripcional del blastema del ajolote que coincide con lo descrito en el estudio original.
2. El blastema mostró enriquecimiento de procesos de remodelación de matriz extracelular, migración y proliferación celular, con participación de las vías BMP y Wnt, y menor expresión de genes musculares, en parte atribuible a diferencias de composición del tejido.
3. *KAZALD1* fue la señal más fuerte del análisis y su proteína conserva la arquitectura de dominios IGFBP–Kazal–Ig, aunque con divergencia moderada de secuencia respecto a humano.
4. El homeodominio de PRRX1 y el dominio RRM de CIRBP están muy conservados entre especies, lo que sugiere que la diferencia regenerativa podría depender más de la regulación de la expresión que de la secuencia de estas proteínas (hipótesis que requiere validación experimental).
5. Los resultados deben considerarse **preliminares** por el posible efecto de lote, el número reducido de réplicas y la anotación por homología.

**Trabajo futuro:** cuantificar desde lecturas crudas con kallisto o salmon; incorporar datos de series temporales de regeneración; comparar con pez cebra; analizar datos de célula única; esclarecer la relación entre las entradas de la familia KAZALD; y validar experimentalmente los genes candidatos.

---

## 9. Reproducibilidad y estructura del repositorio

```text
.
├── README.md
├── cuaderno_ajolote.ipynb        # cuaderno de Colab con todo el flujo
├── figuras/                      # PCA, volcano, enriquecimiento, heatmap, conservación
├── resultados/                   # tablas de resultados (CSV)
└── .gitignore                    # no subir los archivos .gz de los datos
```

**Cómo reproducir:**

1. Abre el cuaderno en Colab con el botón del inicio de este documento.
2. Descarga los datos con `wget` desde los enlaces de la sección 10. Los archivos pesan unos 20-30 MB cada uno y **no se incluyen en el repositorio**.
3. Ejecuta las celdas en orden.

**Notas técnicas:**

- Si instalas `pydeseq2` y aparece un error de `pandas` (`StringDtype`), **reinicia el entorno de ejecución** después de instalar y vuelve a importar, sin repetir el `pip install`.
- MAFFT se instala con `apt-get update` y `apt-get install -y mafft`.
- Las consultas masivas a la API de UniProt pueden devolver el error 429; el código espera y reintenta.
- Paquetes principales: `pandas`, `numpy`, `scipy`, `statsmodels`, `scikit-learn`, `seaborn`, `matplotlib`, `pydeseq2`, `biopython`, `requests`.

---

## 10. Recursos de datos

### Datos del ajolote (Axolotl-Omics, Broad Institute; Bryant et al., 2017)

| Recurso | Uso | Enlace |
|---|---|---|
| Portal Axolotl-Omics | Origen de los archivos | https://portals.broadinstitute.org/axolotlomics/ |
| Conteos por gen con anotación | **Dataset principal** | https://data.broadinstitute.org/Trinity/SalamanderWeb/Axo.rsem.genes.counts.wAnnot.matrix.gz |
| Conteos por gen sin anotación | Alternativa | https://data.broadinstitute.org/Trinity/SalamanderWeb/Axo.rsem.genes.counts.matrix.gz |
| FPKM por gen con anotación | Visualización | https://data.broadinstitute.org/Trinity/SalamanderWeb/Axo.rsem.genes.TMM.fpkm.wAnnot.matrix.gz |
| Anotación funcional (Trinotate) | Dominios y homologías | https://data.broadinstitute.org/Trinity/SalamanderWeb/Axolotl.Trinity.CellReports2017.Trinotate.xls.gz |
| Proteínas predichas (TransDecoder) | Parte comparativa | https://data.broadinstitute.org/Trinity/SalamanderWeb/Axolotl.Trinity.CellReports2017.transdecoder.pep.gz |
| Transcriptoma ensamblado (FASTA) | Referencia (archivo pesado, no usado) | https://data.broadinstitute.org/Trinity/SalamanderWeb/Axolotl.Trinity.CellReports2017.fasta.gz |
| Heatmap interactivo de genes enriquecidos en blastema | Contraste | https://portals.broadinstitute.org/axolotlomics/morpheus.blastema_enriched.html |
| Registro GEO del estudio (GSE92429) | Referencia de las muestras | https://www.omicsdi.org/dataset/geo/GSE92429 |

### Genoma del ajolote (solo referencia; no se usa en este proyecto)

| Ensamblaje | Enlace |
|---|---|
| GCA_002915635.3 (AmbMex60DD), UCSC GenArk | https://hgdownload.soe.ucsc.edu/hubs/GCA/002/915/635/GCA_002915635.3/html/GCA_002915635.3_AmbMex60DD.description.html |
| Recursos genómicos del Ambystoma Genetic Stock Center | https://ambystoma.uky.edu/genome-resources |

### Proteomas de especies comparativas (UniProt)

| Especie | ID de proteoma | Enlace |
|---|---|---|
| Humano (*Homo sapiens*) | UP000005640 | https://www.uniprot.org/proteomes/UP000005640 |
| Ratón (*Mus musculus*) | UP000000589 | https://www.uniprot.org/proteomes/UP000000589 |
| Pez cebra (*Danio rerio*) | UP000000437 | https://www.uniprot.org/proteomes/UP000000437 |
| *Xenopus tropicalis* | UP000008143 | https://www.uniprot.org/proteomes/UP000008143 |

Proteínas de referencia utilizadas: KAZALD1 humano Q96I82 y ratón Q8BJ66; CIRBP humano Q14011 y ratón P60824; PRRX1 humano P54821 y ratón P63013; RBM3 humano P98179; IGFBP7 humano Q16270.

---

## 11. Referencias

1. Bryant, D. M., Johnson, K., DiTommaso, T., et al. (2017). A tissue-mapped axolotl de novo transcriptome enables identification of limb regeneration factors. *Cell Reports*, 18, 762-776. https://doi.org/10.1016/j.celrep.2016.12.063
2. Leigh, N. D., Dunlap, G. S., Johnson, K., et al. (2018). Transcriptomic landscape of the blastema niche in regenerating adult axolotl limbs at single-cell resolution. *Nature Communications*, 9, 5153. https://doi.org/10.1038/s41467-018-07604-0
3. Stewart, R., Rascón, C. A., Tian, S., et al. (2013). Comparative RNA-seq analysis in the unsequenced axolotl: the oncogene burst highlights early gene expression in the blastema. *PLoS Computational Biology*, 9(3), e1002936.
4. Smith, J. J., Timoshevskaya, N., Ye, C., et al. (2019). The axolotl genome and the evolution of key tissue formation regulators. *Genome Research*, 29, 317-324. https://genome.cshlp.org/content/29/2/317.abstract
5. Grabherr, M. G., Haas, B. J., Yassour, M., et al. (2011). Full-length transcriptome assembly from RNA-seq data without a reference genome. *Nature Biotechnology*, 29, 644-652.
6. Li, B., & Dewey, C. N. (2011). RSEM: accurate transcript quantification from RNA-seq data with or without a reference genome. *BMC Bioinformatics*, 12, 323.
7. Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550.
8. Muzellec, B., Teleńczuk, M., Cabeli, V., & Andreux, M. (2023). PyDESeq2: a Python package for bulk RNA-seq differential expression analysis. *Bioinformatics*, 39(9), btad547.
9. Chen, E. Y., Tan, C. M., Kou, Y., et al. (2013). Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool. *BMC Bioinformatics*, 14, 128.
10. Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7. *Molecular Biology and Evolution*, 30, 772-780.
11. Cock, P. J. A., Antao, T., Chang, J. T., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25, 1422-1423.
12. The UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase in 2023. *Nucleic Acids Research*, 51, D523-D531.

> **Nota:** verifica los datos bibliográficos completos (autores, volumen, páginas y DOI) de cada referencia antes de publicar el repositorio. Las referencias 1 y 3 se contrastaron con el propio artículo de Bryant et al.; el resto debe comprobarse.

---

## Licencia y cómo citar

Licencia del código: [elige una, por ejemplo MIT]. Los datos originales pertenecen a sus autores y se rigen por los términos de uso de los portales indicados.

Si utilizas este material, cita el artículo de Bryant et al. (2017) y este repositorio: *[Autor(es)] (2026). Análisis bioinformático de la regeneración de extremidades en el ajolote. GitHub: https://github.com/TU_USUARIO/TU_REPOSITORIO*
