# Explicacion de Codigos paso a paso

## Archivos de entrada

**`output_test.root`**

Lo que salio de correr Delphes (senal) con nuestros datos de Whizard a 280 GeV. Lo corrimos con una card generica de Delphes, no la de Umar.

**`GammaGammaHHESpreadAllSiD2024XCC.root`**

Es el output de Delphes de Santiago que corrio a 380 GeV. Estos archivos estan en el Drive ya que son muy grandes. Segun entiendo, Santiago si corrio con la card especifica.

Explico lo de la card porque este archivo tiene algo que nuestra salida de Delphes no tiene: los jets medidos con AntiKt. Esto hace que no pueda comparar las distribuciones de AntiKt.

---

## Corridas para extraer variables

### 280 GeV

**`ExtractVariables_Etapa1.C`** → **`variables_etapa1.root`**

Primeras pruebas en las que solo extraigo informacion de los archivos de entrada.

**`ExtractVariables_Etapa2.C`** → **`variables_etapa2.root`**

Tiene todo lo que tiene el codigo **`FSRGammaGammaHHbbbb.C`** del repositorio `XCC_HH_bbbb`, excepto:
- **Kinematic Fit**: me falta el archivo `f2c.h` que entiendo que Santiago no subio al repositorio, o tal vez este dentro del `epConstrainHH.so` en `XCC_HH_bbbb/delphes_config/modules` y yo no lo se leer.
- **Variables de AntiKt**: no disponibles porque nuestro archivo de Delphes no tiene el branch `JetAntiKt`.

Yo compare esta salida con la de 380 GeV de Santiago y me dio los graficos que estan en la carpeta `Graficos Comparacion/280`.

---

### 380 GeV

**`ExtractVariables_HHbbbb_380GeV.C`** → **`variables_HHbbbb_380GeV.root`**

Usa el output de Delphes de Santiago y extrae las variables. Es el archivo que tiene las variables que buscamos con el output de Delphes de Santiago. Esto nos deja comparar mejor y ver si estamos teniendo algun error en la extraccion de variables.

Yo compare esta salida con la de 380 GeV de Santiago y me dio los graficos que estan en la carpeta `Graficos Comparacion/380`.

---

## Comparacion

**`outputTreeSHHbbbbESpreadDurham1034BSplitSampleN.root`**

Archivo con variables definidas por Santiago en su codigo **`FSRGammaGammaHHbbbb.C`**. Es con el que vamos a comparar nuestros outputs y variables.

Para compararlos usamos estos codigos:

**`CompararVariables.C`**

Compara las variables a 280 GeV nuestras con las de 380 GeV de Santiago.

**`CompararVariables_380.C`**

Compara las variables a 380 GeV "nuestras" con las de 380 GeV de Santiago. Esta solo compara la extraccion de variables.

---

## `README_Code_Comparison.md`

Compara las partes compartidas del codigo de Santiago y el mio. Lo hice porque el codigo de Santiago es muy largo.
