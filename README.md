# flexray-physical-decoder

Decodifica tramas FlexRay a partir de la señal eléctrica capturada con un
osciloscopio, y genera señales FlexRay sintéticas para probarlo. Sin librerías
de FlexRay: todo el protocolo está implementado desde cero, del umbral de
tensión a los dos polinomios de CRC.

![tests](https://github.com/andriusvillena-cmd/flexray-physical-decoder/actions/workflows/tests.yml/badge.svg)

---

## Qué hace

Entrada: un CSV exportado de PicoScope con tres columnas — tiempo y los dos
hilos del par diferencial en voltios.

```
python flexray_decode.py FlexRay_Trace.csv
```

```
  paso de muestreo    0.0125 us   (80 MS/s)
  tiempo de bit       100.0 ns
  velocidad del bus   10.00 Mbit/s
  muestras por bit    8.0

  TSS                 11 bits
  bytes leidos        40
  fin de trama        FES correcto

  identificador       1
  tipo                sincronismo, arranque
  contador de ciclo   45
  carga util          16 palabras = 32 bytes

  CRC de cabecera     0x0F2 leido / 0x0F2 calculado   CORRECTO
  CRC de trama        0xA011D3 leido / 0xA011D3 calculado   CORRECTO
  canal               A   (deducido de la semilla del CRC que cuadra)
```

Y en sentido contrario:

```
python flexray_gen.py --id 42 --ciclo 7 --bytes 16 --patron cuenta
```

---

## Por qué FlexRay no se decodifica como CAN

**No hay relleno de bits.** CAN prohíbe más de 5 bits iguales seguidos e
inserta uno contrario para forzar un flanco. FlexRay hace otra cosa: delante de
**cada byte** manda dos bits de servicio, uno alto y uno bajo, llamados BSS. El
flanco de bajada de esos dos bits es un punto de sincronización cada 10 bits.

Eso cambia el diseño del decodificador. En CAN se mide el tiempo de bit una vez
y se muestrea toda la trama contando desde el bit de arranque. Aquí, si haces
eso, la deriva te come: a 10 Mbit/s con 8 muestras por bit, un error del 1 % en
el tiempo de bit desplaza más de medio bit antes de llegar al final. Así que
este decodificador **se recoloca en cada BSS**, exactamente como el receptor
real, y por eso funciona incluso con 4 muestras por bit.

**Hay dos CRC, no uno.** Uno de 11 bits que protege solo la cabecera, y otro de
24 que protege cabecera y carga útil. La cabecera va protegida aparte porque un
nodo tiene que poder fiarse del identificador y de la longitud *antes* de haber
recibido la trama entera.

**Y el CRC de trama no arranca de cero.** Empieza desde una constante distinta
según el canal físico: una para el A y otra para el B. Esto tiene una
consecuencia práctica bonita: si pruebas las dos y solo una cuadra, has
averiguado de qué canal es la captura sin que nadie te lo diga.

---

## Cómo funciona, paso a paso

**1. Diferencial.** Los dos hilos se restan. En reposo están al mismo
potencial y la diferencia es casi cero; eso permite acotar dónde empieza y
acaba la trama sin buscar patrones.

**2. Umbral con histéresis.** Dos rayas a ±0,3 V, y entre ellas se conserva el
estado anterior.

**3. Tiempo de bit.** Se mide de la propia señal, y **no** a partir del hueco
más corto entre flancos. Con 8 muestras por bit, cada flanco se localiza con un
error de hasta una muestra, que es el 12 % de un bit, y el mínimo sale corto por
puro redondeo. Se usa el hueco más *frecuente*, donde el error se reparte a los
dos lados, y luego se afina exigiendo que todos los huecos duren un número
entero de bits.

**4. TSS.** La trama abre con una tirada de nivel bajo, de 3 a 15 bits.

**5. Recorrido byte a byte.** FSS, y después, para cada byte: se comprueba el
BSS, se busca su flanco de bajada real, se ancla ahí y se leen los 8 bits por
su centro. El anclaje es lo que evita la deriva.

**6. FES.** La trama cierra con un bit bajo seguido de uno alto donde tocaría
un BSS. El decodificador lo distingue de un error de sincronización.

**7. Cabecera.** 40 bits: reservado, indicador de preámbulo, indicador de trama
nula, sincronismo, arranque, identificador de 11 bits, longitud de carga útil
de 7 bits en palabras de 2 bytes, CRC de cabecera de 11 bits, y contador de
ciclo de 6 bits.

**8. Los dos CRC.** Se recalculan y se comparan con los transmitidos. El de
cabecera con el polinomio 0x385 desde 0x01A, sobre 20 bits. El de trama con el
polinomio 0x5D6DCB, sobre la cabecera entera más la carga útil, probando las
dos semillas de canal.

---

## Los datos

**El repositorio no incluye ninguna captura de osciloscopio.**

El decodificador se desarrolló y se validó contra una captura real de la
biblioteca de formas de onda de PicoScope: una trama de arranque, canal A,
identificador 1, ciclo 45, 32 bytes de carga útil a cero, tomada a 80 MS/s. Las
condiciones de reutilización de esa biblioteca no están publicadas, así que el
archivo no se redistribuye.

Lo que sí se publica es `flexray_gen.py`, que fabrica señales equivalentes.
Permite probar el decodificador en cualquier máquina y en condiciones que una
sola captura no da: distintos identificadores, longitudes de carga útil de 2 a
254 bytes, los dos canales, velocidades de muestreo más bajas, ruido y
amplitudes reducidas.

Escribir el generador no fue un rodeo. Construir una trama obliga a saber qué va
en cada bit, y deja el protocolo implementado por los dos lados.

---

## El círculo cerrado, y cómo se rompe

Probar un decodificador con datos de tu propio generador tiene un riesgo obvio:
si los dos comparten el mismo malentendido, las pruebas salen verdes y el
resultado es falso.

Por eso hay dos pruebas de **anclaje**, que comparan contra valores medidos en
la captura real:

- Los 5 bytes de cabecera que produce el generador son `38 01 20 3C AD`.
- Los dos CRC valen `0x0F2` y `0xA011D3`.

Esos números salieron del osciloscopio antes de que existiera el generador. Si
alguno de los dos módulos se desvía de la norma, esas pruebas se ponen rojas
aunque el resto siga cuadrando consigo mismo.

---

## Pruebas

```
pytest -v
```

Dieciocho pruebas. Verde significa que el código funciona.

| Grupo | Qué comprueba |
|---|---|
| Anclaje | Cabecera y los dos CRC iguales a los medidos en la captura real |
| Capa física | 10 Mbit/s, TSS de 11 y de 5 bits, FES, histéresis, medida del tiempo de bit frente a un hueco anómalo |
| Robustez | 4 muestras por bit, ruido aumentado, amplitud reducida |
| Ida y vuelta | Identificadores 1, 42 y 2047; ciclos 0, 7 y 63; cargas de 2, 16 y 254 bytes |
| Canal | El canal se deduce de la semilla del CRC que cuadra |
| Detección de fallos | Un bit corrompido en la carga rompe el CRC de trama pero no el de cabecera; uno en la cabecera rompe el suyo |

La prueba de 4 muestras por bit es la que más dice: a 40 MS/s cada bit tiene
cuatro puntos, y el único motivo por el que la trama sigue saliendo entera es la
resincronización en cada BSS.

---

## Instalación

```
pip install -r requirements.txt
```

Python 3.10 o superior. numpy, pandas, matplotlib y pytest.

---

## Relacionado

[can-physical-decoder](https://github.com/andriusvillena-cmd/can-physical-decoder) —
lo mismo para CAN: de la señal de osciloscopio al mensaje, con verificación de
CRC y detección de tramas sin acuse de recibo.

## Normas de referencia

ISO 17458-2 (capa de enlace) · ISO 17458-4 (capa física)
